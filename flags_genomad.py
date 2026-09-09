import csv
import gzip
import os
import shutil
import subprocess
from typing import Dict, List, NamedTuple, Optional, Tuple

from flags_scan import ScanWindow, merge_windows, place, scanned_bases, write_windows


class MgeHit(NamedTuple):
	assembly: str
	contig: str
	start: int
	end: int
	type: str            # the classification itself, not a generic label
	probability: float
	columns: Tuple[str, ...]
	values: Tuple[str, ...]
	coverage: str = "whole-genome"


DB_MARKERS = ("genomad_marker_metadata.tsv", "version.txt")


def resolve_database(path: str):
	if not path:
		return "", "no geNomad database. Run genomad_loader.sh, which downloads it "\
				   "and writes its path into the third column of the genomad row "\
				   "of tools_table.tsv"
	if not os.path.isdir(path):
		return "", "geNomad database not found at {}".format(path)
	if any(os.path.isfile(os.path.join(path, m)) for m in DB_MARKERS):
		return path, ""
	nested = os.path.join(path, "genomad_db")
	if any(os.path.isfile(os.path.join(nested, m)) for m in DB_MARKERS):
		return nested, "nested"
	return "", ("{} does not look like a geNomad database; it has none of {}. "
				"geNomad creates a genomad_db directory inside whatever "
				"destination it is given, so the path wanted here is usually one "
				"level further in".format(path, " or ".join(DB_MARKERS)))


def virus_name(taxonomy: str) -> str:
	ranks = [r.strip() for r in (taxonomy or "").split(";") if r.strip()]
	ranks = [r for r in ranks if r.lower() not in ("viruses", "unclassified")]
	return ranks[-1] if ranks else "virus"


def plasmid_name(conjugation: str, amr: str) -> str:
	parts = []
	if (conjugation or "").strip() not in ("", "0", "NA"):
		parts.append("conjugative")
	parts.append("plasmid")
	if (amr or "").strip() not in ("", "0", "NA"):
		parts.append("(AMR)")
	return " ".join(parts)


class GenomadScanner:
	def __init__(self, out_dir: str, database: str = "", threads: int = 0,
				 windows_dir: str = "", window_key: str = "genome"):
		self.out_dir = out_dir
		self.windows_dir = windows_dir or os.path.join(out_dir, "windows")
		self.window_key = window_key
		self.batches = 0
		self.database = database
		self.threads = threads
		os.makedirs(out_dir, exist_ok=True)
		self._cache: Dict[str, List[MgeHit]] = {}
		self.scanned_bases: Dict[str, int] = {}
		self.window_count: Dict[str, int] = {}

	def available(self) -> Tuple[bool, str]:
		import flags_tools
		cmd, wd = flags_tools.command(
			"genomad", **{"in": "x", "out": "y", "db": self.database or "db",
						  "threads": self.threads or 1})
		if not cmd:
			return False, "no genomad row in tools_table.tsv"
		ok, where = flags_tools.locate(cmd)
		if not ok:
			return False, where
		resolved, note = resolve_database(self.database)
		if not resolved:
			return False, note
		if note == "nested":
			print("Note: using the geNomad database at {}, one level inside the "
				  "path in tools_table.tsv.".format(resolved))
		self.database = resolved
		return True, where

	def scan(self, jobs, statuses):
		"""jobs: [(assembly, genome_path, windows)]. geNomad's cost is dominated
		by loading its database and model, so every assembly goes into as few
		invocations as possible rather than one each."""
		import flags_scan
		hits = []
		batch = 0
		for fasta, offsets in flags_scan.shared_batches(
				self.window_key, jobs, self.windows_dir):
			batch += 1
			here = {a for a, _, _ in offsets.values()}
			for assembly in here:
				mine = [w for a, _, w in offsets.values() if a == assembly]
				self.scanned_bases[assembly] = sum(
					w.slice_end - w.slice_start + 1 for w in mine)
				self.window_count[assembly] = len(mine)
			try:
				found = self._run("batch {}".format(batch), fasta,
								  os.path.join(self.out_dir, "batch{:03d}".format(batch)),
								  offsets)
			except Exception as e:
				for assembly in here:
					statuses[assembly] = "error: {}".format(e)
				continue
			hits.extend(found)
			counted = {}
			for h in found:
				counted[h.assembly] = counted.get(h.assembly, 0) + 1
			for assembly in here:
				statuses[assembly] = "{}, {}".format(
					"{} window(s) scanned".format(self.window_count.get(assembly, 0)),
					"{} mobile element(s) predicted".format(counted.get(assembly, 0))
					if counted.get(assembly) else "no mobile element predicted")
		self.batches = batch
		return hits

	def _run(self, assembly: str, fasta: str, asm_dir: str, offsets) -> List[MgeHit]:
		os.makedirs(asm_dir, exist_ok=True)
		out_dir = os.path.join(asm_dir, "genomad")
		import flags_tools
		values = {"in": fasta, "out": out_dir, "db": self.database,
				  "threads": self.threads or 1}
		empty = flags_tools.missing_values("genomad", values)
		if empty:
			raise RuntimeError(
				"the genomad command in tools_table.tsv needs {} but it is empty; "
				"the tool would be run without it".format(
					" and ".join("{" + e + "}" for e in empty)))
		cmd, wd = flags_tools.command("genomad", **values)
		try:
			proc = subprocess.run(cmd, cwd=wd or None, capture_output=True, text=True,
								  env=flags_tools.env_for(cmd))
			try:
				import flags_log
				flags_log.record_command(cmd, proc.returncode, proc.stdout, proc.stderr)
			except ImportError:
				pass
		except FileNotFoundError:
			raise FileNotFoundError(
				"genomad not found; run genomad_installer.sh.")
		if proc.returncode != 0:
			raise RuntimeError("genomad exited {} on {}: {}".format(
				proc.returncode, assembly,
				flags_tools.brief(proc.stderr or proc.stdout)))
		return self._collect(assembly, out_dir, offsets)

	def _collect(self, assembly: str, out_dir: str, offsets) -> List[MgeHit]:  # noqa
		hits = []
		for path, kind in self._summaries(out_dir):
			hits.extend(self._parse(path, kind, assembly, offsets))
		return hits

	@staticmethod
	def _summaries(out_dir: str):
		found = []
		for base, _, names in os.walk(out_dir):
			for name in names:
				if name.endswith("_virus_summary.tsv"):
					found.append((os.path.join(base, name), "virus"))
				elif name.endswith("_plasmid_summary.tsv"):
					found.append((os.path.join(base, name), "plasmid"))
		return sorted(found)

	def _parse(self, path: str, kind: str, assembly: str, offsets) -> List[MgeHit]:
		hits = []
		with open(path, newline="") as fh:
			reader = csv.DictReader(fh, delimiter="\t")
			header = tuple(reader.fieldnames or ())
			for row in reader:
				record, start, end = self._span(row)
				if record is None:
					continue
				if kind == "virus":
					label = virus_name(row.get("taxonomy", ""))
					score = row.get("virus_score")
				else:
					label = plasmid_name(row.get("conjugation_genes", ""),
										 row.get("amr_genes", ""))
					score = row.get("plasmid_score")
				import flags_scan
				placed = flags_scan.place_batched(record, start, end, offsets)
				if placed is None:
					continue
				assembly, contig, start, end, coverage = placed
				hits.append(MgeHit(
					assembly=assembly, contig=contig, start=start, end=end,
					type=label, probability=float(score or 0.0),
					columns=header,
					values=tuple(row.get(c, "") or "" for c in header),
					coverage=coverage))
		return hits

	@staticmethod
	def _span(row):
		name = (row.get("seq_name") or "").strip()
		if not name:
			return None, 0, 0
		record = name.split("|")[0]
		coords = (row.get("coordinates") or "").strip()
		if coords and coords.upper() != "NA" and "-" in coords:
			lo, _, hi = coords.partition("-")
			try:
				return record, int(lo), int(hi)
			except ValueError:
				pass
		try:
			return record, 1, int(row.get("length") or 0)
		except ValueError:
			return record, 1, 0


def match_rows(hits, rows) -> Dict[int, List[str]]:
	index: Dict[tuple, list] = {}
	for row_id, (assembly, contig, lo, hi) in rows.items():
		index.setdefault((assembly, contig), []).append((lo, hi, row_id))
	matches: Dict[int, List[str]] = {}
	for i, h in enumerate(hits):
		for lo, hi, row_id in index.get((h.assembly, h.contig), ()):
			if h.start <= hi and h.end >= lo:
				matches.setdefault(i, []).append(row_id)
	return matches


def write_report(hits: List[MgeHit], matches: Dict[int, List[str]], path: str):
	with open(path, "w") as out:
		if not hits:
			out.write("#assembly\t(no mobile elements predicted)\n")
			return
		header = hits[0].columns
		out.write("#assembly\tcontig\tstart\tend\ttype\tprobability\tcoverage\t{}"
				  "\toverlapping_rows\n".format(
					  "\t".join("genomad_" + h for h in header)))
		for i, h in enumerate(hits):
			values = h.values
			if len(values) != len(header):
				values = (values + ("",) * len(header))[:len(header)]
			rows = ",".join(matches.get(i, [])) or "-"
			out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
				h.assembly, h.contig, h.start, h.end, h.type, h.probability,
				h.coverage, "\t".join(values), rows))


def write_diagnostics(statuses: Dict[str, str], path: str,
					  scanned: Optional[Dict[str, int]] = None,
					  windows: Optional[Dict[str, int]] = None,
					  genome_size: Optional[Dict[str, int]] = None):
	scanned = scanned or {}
	windows = windows or {}
	genome_size = genome_size or {}
	with open(path, "w") as out:
		out.write("#assembly\tstatus\twindows\tbases_scanned\tbases_in_genome\t"
				  "fraction_scanned\n")
		for assembly in sorted(statuses):
			bases = scanned.get(assembly)
			total = genome_size.get(assembly)
			fraction = "{:.4f}".format(bases / total) if bases and total else "-"
			out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
				assembly, " ".join(statuses[assembly].split()),
				windows.get(assembly, "-"),
				bases if bases is not None else "-",
				total if total is not None else "-", fraction))
