import gzip
import os
import shutil
import subprocess
from typing import Dict, List, NamedTuple, Optional, Tuple

from flags_scan import ScanWindow, merge_windows, place, scanned_bases, write_windows


class SecretionHit(NamedTuple):
	assembly: str
	contig: str
	start: int          
	end: int
	type: str           
	probability: float  
	columns: Tuple[str, ...]   
	values: Tuple[str, ...]    
	coverage: str = "whole-genome"   # full | partial | whole-genome


class SismisScanner: 
	CONTIG_COLS = ("sequence_id", "seq_id", "contig", "scaffold", "record_id")

	def __init__(self, out_dir: str):
		self.out_dir = out_dir
		os.makedirs(out_dir, exist_ok=True)
		self._cache: Dict[str, List[SecretionHit]] = {}   
		self.scanned_bases: Dict[str, int] = {}
		self.window_count: Dict[str, int] = {}
		self.windows_dir = os.path.join(out_dir, os.pardir, "windows")
		self.window_key = "genome"
		self.batches = 0

	def scan(self, jobs, statuses):
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
			asm_dir = os.path.join(self.out_dir, "batch{:03d}".format(batch))
			os.makedirs(asm_dir, exist_ok=True)
			try:
				raw = self._run("batch {}".format(batch), fasta, asm_dir)
			except Exception as e:
				for assembly in here:
					statuses[assembly] = "error: {}".format(e)
				continue
			counted = {}
			for hit in raw:
				placed = flags_scan.place_batched(
					hit.contig, hit.start, hit.end, offsets)
				if placed is None:
					continue
				assembly, contig, start, end, coverage = placed
				hits.append(hit._replace(assembly=assembly, contig=contig,
										 start=start, end=end, coverage=coverage))
				counted[assembly] = counted.get(assembly, 0) + 1
			for assembly in here:
				statuses[assembly] = "{} window(s) scanned, {}".format(
					self.window_count.get(assembly, 0),
					"{} secretion system(s) predicted".format(counted[assembly])
					if counted.get(assembly) else "no secretion system predicted")
		self.batches = batch
		return hits

	def scan_windows(self, assembly: str, genome_path: str,
					 windows: List[ScanWindow]) -> List[SecretionHit]:
		if assembly in self._cache:
			return self._cache[assembly]
		if not windows:
			return []
		asm_dir = os.path.join(self.out_dir, assembly)
		os.makedirs(asm_dir, exist_ok=True)
		fasta = os.path.join(asm_dir, "windows.fna")
		offsets = write_windows(genome_path, windows, fasta)
		if not offsets:
			raise RuntimeError(
				"none of the {} scan window(s) matched a contig in {}".format(
					len(windows), os.path.basename(genome_path)))
		self.scanned_bases[assembly] = scanned_bases(offsets)
		self.window_count[assembly] = len(offsets)
		hits = self._run(assembly, fasta, asm_dir)
		mapped = []
		for hit in hits:
			placed = place(hit.contig, hit.start, hit.end, offsets)
			if placed is None:
				continue
			contig, start, end, coverage = placed
			mapped.append(hit._replace(contig=contig, start=start, end=end,
									   coverage=coverage))
		self._cache[assembly] = mapped
		return mapped

	def scan_assembly(self, assembly: str, genome_path: str) -> List[SecretionHit]:
		if assembly in self._cache:
			return self._cache[assembly]

		asm_dir = os.path.join(self.out_dir, assembly)
		os.makedirs(asm_dir, exist_ok=True)
		fasta = self._decompressed(genome_path, asm_dir)

		hits = self._run(assembly, fasta, asm_dir)
		self._cache[assembly] = hits
		return hits

	def _run(self, assembly: str, fasta: str, asm_dir: str) -> List[SecretionHit]:
		sismis_out = os.path.join(asm_dir, "sismis")
		import flags_tools
		cmd, wd = flags_tools.command("sismis", **{"in": fasta, "out": sismis_out})
		try:
			proc = subprocess.run(cmd, cwd=wd or None, capture_output=True, text=True,
								  env=flags_tools.env_for(cmd))
			try:
				import flags_log
				flags_log.record_command(cmd, proc.returncode, proc.stdout, proc.stderr)
			except ImportError:
				pass
		except FileNotFoundError:
			raise FileNotFoundError("sismis not found on PATH; install with 'pip install sismis'.")
		if proc.returncode != 0:
			raise RuntimeError("sismis exited {} on {}: {}".format(
				proc.returncode, assembly,
				flags_tools.brief(proc.stderr or proc.stdout)))
		return self._parse_clusters(sismis_out, assembly)

	@staticmethod
	def _decompressed(genome_path: str, asm_dir: str) -> str:
		if not genome_path.endswith(".gz"):
			return genome_path
		local = os.path.join(asm_dir, "genome.fna")
		with gzip.open(genome_path, "rt", encoding="utf-8", errors="replace") as fin, \
			 open(local, "w") as fout:
			shutil.copyfileobj(fin, fout)
		return local

	def _parse_clusters(self, sismis_out: str, assembly: str) -> List[SecretionHit]:
		clusters_tsv = self._find_clusters_tsv(sismis_out)
		if clusters_tsv is None:
			return []
		hits = []
		with open(clusters_tsv) as fh:
			header = fh.readline().rstrip("\n").split("\t")
			lower = [h.lower() for h in header]
			contig_i = next((lower.index(c) for c in self.CONTIG_COLS if c in lower), None)
			start_i = lower.index("start") if "start" in lower else None
			end_i = lower.index("end") if "end" in lower else None
			type_i = lower.index("type") if "type" in lower else None
			prob_i = lower.index("max_p") if "max_p" in lower else None
			missing = [name for name, i in (("contig", contig_i), ("start", start_i),
											 ("end", end_i), ("type", type_i),
											 ("max_p", prob_i)) if i is None]
			if missing:
				raise RuntimeError(
					"{} is missing column(s) {} needed to use its output; "
					"columns seen: {}".format(clusters_tsv, missing, header))
			for line in fh:
				if not line.strip():
					continue
				values = line.rstrip("\n").split("\t")
				hits.append(SecretionHit(
					assembly=assembly, contig=values[contig_i],
					start=int(values[start_i]), end=int(values[end_i]),
					type=values[type_i], probability=float(values[prob_i]),
					columns=tuple(header), values=tuple(values)))
		return hits

	@staticmethod
	def _find_clusters_tsv(sismis_out: str) -> Optional[str]:
		if not os.path.isdir(sismis_out):
			return None
		for name in os.listdir(sismis_out):
			if name.endswith(".clusters.tsv"):
				return os.path.join(sismis_out, name)
		return None


def match_rows(hits: List[SecretionHit], rows: Dict[str, Tuple[str, str, int, int]]
			   ) -> Dict[int, List[str]]:
	index: Dict[tuple, list] = {}
	for row_id, (assembly, contig, lo, hi) in rows.items():
		index.setdefault((assembly, contig), []).append((lo, hi, row_id))
	matches: Dict[int, List[str]] = {}
	for i, h in enumerate(hits):
		for lo, hi, row_id in index.get((h.assembly, h.contig), ()):
			if h.start <= hi and h.end >= lo:
				matches.setdefault(i, []).append(row_id)
	return matches


def write_report(hits: List[SecretionHit], matches: Dict[int, List[str]], path: str):
	with open(path, "w") as out:
		if not hits:
			out.write("#assembly\t(no secretion systems predicted)\n")
			return
		header = hits[0].columns
		out.write("#assembly\tcontig\tstart\tend\ttype\tprobability\tcoverage\t{}"
				  "\toverlapping_rows\n".format(
					  "\t".join("sismis_" + h for h in header)))
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
			fraction = ("{:.4f}".format(bases / total)
						if bases and total else "-")
			out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
				assembly, " ".join(statuses[assembly].split()),
				windows.get(assembly, "-"),
				bases if bases is not None else "-",
				total if total is not None else "-", fraction))