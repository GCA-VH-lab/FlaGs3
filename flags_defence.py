import csv
import os
import re
import shutil
import subprocess
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple


class DefenceHit(NamedTuple):
	assembly: str
	contig: str
	start: int
	end: int
	type: str            # the system name, which is what the band is labelled with
	probability: float
	tools: str
	genes: Tuple[str, ...]


class Replicon(NamedTuple):
	name: str            # synthetic contig id, safe for filenames and tool ids
	row: str
	assembly: str
	contig: str
	genes: Tuple[dict, ...]


def build_replicons(neighborhoods) -> List[Replicon]:
	by_row: Dict[str, List] = {}
	for g in neighborhoods:
		if g.is_rna:
			continue
		by_row.setdefault(g.query, []).append(g)
	out = []
	for index, row in enumerate(sorted(by_row)):
		genes = sorted(by_row[row], key=lambda g: g.start)
		if not genes:
			continue
		out.append(Replicon(
			name="r{}".format(index), row=row,
			assembly=row.rsplit("|", 1)[-1], contig=genes[0].contig,
			genes=tuple({"accession": g.accession, "start": g.start, "end": g.end,
						 "strand": g.strand, "product": g.product} for g in genes)))
	return out


def write_inputs(replicons: List[Replicon], sequences: Dict[str, str],
				 gff_path: str, faa_path: str) -> Dict[str, Tuple[str, dict]]:
	index: Dict[str, Tuple[str, dict]] = {}
	with open(gff_path, "w") as gff, open(faa_path, "w") as faa:
		gff.write("##gff-version 3\n")
		for rep in replicons:
			lo = min(g["start"] for g in rep.genes)
			hi = max(g["end"] for g in rep.genes)
			gff.write("##sequence-region {} {} {}\n".format(rep.name, lo, hi))
		for rep in replicons:
			for order, gene in enumerate(rep.genes, start=1):
				tag = "{}_{}".format(rep.name, order)
				seq = sequences.get(gene["accession"])
				if not seq:
					continue
				gff.write("{}\tFlaGs3\tCDS\t{}\t{}\t.\t{}\t0\t"
						  "ID={};locus_tag={};product={}\n".format(
							  rep.name, gene["start"], gene["end"],
							  "-" if str(gene["strand"]).startswith("-") else "+",
							  tag, tag, gene["product"] or "hypothetical protein"))
				faa.write(">{}\n{}\n".format(tag, seq))
				index[tag] = (rep.name, gene)
	return index


def _span(tags: Iterable[str], index) -> Optional[Tuple[str, str, int, int, tuple]]:
	genes, replicon = [], None
	for tag in tags:
		entry = index.get(tag)
		if entry:
			replicon = entry[0]
			genes.append(entry[1])
	if not genes:
		return None
	return (replicon, min(g["start"] for g in genes), max(g["end"] for g in genes),
			tuple(g["accession"] for g in genes))


class DefenceScanner:
	def __init__(self, out_dir: str, threads: int = 0):
		self.out_dir = out_dir
		self.threads = threads
		os.makedirs(out_dir, exist_ok=True)
		self.statuses: Dict[str, str] = {}

	@staticmethod
	def available(tool: str) -> Tuple[bool, str]:
		import flags_tools
		cmd, _ = flags_tools.command(
			tool, **{"in": "x", "gff": "g", "faa": "f", "out": "y", "threads": 1})
		if not cmd:
			return False, "no {} row in tools_table.tsv".format(tool)
		return flags_tools.locate(cmd)

	def run(self, tool: str, replicons: List[Replicon],
			sequences: Dict[str, str]) -> List[DefenceHit]:
		work = os.path.join(self.out_dir, tool)
		os.makedirs(work, exist_ok=True)
		gff = os.path.join(work, "neighbourhoods.gff")
		faa = os.path.join(work, "neighbourhoods.faa")
		index = write_inputs(replicons, sequences, gff, faa)
		if not index:
			self.statuses[tool] = "skipped: no protein sequences to scan"
			return []
		out_dir = os.path.join(work, "out")
		shutil.rmtree(out_dir, ignore_errors=True)
		os.makedirs(out_dir, exist_ok=True)

		import flags_tools
		cmd, wd = flags_tools.command(
			tool, **{"in": faa, "faa": faa, "gff": gff, "out": out_dir,
					 "threads": self.threads or 1})
		try:
			proc = subprocess.run(cmd, cwd=wd or None, capture_output=True, text=True,
								  env=flags_tools.env_for(cmd))
			try:
				import flags_log
				flags_log.record_command(cmd, proc.returncode, proc.stdout, proc.stderr)
			except ImportError:
				pass
		except FileNotFoundError:
			self.statuses[tool] = "error: {} not found".format(cmd[0])
			return []
		if proc.returncode != 0:
			self.statuses[tool] = "error: exited {} ({})".format(
				proc.returncode, flags_tools.brief(proc.stderr or proc.stdout))
			return []

		reader = self._read_defensefinder if tool == "defensefinder" else self._read_padloc
		by_replicon = {r.name: r for r in replicons}
		hits = reader(out_dir, index, by_replicon)
		self.statuses[tool] = ("{} system(s) predicted".format(len(hits)) if hits
							   else "no defence system predicted")
		return hits

	@staticmethod
	def _find(out_dir: str, *endings) -> List[str]:
		found = []
		for base, _, names in os.walk(out_dir):
			for name in names:
				if name.endswith(endings):
					found.append(os.path.join(base, name))
		return sorted(found)

	@staticmethod
	def _emit(replicon, lo, hi, accessions, label, score, tool, by_replicon):
		rep = by_replicon.get(replicon)
		if rep is None:
			return None
		return DefenceHit(assembly=rep.assembly, contig=rep.contig,
						  start=lo, end=hi, type=label,
						  probability=score, tools=tool, genes=accessions)

	def _read_defensefinder(self, out_dir, index, by_replicon) -> List[DefenceHit]:
		hits = []
		for path in self._find(out_dir, "_systems.tsv", "systems.tsv"):
			with open(path, newline="") as fh:
				for row in csv.DictReader(fh, delimiter="\t"):
					tags = re.split(r"[,\s]+", (row.get("protein_in_syst") or "").strip())
					placed = _span([t for t in tags if t], index)
					if not placed:
						continue
					replicon, lo, hi, accessions = placed
					label = (row.get("subtype") or row.get("type") or "defence system")
					hit = self._emit(replicon, lo, hi, accessions, label, 0.0,
									 "defensefinder", by_replicon)
					if hit:
						hits.append(hit)
		return hits

	def _read_padloc(self, out_dir, index, by_replicon) -> List[DefenceHit]:
		groups: Dict[Tuple[str, str, str], List[str]] = {}
		labels: Dict[Tuple[str, str, str], str] = {}
		for path in self._find(out_dir, "_padloc.csv", "padloc.csv"):
			with open(path, newline="") as fh:
				for row in csv.DictReader(fh):
					tag = (row.get("target.name") or "").strip()
					system = (row.get("system") or "").strip()
					number = (row.get("system.number") or "").strip()
					seqid = (row.get("seqid") or "").strip()
					if not tag or not system:
						continue
					key = (seqid, system, number)
					groups.setdefault(key, []).append(tag)
					labels[key] = system
		hits = []
		for key, tags in groups.items():
			placed = _span(tags, index)
			if not placed:
				continue
			replicon, lo, hi, accessions = placed
			hit = self._emit(replicon, lo, hi, accessions, labels[key], 0.0,
							 "padloc", by_replicon)
			if hit:
				hits.append(hit)
		return hits


def merge_calls(hits: List[DefenceHit]) -> List[DefenceHit]:
	seen: Dict[Tuple[str, str, int, int, str], DefenceHit] = {}
	for h in hits:
		key = (h.assembly, h.contig, h.start, h.end, h.type)
		if key in seen:
			tools = sorted(set(seen[key].tools.split(",")) | {h.tools})
			seen[key] = seen[key]._replace(tools=",".join(tools))
		else:
			seen[key] = h
	return sorted(seen.values(),
				  key=lambda h: (h.assembly, h.contig, h.start, h.type))


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


def write_report(hits: List[DefenceHit], matches: Dict[int, List[str]], path: str):
	with open(path, "w") as out:
		out.write("#assembly\tcontig\tstart\tend\ttype\tprobability\tcalled_by\t"
				  "n_genes\tgenes\toverlapping_rows\n")
		for i, h in enumerate(hits):
			out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
				h.assembly, h.contig, h.start, h.end, h.type, h.probability,
				h.tools, len(h.genes), ",".join(h.genes),
				",".join(matches.get(i, [])) or "-"))


def write_diagnostics(statuses: Dict[str, str], replicons: int, path: str):
	with open(path, "w") as out:
		out.write("# {} neighbourhood(s) written as synthetic replicons\n".format(
			replicons))
		out.write("#tool\tstatus\n")
		for tool in sorted(statuses):
			out.write("{}\t{}\n".format(tool, " ".join(statuses[tool].split())))
