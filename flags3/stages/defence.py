import csv
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from flags3 import fasta, windows
from flags3.log import debug, note, record_command
from flags3.schema import MISSING, Annotation, Gene, Row, Window
from flags3.stage import Stage
from flags3.tools import Tool, Tools, brief

TOOLS = ("defensefinder", "padloc")


@dataclass(frozen=True)
class DefenceCall(Row):
	FILE: ClassVar[str] = "defence.tsv"
	assembly: str
	contig: str
	start: int
	end: int
	type: str
	called_by: str
	genes: str
	rows: str


@dataclass(frozen=True)
class ToolStatus(Row):
	FILE: ClassVar[str] = "diagnostics.tsv"
	tool: str
	systems: int
	status: str


class Replicons:
	def __init__(self, genes: list[Gene], sequences: dict[str, str]):
		self.rows: dict[str, list[Gene]] = {}
		for g in genes:
			if not g.is_rna and g.accession in sequences:
				self.rows.setdefault(g.row_id, []).append(g)
		self.names = {row: "r{}".format(i) for i, row in enumerate(sorted(self.rows))}
		self.sequences = sequences
		self.tags: dict[str, Gene] = {}

	def write(self, gff: Path, faa: Path) -> int:
		with open(gff, "w") as g_out, open(faa, "w") as f_out:
			g_out.write("##gff-version 3\n")
			for row, genes in sorted(self.rows.items()):
				g_out.write("##sequence-region {} {} {}\n".format(self.names[row], min(g.start for g in genes), max(g.end for g in genes)))
			for row, genes in sorted(self.rows.items()):
				for order, gene in enumerate(sorted(genes, key=lambda g: g.start), start=1):
					tag = "{}_{}".format(self.names[row], order)
					g_out.write("{}\tFlaGs3\tCDS\t{}\t{}\t.\t{}\t0\tID={};locus_tag={};product={}\n".format(
						self.names[row], gene.start, gene.end, gene.strand, tag, tag,
						gene.product if gene.product != MISSING else "hypothetical protein"))
					f_out.write(">{}\n{}\n".format(tag, self.sequences[gene.accession]))
					self.tags[tag] = gene
		return len(self.tags)

	def span(self, tags) -> tuple[Gene, ...]:
		return tuple(self.tags[t] for t in tags if t in self.tags)


def read_defensefinder(out_dir: Path):
	for path in sorted(p for p in out_dir.rglob("*systems.tsv")):
		with open(path, newline="") as handle:
			for row in csv.DictReader(handle, delimiter="\t"):
				tags = [t for t in re.split(r"[,\s]+", (row.get("protein_in_syst") or "").strip()) if t]
				yield (row.get("subtype") or row.get("type") or "defence system"), tags


def read_padloc(out_dir: Path):
	groups: dict[tuple, list[str]] = {}
	for path in sorted(p for p in out_dir.rglob("*padloc.csv")):
		with open(path, newline="") as handle:
			for row in csv.DictReader(handle):
				tag, system = (row.get("target.name") or "").strip(), (row.get("system") or "").strip()
				if tag and system:
					groups.setdefault(((row.get("seqid") or "").strip(), system, (row.get("system.number") or "").strip()), []).append(tag)
	for (_, system, _), tags in groups.items():
		yield system, tags


READERS = {"defensefinder": read_defensefinder, "padloc": read_padloc}


class Defence(Stage):
	name = "defence"
	requires = ("extract",)
	optional = True

	def wanted(self, config) -> bool:
		return any(config.flag(t) for t in TOOLS)

	def run(self, run, config, out: Path) -> None:
		tools = Tools.load(config.path("tools"))
		genes = Gene.read(run.stage_file("extract", Gene.FILE))
		rows = Window.read(run.stage_file("extract", Window.FILE))
		sequences = {n: s for n, _, s in fasta.read(run.stage_file("extract", "proteins.faa"))}
		replicons = Replicons(genes, sequences)
		raw = out / "raw"
		raw.mkdir(parents=True, exist_ok=True)
		count = replicons.write(raw / "neighbourhoods.gff", raw / "neighbourhoods.faa")
		note("defence: {} proteins in {} synthetic replicons".format(count, len(replicons.rows)))
		calls: dict[tuple, dict] = {}
		statuses = []
		for name in (t for t in TOOLS if config.flag(t)):
			tool = tools[name]
			found, where = tool.locate()
			if not found:
				statuses.append(ToolStatus(name, 0, "not run: {}; flags3 install {} sets it up".format(where, name)))
				print("Warning: {} not run ({}).".format(name, where))
				continue
			try:
				systems = self._call(tool, name, raw, replicons, config)
			except RuntimeError as error:
				statuses.append(ToolStatus(name, 0, "failed: {}".format(error)))
				print("Warning: {} failed ({}).".format(name, error))
				continue
			for label, members in systems:
				key = (members[0].assembly, members[0].contig, min(g.start for g in members), max(g.end for g in members), label)
				entry = calls.setdefault(key, {"tools": set(), "genes": tuple(g.accession for g in members)})
				entry["tools"].add(name)
			statuses.append(ToolStatus(name, len(systems), "ok"))
		ToolStatus.write(out / ToolStatus.FILE, statuses)
		placed = []
		for (assembly, contig, start, end, label), entry in sorted(calls.items()):
			overlap = windows.overlapping_rows(rows, assembly, contig, start, end)
			placed.append(DefenceCall(assembly, contig, start, end, label, ",".join(sorted(entry["tools"])),
				",".join(entry["genes"]), ",".join(overlap) or MISSING))
		DefenceCall.write(out / DefenceCall.FILE, placed)
		Annotation.write(out / Annotation.FILE, (
			Annotation("{}|{}".format(c.assembly, c.contig), "bp", c.start, c.end, "band", c.type, c.type, c.called_by, None)
			for c in placed))
		note("defence: {} systems".format(len(placed)))
		if statuses and all(s.status != "ok" for s in statuses):
			raise RuntimeError("no defence tool ran; see defence/{}".format(ToolStatus.FILE))

	def _call(self, tool: Tool, name: str, raw: Path, replicons: Replicons, config) -> list:
		out_dir = raw / name
		out_dir.mkdir(exist_ok=True)
		argv = tool.argv(**{"in": raw / "neighbourhoods.faa", "faa": raw / "neighbourhoods.faa",
			"gff": raw / "neighbourhoods.gff", "out": out_dir, "threads": config.workers()})
		debug("running: " + " ".join(argv))
		done = subprocess.run(argv, capture_output=True, text=True, cwd=tool.directory or None, env=tool.environment())
		record_command(argv, done.returncode, done.stdout, done.stderr)
		if done.returncode != 0:
			raise RuntimeError("exited {}: {}".format(done.returncode, brief(done.stderr or done.stdout)))
		systems = []
		for label, tags in READERS[name](out_dir):
			members = replicons.span(tags)
			if members:
				systems.append((label, members))
		return systems
