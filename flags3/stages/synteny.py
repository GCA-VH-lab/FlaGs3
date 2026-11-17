import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Optional

from flags3 import fasta
from flags3.genome import GenomeFasta
from flags3.log import debug, note, record_command
from flags3.schema import GenomeFiles, Row, Window
from flags3.stage import Stage
from flags3.stages.tree import Leaf
from flags3.tools import Tool, Tools

COLUMNS = ("qstart", "qend", "sstart", "send", "pident", "length")


class SyntenyError(RuntimeError):
	pass


@dataclass(frozen=True)
class Link(Row):
	FILE: ClassVar[str] = "links.tsv"
	row_a: str
	start_a: int
	end_a: int
	row_b: str
	start_b: int
	end_b: int
	strand: str
	identity: float
	length: int


class Pair:
	def __init__(self, upper: Window, lower: Window, query: Path, subject: Path, table: Path):
		self.upper = upper
		self.lower = lower
		self.query = query
		self.subject = subject
		self.table = table

	def links(self, min_identity: float, min_length: int) -> list[Link]:
		found = []
		for line in self.table.read_text().splitlines():
			hit = dict(zip(COLUMNS, line.split("\t")))
			identity, length = float(hit["pident"]), int(hit["length"])
			if identity < min_identity * 100.0 or length < min_length:
				continue
			sstart, send = int(hit["sstart"]), int(hit["send"])
			found.append(Link(
				self.upper.row_id, self.upper.lo + int(hit["qstart"]) - 1, self.upper.lo + int(hit["qend"]) - 1,
				self.lower.row_id, self.lower.lo + min(sstart, send) - 1, self.lower.lo + max(sstart, send) - 1,
				"+" if sstart <= send else "-", identity, length))
		return found


class Synteny(Stage):
	name = "synteny"
	requires = ("extract", "tree")
	needs = ("genome",)
	optional = True
	tool_name = "synteny"

	def wanted(self, config) -> bool:
		return config.flag("synteny")

	def run(self, run, config, out: Path) -> None:
		tool = Tools.load(config.tools_table()).tools[self.tool_name]
		windows = {w.row_id: w for w in Window.read(run.stage_file("extract", Window.FILE))}
		order = [l.row_id for l in Leaf.read(run.stage_file("tree", Leaf.FILE)) if l.row_id in windows]
		raw = out / "raw"
		raw.mkdir(parents=True, exist_ok=True)
		cut = self._cut(run, [windows[r] for r in order], raw)
		pairs = [Pair(windows[a], windows[b], cut[a], cut[b], raw / "{}.tsv".format(i + 1))
			for i, (a, b) in enumerate(zip(order, order[1:])) if a in cut and b in cut]
		note("comparing {} pairs of neighbouring rows".format(len(pairs)))
		evalue = tool.options.get("evalue", "1e-5")
		with ThreadPoolExecutor(max_workers=config.workers()) as pool:
			done = list(pool.map(lambda p: self._compare(tool, p, evalue), pairs))
		for argv, result in done:
			record_command(argv, result.returncode, stderr=result.stderr)
			if result.returncode != 0:
				raise SyntenyError("{} exited {}: {}".format(argv[0], result.returncode, (result.stderr or "").strip()[-300:]))
		identity = float(tool.options.get("identity", "0.6"))
		length = int(tool.options.get("length", "100"))
		links = [link for pair in pairs for link in pair.links(identity, length)]
		Link.write(out / Link.FILE, links)
		skipped = len(order) - 1 - len(pairs)
		note("{} synteny blocks between {} pairs{}".format(len(links), len(pairs),
			", {} pairs without a genome sequence".format(skipped) if skipped > 0 else ""))

	def _cut(self, run, windows: list[Window], raw: Path) -> dict[str, Path]:
		files = {g.assembly: g.path("genome") for g in GenomeFiles.read(run.stage_file("fetch", GenomeFiles.FILE))}
		by_assembly: dict[str, list[tuple[int, Window]]] = {}
		for n, window in enumerate(windows, 1):
			by_assembly.setdefault(window.assembly, []).append((n, window))
		cut = {}
		for assembly, members in by_assembly.items():
			path: Optional[Path] = files.get(assembly)
			if path is None or not path.is_file():
				continue
			contigs = GenomeFasta(path).contigs
			for n, window in members:
				sequence = contigs.get(window.contig, "")[window.lo - 1:window.hi]
				if sequence:
					cut[window.row_id] = raw / "{}.fna".format(n)
					fasta.write(cut[window.row_id], [(str(n), sequence)])
		return cut

	@staticmethod
	def _compare(tool: Tool, pair: Pair, evalue: str) -> tuple[list[str], subprocess.CompletedProcess]:
		argv = tool.argv(**{"in": pair.query, "subject": pair.subject, "out": pair.table, "evalue": evalue})
		debug("running: " + " ".join(argv))
		try:
			return argv, subprocess.run(argv, capture_output=True, text=True, cwd=tool.directory or None, env=tool.environment())
		except OSError as error:
			raise SyntenyError("could not run {}: {}".format(argv[0], error))
