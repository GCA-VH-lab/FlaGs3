import csv
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Optional

from flags3 import windows
from flags3.log import debug, note, record_command
from flags3.schema import MISSING, Annotation, GenomeFiles, Row, Window
from flags3.stage import Stage
from flags3.tools import Tool, Tools, brief


class ScanError(RuntimeError):
	pass


@dataclass
class Hit:
	record: str
	start: int
	end: int
	type: str
	score: Optional[float]
	label: str = ""
	extra: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Placed:
	assembly: str
	contig: str
	start: int
	end: int
	type: str
	score: Optional[float]
	label: str
	coverage: str
	rows: str
	extra: dict


@dataclass(frozen=True)
class Diagnostic(Row):
	FILE: ClassVar[str] = "diagnostics.tsv"
	assembly: str
	windows: int
	bases_scanned: int
	hits: int
	status: str


class WindowScan(Stage):
	requires = ("extract",)
	needs = ("genome",)
	optional = True
	tool_name = ""
	report_file = ""
	passthrough_prefix = ""
	kind = "band"

	def run(self, run, config, out: Path) -> None:
		tools = Tools.load(config.path("tools"))
		tool = tools[self.tool_name]
		found, where = tool.locate()
		if not found:
			raise ScanError("{}: {}; flags3 install {} sets it up".format(self.tool_name, where, self.tool_name))
		rows = Window.read(run.stage_file("extract", Window.FILE))
		genomes = {g.assembly: g for g in GenomeFiles.read(run.stage_file("fetch", GenomeFiles.FILE))}
		usable, missing = {}, []
		for assembly in sorted({w.assembly for w in rows}):
			path = genomes[assembly].path("genome") if assembly in genomes else None
			if path is None:
				missing.append(assembly)
			else:
				usable[assembly] = path
		if missing:
			print("Warning: {} has no genome FASTA for {} assemblies; they are not scanned.".format(
				self.name, len(missing)))
		cuts = windows.merge([w for w in rows if w.assembly in usable])
		batches = windows.Batches(out / "raw")
		paths = batches.write(usable, cuts)
		note("{}: {} windows in {} batch files".format(self.tool_name, len(batches.cuts), len(paths)))
		by_record = batches.by_record()
		placed: list[Placed] = []
		failed: dict[str, str] = {}
		for path in paths:
			try:
				hits = self.scan_batch(tool, path, out / "raw" / (path.stem + "_" + self.tool_name), config)
			except ScanError as error:
				for c in batches.cuts:
					if c.batch == path.name:
						failed[c.assembly] = str(error)
				print("Warning: {} failed on {}: {}".format(self.tool_name, path.name, error))
				continue
			for hit in hits:
				cut = by_record.get(hit.record) or by_record.get(hit.record.split()[0])
				if cut is None:
					continue
				lo, hi, coverage = cut.place(hit.start, hit.end)
				placed.append(Placed(cut.assembly, cut.contig, lo, hi, hit.type, hit.score, hit.label or hit.type,
					coverage, ",".join(windows.overlapping_rows(rows, cut.assembly, cut.contig, lo, hi)) or MISSING, hit.extra))
		placed.sort(key=lambda p: (p.assembly, p.contig, p.start, p.end))
		self.write_report(out / self.report_file, placed)
		Annotation.write(out / Annotation.FILE, (
			Annotation("{}|{}".format(p.assembly, p.contig), "bp", p.start, p.end, self.kind, p.type, p.label,
				self.tool_name, p.score) for p in placed))
		counts = {}
		for p in placed:
			counts[p.assembly] = counts.get(p.assembly, 0) + 1
		stats = batches.per_assembly()
		Diagnostic.write(out / Diagnostic.FILE, (
			Diagnostic(a, stats.get(a, (0, 0))[0], stats.get(a, (0, 0))[1], counts.get(a, 0),
				failed.get(a) or ("not scanned: no genome FASTA" if a in missing else "ok"))
			for a in sorted(set(usable) | set(missing))))
		note("{}: {} hits on {} assemblies".format(self.tool_name, len(placed), len(counts)))
		if failed and len(failed) == len(usable):
			raise ScanError("{} failed on every batch".format(self.tool_name))

	def scan_batch(self, tool: Tool, fasta: Path, out_dir: Path, config) -> list[Hit]:
		raise NotImplementedError

	def execute(self, tool: Tool, argv: list[str]) -> subprocess.CompletedProcess:
		debug("running: " + " ".join(argv))
		try:
			done = subprocess.run(argv, capture_output=True, text=True, cwd=tool.directory or None, env=tool.environment())
		except OSError as error:
			raise ScanError("could not run {}: {}".format(argv[0], error))
		record_command(argv, done.returncode, done.stdout, done.stderr)
		if done.returncode != 0:
			raise ScanError("{} exited {}: {}".format(argv[0], done.returncode, brief(done.stderr or done.stdout)))
		return done

	def write_report(self, path: Path, placed: list[Placed]) -> None:
		extra_columns: list[str] = []
		for p in placed:
			for key in p.extra:
				if key not in extra_columns:
					extra_columns.append(key)
		with open(path, "w", newline="", encoding="utf-8") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["assembly", "contig", "start", "end", "type", "score", "coverage", "rows"]
				+ [self.passthrough_prefix + c for c in extra_columns])
			for p in placed:
				writer.writerow([p.assembly, p.contig, p.start, p.end, p.type,
					MISSING if p.score is None else format(p.score, ".4g"), p.coverage, p.rows]
					+ [p.extra.get(c, MISSING) for c in extra_columns])
