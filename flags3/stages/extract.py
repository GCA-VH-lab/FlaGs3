from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from flags3 import fasta
from flags3.genome import Feature, GeneTable, GenomeFasta, ProteinFasta, RnaFasta
from flags3.log import note
from flags3.schema import (MISSING, Gene, GenomeFiles, QueryTarget, RangeReport, RowInfo,
	Unmatched, Window, row_id)
from flags3.stage import Stage

SCAN_GENES = "scan_genes.tsv"
PROTEINS = "proteins.faa"
QUERIES = "queries.faa"
RNA = "rna.fna"
SCAN_PROTEINS = "scan_proteins.faa"


@dataclass(frozen=True)
class Settings:
	flank: int = 4
	range_bp: Optional[int] = None
	scan_range: Optional[int] = None
	margin: int = 10000

	@classmethod
	def from_config(cls, config) -> "Settings":
		return cls(flank=config.integer("gene", 4), range_bp=config.integer("range"),
			scan_range=config.integer("scan_range"), margin=config.integer("scan_margin", 10000))


class Neighbourhood:
	def __init__(self, table: GeneTable, index: int, settings: Settings):
		self.table = table
		self.index = index
		self.settings = settings
		self.query = table.genes[index]
		self.contig = table.contig_of(index)
		self.lo, self.hi = self._window(settings.range_bp)

	def _window(self, reach_bp: Optional[int]) -> tuple[int, int]:
		genes, contig, idx = self.table.genes, self.contig, self.index
		if not reach_bp:
			return max(contig.first, idx - self.settings.flank), min(contig.last, idx + self.settings.flank + 1)
		lo_bp = self.query.start - reach_bp
		hi_bp = self.query.end + reach_bp
		lo = idx
		while lo > contig.first and contig.reach[lo - 1 - contig.first] >= lo_bp:
			lo -= 1
		hi = idx + 1
		while hi < contig.last and genes[hi].start <= hi_bp:
			hi += 1
		return lo, hi

	def offset(self, index: int) -> int:
		delta = index - self.index
		return -delta if self.query.strand == "-" else delta

	def genes(self) -> list[tuple[int, Feature]]:
		return [(j, self.table.genes[j]) for j in range(self.lo, self.hi)]

	def scan_genes(self) -> list[tuple[int, Feature]]:
		if not self.settings.scan_range:
			return []
		lo, hi = self._window(self.settings.scan_range)
		return [(j, self.table.genes[j]) for j in range(lo, hi) if not self.table.genes[j].is_rna]

	def scan_span(self, margin: int = 0) -> tuple[int, int]:
		length = self.contig.length
		if not self.settings.scan_range:
			return 1, length
		reach = self.settings.scan_range + margin
		return max(1, self.query.start - reach), min(length, self.query.end + reach)

	def span(self) -> tuple[int, int]:
		members = [g for _, g in self.genes()]
		return min(g.start for g in members), max(g.end for g in members)

	def range_report(self, row: str) -> RangeReport:
		q, length = self.query, self.contig.length
		lo_bp, hi_bp = self.span()
		left_avail, right_avail = q.start - 1, max(0, length - q.end)
		left_reach, right_reach = q.start - lo_bp, hi_bp - q.end
		before = sum(1 for j, _ in self.genes() if self.offset(j) < 0)
		after = sum(1 for j, _ in self.genes() if self.offset(j) > 0)
		if q.strand == "-":
			left_avail, right_avail = right_avail, left_avail
			left_reach, right_reach = right_reach, left_reach
		return RangeReport(row, length, q.start, q.end, left_avail, right_avail,
			left_reach, right_reach, before, after)


class Outputs:
	def __init__(self, out: Path):
		self.out = out
		self.genes: list[Gene] = []
		self.scan_genes: list[Gene] = []
		self.windows: list[Window] = []
		self.rows: list[RowInfo] = []
		self.unmatched: list[Unmatched] = []
		self.ranges: list[RangeReport] = []
		self.proteins: dict[str, str] = {}
		self.scan_proteins: dict[str, str] = {}
		self.queries: dict[str, str] = {}
		self.rna: dict[str, str] = {}

	def write(self, scan_range: bool, order: dict[str, int]) -> None:
		rank = lambda row: order.get(row, len(order))
		self.genes.sort(key=lambda g: (rank(g.row_id), g.start))
		self.scan_genes.sort(key=lambda g: (rank(g.row_id), g.start))
		self.windows.sort(key=lambda w: rank(w.row_id))
		self.rows.sort(key=lambda r: rank(r.row_id))
		self.ranges.sort(key=lambda r: rank(r.row_id))
		Gene.write(self.out / Gene.FILE, self.genes)
		Window.write(self.out / Window.FILE, self.windows)
		RowInfo.write(self.out / RowInfo.FILE, self.rows)
		Unmatched.write(self.out / Unmatched.FILE, self.unmatched)
		RangeReport.write(self.out / RangeReport.FILE, self.ranges)
		fasta.write(self.out / PROTEINS, sorted(self.proteins.items()))
		fasta.write(self.out / QUERIES, self.queries.items())
		fasta.write(self.out / RNA, sorted(self.rna.items()))
		if scan_range:
			Gene.write(self.out / SCAN_GENES, self.scan_genes)
			fasta.write(self.out / SCAN_PROTEINS, sorted(self.scan_proteins.items()))


class GenomeSource:
	def __init__(self, files: GenomeFiles):
		self.files = files
		self.table = GeneTable(files.path("gff"))
		self.proteins = ProteinFasta(files.path("faa"))
		self.rna = RnaFasta(files.path("rna")) if files.rna != MISSING else None
		self._genome: Optional[GenomeFasta] = None

	def rna_sequence(self, gene: Feature) -> Optional[str]:
		seq = self.rna.get(gene) if self.rna else None
		if seq is None and self.files.genome != MISSING:
			if self._genome is None:
				self._genome = GenomeFasta(self.files.path("genome"))
			seq = self._genome.slice(gene)
		return seq


class Extract(Stage):
	name = "extract"
	requires = ("fetch",)

	def run(self, run, config, out: Path) -> None:
		settings = Settings.from_config(config)
		genomes = {g.assembly: g for g in GenomeFiles.read(run.stage_file("fetch", GenomeFiles.FILE))}
		targets = QueryTarget.read(run.stage_file("fetch", QueryTarget.FILE))
		outputs = Outputs(out)
		grouped: dict[str, list[QueryTarget]] = {}
		for target in targets:
			if target.status != "ok":
				outputs.unmatched.append(Unmatched(target.query, target.assembly, target.status))
				continue
			grouped.setdefault(target.assembly, []).append(target)
		for assembly in sorted(grouped):
			files = genomes.get(assembly)
			if files is None or not files.usable:
				for target in grouped[assembly]:
					outputs.unmatched.append(Unmatched(target.query, assembly, "genome files missing"))
				continue
			source = GenomeSource(files)
			for target in grouped[assembly]:
				index = source.table.find(target.acceptable)
				if index is None:
					outputs.unmatched.append(Unmatched(target.query, assembly, "not in the annotation"))
					continue
				self._one(source, target, index, settings, outputs)
		order = {row_id(t.query, t.assembly): i for i, t in enumerate(targets)}
		outputs.write(bool(settings.scan_range), order)
		note("extracted {} rows, {} genes; {} queries unmatched".format(
			len(outputs.rows), len(outputs.genes), len(outputs.unmatched)))
		if not outputs.rows:
			raise RuntimeError("no neighbourhood could be extracted for any query")

	def _one(self, source: GenomeSource, target: QueryTarget, index: int,
			settings: Settings, outputs: Outputs) -> None:
		hood = Neighbourhood(source.table, index, settings)
		row = row_id(target.query, target.assembly)
		q = hood.query
		for j, g in hood.genes():
			outputs.genes.append(Gene(row, target.assembly, g.contig, g.accession, g.start, g.end,
				g.strand, g.product or MISSING, g.is_rna, hood.offset(j)))
			if g.is_rna:
				seq = source.rna_sequence(g)
				if seq:
					outputs.rna.setdefault(g.accession, seq)
			else:
				seq = source.proteins.get(g.accession)
				if seq:
					outputs.proteins.setdefault(g.accession, seq)
		query_seq = source.proteins.get(q.accession)
		if query_seq:
			outputs.queries[row] = query_seq
		for j, g in hood.scan_genes():
			seq = source.proteins.get(g.accession)
			if seq:
				outputs.scan_proteins.setdefault(g.accession, seq)
				outputs.scan_genes.append(Gene(row, target.assembly, g.contig, g.accession, g.start,
					g.end, g.strand, g.product or MISSING, False, hood.offset(j)))
		lo, hi = hood.span()
		scan_lo, scan_hi = hood.scan_span()
		cut_lo, cut_hi = hood.scan_span(settings.margin)
		outputs.windows.append(Window(row, target.assembly, q.contig, hood.contig.length,
			lo, hi, scan_lo, scan_hi, cut_lo, cut_hi))
		outputs.rows.append(RowInfo(row, target.query, target.assembly, q.accession, q.contig,
			q.strand, source.proteins.organisms.get(q.accession, MISSING)))
		outputs.ranges.append(hood.range_report(row))
