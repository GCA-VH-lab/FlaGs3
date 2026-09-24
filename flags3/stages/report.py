import time
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from flags3 import fasta
from flags3.log import note
from flags3.schema import (MISSING, Family, Gene, GenomeFiles, QueryTarget, RangeReport, Row, RowInfo,
	Unmatched, Window, split_row_id)
from flags3.stage import Stage
from flags3.stages.extract import PROTEINS, QUERIES
from flags3.stages.domains import Domain
from flags3.inputs import InputList

LEGACY = "legacy"


@dataclass(frozen=True)
class NeighbourhoodRow(Row):
	FILE: ClassVar[str] = "neighbourhoods.tsv"
	row_id: str
	query: str
	assembly: str
	species: str
	offset: int
	accession: str
	family: str
	contig: str
	start: int
	end: int
	strand: str
	length: int
	is_rna: bool
	product: str
	domains: str


@dataclass(frozen=True)
class QueryStatus(Row):
	FILE: ClassVar[str] = "queries.tsv"
	query: str
	assemblies: str
	rows: int
	status: str


class RunTables:
	def __init__(self, run):
		self.run = run
		self.genes = Gene.read(run.stage_file("extract", Gene.FILE))
		self.rows = {r.row_id: r for r in RowInfo.read(run.stage_file("extract", RowInfo.FILE))}
		self.windows = {w.row_id: w for w in Window.read(run.stage_file("extract", Window.FILE))}
		self.ranges = {r.row_id: r for r in RangeReport.read(run.stage_file("extract", RangeReport.FILE))}
		self.unmatched = Unmatched.read(run.stage_file("extract", Unmatched.FILE))
		self.targets = QueryTarget.read(run.stage_file("fetch", QueryTarget.FILE))
		self.by_row: dict[str, list[Gene]] = {}
		for g in self.genes:
			self.by_row.setdefault(g.row_id, []).append(g)
		for genes in self.by_row.values():
			genes.sort(key=lambda g: g.offset)
		self.occurrences: dict[str, int] = {}
		self.products: dict[str, str] = {}
		for g in self.genes:
			self.occurrences[g.accession] = self.occurrences.get(g.accession, 0) + 1
			self.products.setdefault(g.accession, "" if g.product == MISSING else g.product)
		self.families: list[Family] = []
		for stage in ("cluster", "cluster_rna"):
			if run.has(stage):
				self.families += Family.read(run.stage_file(stage, Family.FILE))
		self.label_of = {a: f.label for f in self.families if f.label != MISSING for a in f.accessions}
		self.domains: dict[str, list[str]] = {}
		if run.has("domains"):
			for d in Domain.read(run.stage_file("domains", Domain.FILE)):
				self.domains.setdefault(d.protein, []).append(d.domain)
		self.proteins = {n: s for n, _, s in fasta.read(run.stage_file("extract", PROTEINS))}
		self.queries = {n: s for n, _, s in fasta.read(run.stage_file("extract", QUERIES))}
		self.query_accessions = {r.accession for r in self.rows.values()}

	def family(self, accession: str) -> str:
		return self.label_of.get(accession, MISSING)

	def species(self, row_id: str) -> str:
		s = self.rows[row_id].species if row_id in self.rows else MISSING
		return "" if s == MISSING else s

	def normalised_strand(self, row_id: str, gene: Gene) -> str:
		if self.rows[row_id].strand == "-":
			return "-" if gene.strand == "+" else "+"
		return gene.strand

	def input_queries(self) -> list[str]:
		seen, out = set(), []
		for t in self.targets:
			if t.query not in seen:
				seen.add(t.query)
				out.append(t.query)
		return out


class Report(Stage):
	name = "report"
	requires = ("extract",)
	optional = True
	barrier = True

	def run(self, run, config, out: Path) -> None:
		t = RunTables(run)
		self.neighbourhoods(t, out)
		self.queries(t, out)
		self.families(t, out)
		self.fastas(t, out)
		self.summary(run, t, out)
		legacy = out / LEGACY
		legacy.mkdir(exist_ok=True)
		prefix = run.name
		Legacy(t, legacy, prefix).write_all()
		note("report: {} rows, {} queries, legacy tables under report/{}".format(len(t.rows), len(t.input_queries()), LEGACY))

	def neighbourhoods(self, t: RunTables, out: Path) -> None:
		rows = []
		for row_id, genes in t.by_row.items():
			query, assembly = split_row_id(row_id)
			for g in genes:
				rows.append(NeighbourhoodRow(row_id, query, assembly, t.species(row_id) or MISSING, g.offset, g.accession,
					t.family(g.accession), g.contig, g.start, g.end, g.strand, g.end - g.start + 1, g.is_rna,
					g.product, ";".join(t.domains.get(g.accession, [])) or MISSING))
		NeighbourhoodRow.write(out / NeighbourhoodRow.FILE, rows)

	def queries(self, t: RunTables, out: Path) -> None:
		rows = []
		for query in t.input_queries():
			targets = [x for x in t.targets if x.query == query]
			assemblies = [x.assembly for x in targets if x.assembly != MISSING]
			found = [r for r in t.rows if split_row_id(r)[0] == query]
			if found:
				status = "ok"
			else:
				reasons = [u.reason for u in t.unmatched if u.query == query] or [x.status for x in targets if x.status != "ok"]
				status = reasons[0] if reasons else "no neighbourhood"
			rows.append(QueryStatus(query, ";".join(assemblies) or MISSING, len(found), status))
		QueryStatus.write(out / QueryStatus.FILE, rows)

	def families(self, t: RunTables, out: Path) -> None:
		with open(out / "families.tsv", "w") as handle:
			handle.write("family\tlabel\tsize\toccurrences\tproduct\tmembers\n")
			for f in t.families:
				product = next((t.products[a] for a in f.accessions if t.products.get(a)), MISSING)
				handle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(f.family, f.label, f.size, f.occurrences, product, f.members))

	def fastas(self, t: RunTables, out: Path) -> None:
		def describe(accession: str) -> str:
			product = t.products.get(accession, "")
			return "{} {}".format(accession, product) if product else accession

		flanking = [(describe(a), s) for a, s in sorted(t.proteins.items()) if a not in t.query_accessions]
		queries = [("{} {}".format(row, t.species(row)).rstrip(), s) for row, s in t.queries.items()]
		fasta.write(out / "flanking.faa", flanking)
		fasta.write(out / "queries.faa", queries)
		fasta.write(out / "all.faa", queries + flanking)

	def summary(self, run, t: RunTables, out: Path) -> None:
		info = run.info()
		started = time.strptime(info.get("started"), "%Y-%m-%d %H:%M:%S") if info.get("started") else None
		elapsed = time.time() - time.mktime(started) if started else None
		lines = ["FlaGs3 {}".format(info.get("version")), "started: {}".format(info.get("started")),
			"command: {}".format(info.get("command"))]
		if elapsed is not None:
			lines.append("elapsed: {} (wall time to this report; figures follow)".format(_human(elapsed)))
		lines.append("")
		lines.append("{} queries, {} rows, {} genes, {} unmatched".format(
			len(t.input_queries()), len(t.rows), len(t.genes), len(t.unmatched)))
		lines.append("{} families, {} labelled".format(len(t.families), len({f.label for f in t.families if f.label != MISSING})))
		lines.append("")
		lines.append("stage        status    seconds")
		total = 0.0
		for stage in sorted(run.stages(), key=lambda s: run.status(s).get("started", "")):
			status = run.status(stage)
			if stage == self.name:
				continue
			total += float(status.get("seconds") or 0)
			lines.append("{:12} {:9} {:>7}{}".format(stage, status.state or "-", status.get("seconds", "-"),
				"   " + status.get("error") if status.get("error") else ""))
		lines.append("{:12} {:9} {:>7}".format("stages total", "", "{:.2f}".format(total)))
		(out / "run_summary.txt").write_text("\n".join(lines) + "\n")


def _human(seconds: float) -> str:
	seconds = int(seconds)
	if seconds < 60:
		return "{}s".format(seconds)
	if seconds < 3600:
		return "{}m {:02d}s".format(seconds // 60, seconds % 60)
	return "{}h {:02d}m {:02d}s".format(seconds // 3600, seconds % 3600 // 60, seconds % 60)


class Legacy:
	def __init__(self, t: RunTables, directory: Path, prefix: str):
		self.t = t
		self.directory = directory
		self.prefix = prefix

	def path(self, suffix: str) -> Path:
		return self.directory / (self.prefix + suffix)

	def write_all(self) -> None:
		self.operon()
		self.clusters()
		self.outdesc()
		self.species_info()
		self.query_status()
		self.flankgene_report()
		self.fastas()
		self.range_report()
		self.accession_issues()

	def operon(self) -> None:
		t = self.t
		with open(self.path("_operon.tsv"), "w") as out:
			out.write("#query\tassembly\tspecies\tfamily\tstrand\toffset\tstart\tend\tlength\tcontig\tis_rna\taccession\tproduct\n")
			for row_id, genes in t.by_row.items():
				query, assembly = split_row_id(row_id)
				for g in genes:
					out.write("\t".join(str(v) for v in (query, assembly, t.species(row_id), t.family(g.accession),
						t.normalised_strand(row_id, g), g.offset, g.start, g.end, g.end - g.start + 1, g.contig,
						"True" if g.is_rna else "False", g.accession, "" if g.product == MISSING else g.product)) + "\n")

	def clusters(self) -> None:
		with open(self.path("_clusters.tsv"), "w") as out:
			out.write("#family\tsize\tmembers\n")
			for f in self.t.families:
				out.write("{}\t{}\t{}\n".format(f.label, f.size, f.members))

	def outdesc(self) -> None:
		t = self.t
		shared = [f for f in t.families if f.label != MISSING]
		shared.sort(key=lambda f: (-f.occurrences, f.accessions[0]))
		with open(self.path("_outdesc.txt"), "w") as out:
			for f in shared:
				for acc in f.accessions:
					out.write("{}({})\t{}\t{}\n".format(f.label, t.occurrences.get(acc, 0), acc, t.products.get(acc, "")))
				out.write("\n\n")

	def species_info(self) -> None:
		with open(self.path("_speciesInfo.txt"), "w") as out:
			out.write("#query\tassembly\tspecies\n")
			for row_id in sorted(self.t.rows):
				query, assembly = split_row_id(row_id)
				out.write("{}\t{}\t{}\n".format(query, assembly, self.t.species(row_id)))

	def query_status(self) -> None:
		t = self.t
		with open(self.path("_QueryStatus.txt"), "w") as out:
			out.write("#query\tassemblies\tflanking_genes_found\n")
			for query in t.input_queries():
				assemblies = [x.assembly for x in t.targets if x.query == query and x.assembly != MISSING]
				found = any(split_row_id(r)[0] == query for r in t.rows)
				out.write("{}\t{}\t{}\n".format(query, ";".join(assemblies) or "-", "Yes" if found else "No"))

	def flankgene_report(self) -> None:
		t = self.t
		with open(self.path("_flankgene_Report.log"), "w") as out:
			for row_id, genes in t.by_row.items():
				chain = " ".join("{}({})".format(g.accession, t.family(g.accession)) for g in genes)
				out.write("{}\t{}\n".format(row_id, chain))

	def fastas(self) -> None:
		t = self.t
		label = lambda acc: "{}|{}".format(acc, t.products[acc]) if t.products.get(acc) else acc
		rows_out = [(row, t.queries[row]) for row in sorted(t.queries)]
		flank_out = [(label(a), t.proteins[a]) for a in sorted(t.proteins) if a not in t.query_accessions]
		fasta.write(self.path("_tree.fasta"), rows_out, width=0)
		fasta.write(self.path("_flankgene.fasta"), flank_out, width=0)
		fasta.write(self.path("_all.fasta"), rows_out + flank_out, width=0)

	def range_report(self) -> None:
		t = self.t
		with open(self.path("_rangeReport.tsv"), "w") as out:
			out.write("#query\tassembly\tcontig\tcontig_length\tquery_start\tquery_end\tquery_strand\trequested_bp\t"
				"up_available\tdown_available\tup_reached\tdown_reached\ttruncated\tgenes_up\tgenes_down\tgenes_total\t"
				"scan_start\tscan_end\tscan_span\n")
			for row_id in t.by_row:
				r, w = t.ranges.get(row_id), t.windows.get(row_id)
				if r is None or w is None:
					continue
				query, assembly = split_row_id(row_id)
				out.write("\t".join(str(v) for v in (query, assembly, w.contig, r.contig_length, r.query_start, r.query_end,
					t.rows[row_id].strand, "-", r.up_available, r.down_available, r.up_reached, r.down_reached, "-",
					r.genes_up, r.genes_down, r.genes_up + r.genes_down + 1, w.scan_lo, w.scan_hi, w.scan_hi - w.scan_lo + 1)) + "\n")

	def accession_issues(self) -> None:
		t = self.t
		with open(self.path("_accessionIssues.txt"), "w") as out:
			out.write("#query\tissue\n")
			for query in t.input_queries():
				targets = [x for x in t.targets if x.query == query]
				if any(split_row_id(r)[0] == query for r in t.rows):
					continue
				if not any(x.assembly != MISSING for x in targets):
					out.write("{}\tno assembly resolved (not found locally or via NCBI)\n".format(query))
				else:
					out.write("{}\tassembly resolved but no flanking neighborhood extracted\n".format(query))
