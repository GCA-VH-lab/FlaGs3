from pathlib import Path
from typing import Optional

from flags3 import cluster, fasta
from flags3.log import note
from flags3.schema import MISSING, Annotation, ClusterHit, Family, Gene, RowInfo
from flags3.stage import Stage
from flags3.tools import Tools


class FamilyTables:
	def __init__(self, families: list[list[str]], adjacency: dict[str, set],
			occurrences: dict[str, int], queries: set[str], prefix: str, tool: str,
			full: Optional[dict[str, set]] = None, contains: Optional[dict[str, set]] = None,
			edge: Optional[set] = None):
		self.families = families
		self.adjacency = adjacency
		self.full = full or {}
		self.contains = contains or {}
		self.edge = edge or set()
		self.occurrences = occurrences
		self.queries = queries
		self.prefix = prefix
		self.tool = tool
		self.labels = self._labels()
		self.sub: dict[int, tuple[dict[str, str], list[list[str]]]] = {}
		if self.full:
			for i, fam in enumerate(self.families):
				letters, groups = cluster.subfamilies(fam, self.full, self.contains)
				if groups:
					self.sub[i] = (letters, groups)

	def letters(self, index: int, accession: str) -> str:
		if index not in self.sub:
			return ""
		if accession in self.edge:
			return "?"
		return self.sub[index][0].get(accession, "")

	def member_label(self, index: int, accession: str) -> str:
		base = self.labels.get(index, MISSING)
		if base == MISSING:
			return base
		return base + self.letters(index, accession)

	def member_category(self, index: int, accession: str) -> str:
		letter = self.letters(index, accession)
		if len(letter) == 1 and letter != "?":
			return "family:{}/{}".format(index + 1, letter)
		return "family:{}".format(index + 1)

	def _labels(self) -> dict[int, str]:
		count = lambda fam: sum(self.occurrences.get(a, 0) for a in fam)
		shared = [i for i, fam in enumerate(self.families) if count(fam) > 1]
		shared.sort(key=lambda i: (-count(self.families[i]), self.families[i][0]))
		labels, plain, query = {}, 0, 0
		for i in shared:
			if self.queries.intersection(self.families[i]):
				query += 1
				labels[i] = "Q{}".format(query)
			elif self.prefix:
				plain += 1
				labels[i] = "{}{}".format(self.prefix, plain)
			else:
				plain += 1
				labels[i] = str(plain)
		return labels

	def write(self, out: Path) -> None:
		count = lambda fam: sum(self.occurrences.get(a, 0) for a in fam)
		def sub_text(i):
			if i not in self.sub:
				return MISSING, MISSING
			letters, groups = self.sub[i]
			groups_text = ";".join("{}:{}".format(cluster.LETTERS[j % 26] * (j // 26 + 1), len(g)) for j, g in enumerate(groups))
			bridges = ";".join("{}:{}".format(m, self.letters(i, m)) for m in sorted(letters)
				if len(self.letters(i, m)) > 1 or self.letters(i, m) == "?")
			return groups_text, bridges or MISSING

		Family.write(out / Family.FILE, (
			Family(i + 1, self.labels.get(i, MISSING), len(fam), count(fam), ",".join(fam), *sub_text(i))
			for i, fam in enumerate(self.families)))
		of = {a: i + 1 for i, fam in enumerate(self.families) for a in fam}
		ClusterHit.write(out / ClusterHit.FILE, (
			ClusterHit(a, of[a], ",".join(sorted(self.adjacency.get(a, ()))) or MISSING,
				",".join(sorted(self.full.get(a, ()))) or MISSING)
			for a in sorted(of)))
		Annotation.write(out / Annotation.FILE, (
			Annotation(a, MISSING, None, None, "fill", self.member_category(i, a), self.member_label(i, a), self.tool, None)
			for i in sorted(self.labels) for a in self.families[i]))


class Cluster(Stage):
	name = "cluster"
	requires = ("extract",)
	method_key = "cluster_method"
	prefix = ""

	def sequences(self, run) -> dict[str, str]:
		return {name: seq for name, _, seq in fasta.read(run.stage_file("extract", "proteins.faa"))}

	def occurrences(self, run) -> dict[str, int]:
		counts: dict[str, int] = {}
		for g in Gene.iterate(run.stage_file("extract", Gene.FILE)):
			if g.is_rna == self.rna:
				counts[g.accession] = counts.get(g.accession, 0) + 1
		return counts

	def edge_accessions(self, run) -> set:
		seen_inside, seen = set(), set()
		for g in Gene.iterate(run.stage_file("extract", Gene.FILE)):
			if g.is_rna != self.rna:
				continue
			seen.add(g.accession)
			if not g.contig_edge:
				seen_inside.add(g.accession)
		return seen - seen_inside

	rna = False

	def run(self, run, config, out: Path) -> None:
		tools = Tools.load(config.path("tools"))
		name = config.text(self.method_key)
		self.override(tools, name, config)
		clusterer = cluster.build(tools, name, config.workers(), out / "raw")
		sequences = self.sequences(run)
		note("clustering {} sequences with {}".format(len(sequences), name))
		families = clusterer.cluster(sequences)
		adjacency = clusterer.adjacency
		families, adjacency = self.extend(run, families, adjacency)
		threshold = config.number("subfamily_coverage")
		if threshold is None:
			threshold = cluster.option(cluster.method(tools, name), "subcov", 0.0)
		full = clusterer.full_length(threshold) if threshold > 0 else {}
		contains = clusterer.contains(threshold) if threshold > 0 else {}
		queries = {r.accession for r in RowInfo.read(run.stage_file("extract", RowInfo.FILE))}
		tables = FamilyTables(families, adjacency, self.occurrences(run), queries, self.prefix, name, full, contains,
			self.edge_accessions(run))
		tables.write(out)
		split = [i for i in tables.sub if i in tables.labels]
		note("{} families, {} of them shared{}".format(len(families), len(tables.labels),
			", {} split into subfamilies at {:.0%} mutual coverage".format(len(split), threshold) if split else ""))

	def extend(self, run, families, adjacency):
		return families, adjacency

	@staticmethod
	def override(tools: Tools, name: str, config) -> None:
		tool = cluster.method(tools, name)
		hmmer = tool.engine in ("jackhmmer", "nhmmer")
		if config.text("iterations"):
			if hmmer:
				tool.options["iterations"] = config.text("iterations")
				note("{}: iterations = {} from the command line".format(name, config.text("iterations")))
			else:
				print("Warning: -n/--iterations has no meaning for {} ({}); ignored.".format(name, tool.engine))
		if config.text("cluster_evalue"):
			tool.options["incE" if hmmer else "evalue"] = config.text("cluster_evalue")
			note("{}: E-value = {} from the command line".format(name, config.text("cluster_evalue")))


class ClusterRna(Cluster):
	name = "cluster_rna"
	needs = ("rna",)
	method_key = "rna_method"
	prefix = "R"
	rna = True

	def wanted(self, config) -> bool:
		return config.flag("cluster_rna")

	def sequences(self, run) -> dict[str, str]:
		return {name: seq for name, _, seq in fasta.read(run.stage_file("extract", "rna.fna"))}

	def extend(self, run, families, adjacency):
		clustered = {a for fam in families for a in fam}
		products = {g.accession: g.product for g in Gene.iterate(run.stage_file("extract", Gene.FILE))
			if g.is_rna and g.accession not in clustered}
		if products:
			print("Warning: {} flanking RNAs had no sequence and were grouped by product name.".format(len(products)))
			extra = cluster.by_name(products)
			families = families + extra
			for fam in extra:
				for a in fam:
					adjacency[a] = set(fam) - {a}
		return families, adjacency
