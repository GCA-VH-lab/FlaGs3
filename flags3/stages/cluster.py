from pathlib import Path

from flags3 import cluster, fasta
from flags3.log import note
from flags3.schema import MISSING, Annotation, ClusterHit, Family, Gene, RowInfo
from flags3.stage import Stage
from flags3.tools import Tools


class FamilyTables:
	def __init__(self, families: list[list[str]], adjacency: dict[str, set],
			occurrences: dict[str, int], queries: set[str], prefix: str, tool: str):
		self.families = families
		self.adjacency = adjacency
		self.occurrences = occurrences
		self.queries = queries
		self.prefix = prefix
		self.tool = tool
		self.labels = self._labels()

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
		Family.write(out / Family.FILE, (
			Family(i + 1, self.labels.get(i, MISSING), len(fam), count(fam), ",".join(fam))
			for i, fam in enumerate(self.families)))
		of = {a: i + 1 for i, fam in enumerate(self.families) for a in fam}
		ClusterHit.write(out / ClusterHit.FILE, (
			ClusterHit(a, of[a], ",".join(sorted(self.adjacency.get(a, ()))) or MISSING)
			for a in sorted(of)))
		Annotation.write(out / Annotation.FILE, (
			Annotation(a, MISSING, None, None, "fill", "family:{}".format(i + 1), self.labels[i], self.tool, None)
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

	rna = False

	def run(self, run, config, out: Path) -> None:
		tools = Tools.load(config.path("tools"))
		name = config.text(self.method_key)
		clusterer = cluster.build(tools, name, config.workers(), out / "raw")
		sequences = self.sequences(run)
		note("clustering {} sequences with {}".format(len(sequences), name))
		families = clusterer.cluster(sequences)
		adjacency = clusterer.adjacency
		families, adjacency = self.extend(run, families, adjacency)
		queries = {r.accession for r in RowInfo.read(run.stage_file("extract", RowInfo.FILE))}
		tables = FamilyTables(families, adjacency, self.occurrences(run), queries, self.prefix, name)
		tables.write(out)
		note("{} families, {} of them shared".format(len(families), len(tables.labels)))

	def extend(self, run, families, adjacency):
		return families, adjacency


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
