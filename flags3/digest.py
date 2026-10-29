from pathlib import Path

from flags3.render.companions import Companions, figure_label
from flags3.render.data import RunData
from flags3.render.style import BAND_ORDER, STAGE_TITLES
from flags3.schema import MISSING, Annotation, Gene

FILE = "digest.md"
PSEUDOGENE = "pseudogene*"
SUBFAMILY = "abcdefghijklmnopqrstuvwxyz?"
MIN_ROWS = 2
MIN_FRACTION = 0.25
RULE_FRACTION = 0.75
OPTIONAL_STAGES = ("cluster", "cluster_rna", "domains", "features", "sismis", "genomad", "defence", "tree")


def offset_text(offset: int) -> str:
	return "{:+d}".format(offset) if offset else "0"


def counted(labels: list) -> str:
	counts: dict[str, int] = {}
	for label in labels:
		counts[label] = counts.get(label, 0) + 1
	return " ".join("{}×{}".format(k, n) if n > 1 else k for k, n in counts.items())


def commonest(values: list) -> str:
	counts: dict[str, int] = {}
	for value in values:
		if value:
			counts[value] = counts.get(value, 0) + 1
	return max(counts, key=lambda v: (counts[v], v)) if counts else ""


def family_rank(code: str):
	head = code.rstrip("0123456789" + SUBFAMILY)
	digits = "".join(ch for ch in code[len(head):] if ch.isdigit())
	return ({"Q": 0, "q": 1, "G": 2, "R": 3}.get(head, 4), int(digits or 0), code)


class Digest:
	def __init__(self, run, config, tables):
		self.run = run
		self.config = config
		self.t = tables
		self.data = RunData.load(run)
		self.companions = Companions(self.data, "occurrence")
		self.info = run.info()
		self.order = self.data.tree_order or self.data.order
		self.family_of = {a.subject: a for stage in ("cluster", "cluster_rna") for a in self.data.annotations.get(stage, [])}
		self.by_protein: dict[str, list[Annotation]] = {}
		for stage in ("domains", "features"):
			for a in self.data.annotations.get(stage, []):
				self.by_protein.setdefault(a.subject, []).append(a)
		self.bands: dict[tuple, list[tuple]] = {}
		self.verdicts: dict[tuple, list[str]] = {}
		for stage in BAND_ORDER:
			for a in self.data.annotations.get(stage, []):
				assembly, _, contig = a.subject.rpartition("|")
				if self.data.placement(stage, a) == "coordinates":
					self.bands.setdefault((assembly, contig), []).append((stage, a))
				else:
					label = self.data.verdict_label(stage, a)
					listed = self.verdicts.setdefault((assembly, contig), [])
					if label not in listed:
						listed.append(label)

	def write(self, path: Path) -> Path:
		lines = self.head() + self.read() + self.families() + self.systems()
		for query in self.queries():
			lines += self.block(query)
		lines += self.missing()
		path.write_text("\n".join(lines).rstrip() + "\n")
		return path

	def queries(self) -> list:
		seen, out = set(), []
		for row_id in self.order:
			query = self.data.rows[row_id].query
			if query not in seen:
				seen.add(query)
				out.append(query)
		return out

	def head(self) -> list:
		ran = [s for s in OPTIONAL_STAGES if s in self.data.present]
		not_run = [s for s in OPTIONAL_STAGES if s not in self.data.present]
		lines = ["# FlaGs3 digest · comparative genomic neighbourhood analysis",
			"run {} · {} queries · {} neighbourhoods · FlaGs3 {}".format(
				self.run.name, len(self.t.input_queries()), len(self.order), self.info.get("version"))]
		command = self.info.get("command")
		if command:
			lines.append("command: {}".format(command))
		lines.append("stages: {} · not run: {} (nothing here is a claim about them)".format(
			" ".join(ran) or "none beyond fetch and extract", " ".join(not_run) or "none"))
		return lines + [""]

	def read(self) -> list:
		c = self.config
		rng = c.integer("range")
		window = "every gene within {} bp of the query".format(rng) if rng else "query + up to {} genes each side".format(c.integer("gene", 4))
		lines = ["## Read",
			"neighbourhood = {}, from one assembly, turned so the query reads +".format(window),
			"gene = offset accession(family)strand · offset counts genes from the query (0) · strand is relative to the query",
			"family = code as on the figure, defined once in Families · pseudogene and rna say what they are",
			"- = no homologue elsewhere in this run; that gene carries (-;product;domains features) inline",
			"gap = bases to the previous gene in the line, one per gene from the second · negative = overlap",
			"systems = code:first..last offset covered · codes in Systems · nothing listed = the scanner called none there",
			"edge = the contig ended inside the window; genes missing there are unknown, not absent",
			"contig = what a scanner called the whole contig the neighbourhood sits on (a plasmid replicon type, or plasmid or "
			"virus); a call names the contig, not any one gene, and two scanners may disagree on the same contig",
			"genera = assemblies per genus of this query",
			"conserved = offset family present/total sameStrand/present · listed at >={} neighbourhoods and >={:.0f} %".format(MIN_ROWS, MIN_FRACTION * 100),
			"dev = exception to a >={:.0f} % pattern: the neighbourhood, and what it holds instead".format(RULE_FRACTION * 100),
			"a family means sequence similarity within this run, not shared function",
			"everything here is computed by FlaGs3 from the run; anything a reader adds is interpretation, not a FlaGs3 result"]
		return lines + [""]

	def families(self) -> list:
		members: dict[str, list] = {}
		for row_id in self.order:
			for g in self.data.genes[row_id]:
				code = self.family_code(g)
				if code not in (MISSING, "pseudogene"):
					members.setdefault(code, []).append(g)
		lines = ["## Families"]
		for code in sorted(members, key=family_rank):
			genes = members[code]
			parts = [code, commonest([g.product for g in genes if g.product != MISSING]) or "unnamed"]
			for extra in (commonest([self.domains(g) for g in genes]), commonest([self.features(g) for g in genes])):
				if extra:
					parts.append(extra)
			parts.append("{} genes".format(len(genes)))
			if code.startswith("Q"):
				parts.append("query family")
			lines.append(" · ".join(parts))
		return lines + [""]

	def systems(self) -> list:
		total = len(self.order)
		lines = ["## Systems"]
		for stage in BAND_ORDER:
			for category, code in self.companions.codes.get(stage, {}).items():
				tools = sorted({a.tool for a in self.data.annotations[stage] if a.category == category})
				touched = {r for r in self.order for s, a in self.row_bands(r) if s == stage and a.category == category}
				lines.append("{} · {} · {} · {} · {}/{} neighbourhoods".format(
					code, category, STAGE_TITLES.get(stage, stage).lower(), " ".join(tools), len(touched), total))
		domains = self.companions.codes.get("domains")
		if len(lines) == 1:
			lines.append("none called")
		if domains:
			lines.append("domain numbers on the figure are not repeated here; domains are named in full above and in each gene")
		return lines + [""]

	def block(self, query: str) -> list:
		rows = [r for r in self.order if self.data.rows[r].query == query]
		genera: dict[str, int] = {}
		for row_id in rows:
			species = self.t.species(row_id)
			genera[species.split()[0] if species else "unknown"] = genera.get(species.split()[0] if species else "unknown", 0) + 1
		lines = ["## {} · {} neighbourhood{} · genera {}".format(query, len(rows), "" if len(rows) == 1 else "s",
			" ".join("{} {}".format(g, n) for g, n in sorted(genera.items(), key=lambda kv: (-kv[1], kv[0]))))]
		if len(rows) >= MIN_ROWS:
			lines += self.conserved(rows) + self.deviations(rows)
		for index, row_id in enumerate(rows, 1):
			lines += self.neighbourhood(index, row_id)
		return lines + [""]

	def held(self, rows: list) -> dict:
		out: dict[tuple, dict] = {}
		for row_id in rows:
			for g in self.data.genes[row_id]:
				code = self.family_code(g)
				if g.offset == 0 or code in (MISSING, "pseudogene"):
					continue
				out.setdefault((g.offset, code.rstrip(SUBFAMILY)), {})[row_id] = self.t.normalised_strand(row_id, g)
		return out

	def conserved(self, rows: list) -> list:
		kept = [(k, v) for k, v in sorted(self.held(rows).items()) if len(v) >= MIN_ROWS and len(v) / len(rows) >= MIN_FRACTION]
		if not kept:
			return ["conserved: none at >={} neighbourhoods and >={:.0f} %".format(MIN_ROWS, MIN_FRACTION * 100)]
		return ["conserved: " + " ; ".join("{} {} {}/{} {}/{}".format(offset_text(offset), code, len(who), len(rows),
			sum(1 for s in who.values() if s == "+"), len(who)) for (offset, code), who in kept)]

	def deviations(self, rows: list) -> list:
		total = len(rows)
		out = []
		for (offset, code), who in sorted(self.held(rows).items()):
			if len(who) < MIN_ROWS or len(who) / total < RULE_FRACTION:
				continue
			for row_id in rows:
				if row_id not in who:
					out.append("dev: {} {} {}/{} · {} holds {}".format(offset_text(offset), code, len(who), total,
						self.name(row_id), self.holds(row_id, offset)))
			same = [r for r, strand in who.items() if strand == "+"]
			opposite = [r for r, strand in who.items() if strand == "-"]
			for majority, minority, sign in ((same, opposite, "+"), (opposite, same, "-")):
				if minority and len(majority) / len(who) >= RULE_FRACTION:
					out.append("dev: {} {} strand {} in {}/{} · flipped in {}".format(offset_text(offset), code, sign,
						len(majority), len(who), ", ".join(self.name(r) for r in minority)))
		hit: dict[tuple, set] = {}
		for row_id in rows:
			for stage, a in self.row_bands(row_id):
				hit.setdefault((stage, a.category), set()).add(row_id)
		for (stage, category), touched in sorted(hit.items(), key=lambda kv: (BAND_ORDER.index(kv[0][0]), kv[0][1])):
			code = self.companions.codes[stage][category]
			if len(touched) < total and len(touched) / total >= RULE_FRACTION:
				out.append("dev: {} {}/{} · absent in {}".format(code, len(touched), total,
					", ".join(self.name(r) for r in rows if r not in touched)))
			elif len(touched) / total <= 1 - RULE_FRACTION:
				out.append("dev: {} only in {}".format(code, ", ".join(self.name(r) for r in sorted(touched))))
		return out

	def name(self, row_id: str) -> str:
		species = self.t.species(row_id)
		assembly = self.data.rows[row_id].assembly
		return "{} {}".format(assembly, species) if species else assembly

	def holds(self, row_id: str, offset: int) -> str:
		gene = next((g for g in self.data.genes[row_id] if g.offset == offset), None)
		return self.token(row_id, gene) if gene is not None else "nothing (the window ends)"

	def token(self, row_id: str, gene: Gene) -> str:
		product = gene.product if gene.product != MISSING else "unnamed"
		code = self.family_code(gene)
		if code == "pseudogene":
			return "pseudogene({})".format(product)
		if gene.is_rna and code == MISSING:
			return "rna({})".format(product)
		if code == MISSING:
			inline = ";".join(x for x in (product, self.domains(gene), self.features(gene)) if x)
			return "{}(-;{})".format(gene.accession, inline)
		return "{}({})".format(gene.accession, code)

	def neighbourhood(self, index: int, row_id: str) -> list:
		info = self.data.rows[row_id]
		genes = sorted(self.data.genes[row_id], key=lambda g: g.offset)
		species = self.t.species(row_id)
		head = "N{} {}{} · contig {} · query on {}{}".format(index, info.assembly, " " + species if species else "",
			info.contig, info.strand, " · edge" if any(g.contig_edge for g in genes) else "")
		tokens, gaps, previous = [], [], None
		for g in genes:
			tokens.append("{} {}{}".format(offset_text(g.offset), self.token(row_id, g), self.t.normalised_strand(row_id, g)))
			if previous is not None:
				earlier, later = sorted((previous, g), key=lambda x: x.start)
				gaps.append(str(later.start - earlier.end - 1))
			previous = g
		called = self.verdicts.get((info.assembly, info.contig), [])
		tail = "  gap {}".format(" ".join(gaps))
		if called:
			tail += " · contig {}".format(" ".join(called))
		covered = []
		for stage, a in self.row_bands(row_id):
			offsets = [g.offset for g in genes if g.end >= a.start and g.start <= a.end]
			covered.append("{}:{}..{}".format(self.companions.codes[stage][a.category], offset_text(min(offsets)), offset_text(max(offsets))))
		if covered:
			tail += " · systems {}".format(" ".join(covered))
		return [head, "  " + " ".join(tokens), tail]

	def row_bands(self, row_id: str) -> list:
		genes = self.data.genes[row_id]
		lo, hi = min(g.start for g in genes), max(g.end for g in genes)
		out = []
		for stage, a in self.bands.get((genes[0].assembly, genes[0].contig), []):
			if a.end >= lo and a.start <= hi and any(g.end >= a.start and g.start <= a.end for g in genes):
				out.append((stage, a))
		return out

	def family_code(self, gene: Gene) -> str:
		if gene.accession == PSEUDOGENE:
			return "pseudogene"
		a = self.family_of.get(gene.accession)
		return figure_label(a.label) if a and a.label != MISSING else MISSING

	def domains(self, gene: Gene) -> str:
		return counted([a.label for a in self.by_protein.get(gene.accession, []) if a.kind == "wedge"])

	def features(self, gene: Gene) -> str:
		return counted([a.label for a in self.by_protein.get(gene.accession, []) if a.kind != "wedge"])

	def missing(self) -> list:
		found = {self.data.rows[r].query for r in self.order}
		lines = []
		for query in self.t.input_queries():
			if query in found:
				continue
			reasons = [u.reason for u in self.t.unmatched if u.query == query] or [
				x.status for x in self.t.targets if x.query == query and x.status != "ok"]
			lines.append("{} · {}".format(query, reasons[0] if reasons else "no neighbourhood"))
		return ["## No neighbourhood"] + lines if lines else []
