import csv
import re
from pathlib import Path

from flags3.render.data import RunData
from flags3.render.style import BAND_CODES, BAND_ORDER, STAGE_TITLES, Colours
from flags3.schema import MISSING, Annotation, Gene

LEGEND = "legend.tsv"
PROTEIN_LEGEND = "protein_clusters_legend.txt"
RNA_LEGEND = "rna_clusters_legend.txt"
SYSTEMS = "systems.tsv"
LETTERS = re.compile(r"^[QqR]?\d+")


def figure_label(label: str) -> str:
	return "G" + label if label[:1].isdigit() else label


class Companions:
	def __init__(self, data: RunData, numbering: str):
		self.data = data
		self.colours = Colours("bright", numbering=numbering)
		self.codes: dict[str, dict[str, str]] = {}
		for stage, annotations in data.annotations.items():
			table = self.colours.assign(stage, annotations)
			if stage in BAND_CODES:
				self.codes[stage] = {c: "{}{}".format(BAND_CODES[stage], i + 1) for i, c in enumerate(table)}
			elif stage == "domains":
				self.codes[stage] = {c: str(i + 1) for i, c in enumerate(table)}

	def write_all(self, out: Path, stamped=lambda name: name) -> list[Path]:
		written = [self.legend(out / stamped(LEGEND))]
		if "cluster" in self.data.present:
			written.append(self.clusters_legend(out / stamped(PROTEIN_LEGEND), "cluster"))
		if "cluster_rna" in self.data.present:
			written.append(self.clusters_legend(out / stamped(RNA_LEGEND), "cluster_rna"))
		if any(s in self.data.annotations for s in BAND_ORDER):
			written.append(self.systems(out / stamped(SYSTEMS)))
		return written

	def legend(self, path: Path) -> Path:
		with open(path, "w", newline="", encoding="utf-8") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["stage", "code", "name", "occurrences"])
			for stage, codes in self.codes.items():
				counts: dict[str, int] = {}
				labels: dict[str, list] = {}
				for a in self.data.annotations[stage]:
					counts[a.category] = counts.get(a.category, 0) + 1
					if a.label not in labels.setdefault(a.category, []):
						labels[a.category].append(a.label)
				for category, code in codes.items():
					names = labels.get(category, [])
					name = category if category in names or len(names) > 3 else "{} ({})".format(category, ", ".join(names)) if names else category
					writer.writerow([stage, code, name, counts.get(category, 0)])
		return path

	def clusters_legend(self, path: Path, stage: str) -> Path:
		occurrences: dict[str, int] = {}
		products: dict[str, str] = {}
		for row in self.data.genes.values():
			for g in row:
				occurrences[g.accession] = occurrences.get(g.accession, 0) + 1
				products.setdefault(g.accession, "" if g.product == MISSING else g.product)
		label_of: dict[str, str] = {}
		groups: dict[str, list[str]] = {}
		for a in self.data.annotations.get(stage, []):
			label_of[a.subject] = a.label
			groups.setdefault(a.category.split("/")[0], []).append(a.subject)

		def letters(acc: str) -> str:
			return LETTERS.sub("", label_of.get(acc, ""))

		def rank(acc: str):
			tail = letters(acc)
			kind = 0 if len(tail) == 1 and tail != "?" else (2 if tail == "?" else (1 if tail else 0))
			return (kind, tail, -occurrences.get(acc, 0), acc)

		ordered = sorted(groups.items(), key=lambda kv: (-sum(occurrences.get(m, 0) for m in kv[1]), kv[1][0]))
		with open(path, "w", encoding="utf-8") as out:
			for _, members in ordered:
				for acc in sorted(members, key=rank):
					out.write("{}({})\t{}\t{}\n".format(label_of.get(acc, "-"), occurrences.get(acc, 0), acc, products.get(acc, "")))
				out.write("\n\n")
		return path

	def systems(self, path: Path) -> Path:
		label_of = {a.subject: a.label for stage in ("cluster", "cluster_rna") for a in self.data.annotations.get(stage, [])}
		by_contig: dict[tuple, list[Gene]] = {}
		for row_id, genes in self.data.genes.items():
			for g in genes:
				by_contig.setdefault((g.assembly, g.contig), []).append(g)
		with open(path, "w", newline="", encoding="utf-8") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["stage", "code", "system", "assembly", "contig", "start", "end", "rows", "genes", "products"])
			for stage in BAND_ORDER:
				for a in self.data.annotations.get(stage, []):
					assembly, _, contig = a.subject.rpartition("|")
					inside = [g for g in by_contig.get((assembly, contig), []) if g.end >= a.start and g.start <= a.end]
					rows = sorted({g.row_id for g in inside})
					if not rows:
						continue
					seen, genes, products = set(), [], []
					for g in sorted(inside, key=lambda g: g.start):
						if g.accession in seen:
							continue
						seen.add(g.accession)
						genes.append("{}({})".format(g.accession, figure_label(label_of.get(g.accession, "-"))))
						products.append("" if g.product == MISSING else g.product)
					writer.writerow([stage, self.codes[stage][a.category], a.label, assembly, contig, a.start, a.end,
						",".join(rows), " ".join(genes), "; ".join(products)])
		return path
