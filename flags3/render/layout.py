from dataclasses import dataclass
from typing import Optional

from flags3.render import text
from flags3.render.data import RunData
from flags3.render.style import BAND_ORDER, FigureSpec
from flags3.schema import Gene


@dataclass
class Arrow:
	x0: float
	x1: float
	cy: float
	strand: str
	points: list[tuple[float, float]]
	clip: str

	@property
	def width(self) -> float:
		return self.x1 - self.x0

	def polygon(self) -> str:
		return " ".join("{:.1f},{:.1f}".format(px, py) for px, py in self.points)


class Layout:
	def __init__(self, data: RunData, spec: FigureSpec, rows: list[str], band_stages: list[str]):
		self.data = data
		self.spec = spec
		self.rows = rows
		self.band_stages = [s for s in BAND_ORDER if s in band_stages]
		self.font = spec.value("font_size")
		self.row_h = spec.value("row_height")
		self.gene_h = spec.value("gene_height")
		self.pad = spec.value("pad")
		self.bp_per_px = spec.value("bases_per_pixel")
		self.band_h = self.gene_h + 4
		self.reversed = {r: data.rows[r].strand == "-" for r in rows}
		self.anchor = {r: self._anchor(r) for r in rows}
		self.labels = data.labels()
		self.tree_w = spec.value("tree_width") if (spec.tree and data.newick) else 0
		self.label_x = self.pad + self.tree_w
		needed = [self.label_x + text.width(self.labels[r], self.font) + 10 - self._reach(r)[0] for r in rows]
		self.center = max(needed, default=self.label_x + 10)
		right = max((self._reach(r)[1] for r in rows), default=0)
		self.genes_right = self.center + right
		self.row_pitch = self.row_h
		self.extra_right = 0.0
		self.rows_height = self.pad + len(rows) * self.row_pitch
		self._clips: dict[tuple[str, str], str] = {}

	@property
	def width(self) -> int:
		return int(self.genes_right + self.extra_right + self.pad)

	def _anchor(self, row: str) -> int:
		query = next((g for g in self.data.genes[row] if g.offset == 0), self.data.genes[row][0])
		return query.end if self.reversed[row] else query.start

	def _reach(self, row: str) -> tuple[float, float]:
		genes = self.data.genes[row]
		lo, hi = min(g.start for g in genes), max(g.end for g in genes)
		if self.reversed[row]:
			return (self.anchor[row] - hi) / self.bp_per_px, (self.anchor[row] - lo) / self.bp_per_px
		return (lo - self.anchor[row]) / self.bp_per_px, (hi - self.anchor[row]) / self.bp_per_px

	def y(self, row: str) -> float:
		return self.pad + self.rows.index(row) * self.row_pitch + self.row_h / 2

	def band_box(self, row: str) -> tuple[float, float]:
		return self.y(row) - self.gene_h / 2 - 3, self.gene_h + 6

	def x(self, row: str, bp: int) -> float:
		if self.reversed[row]:
			return self.center + (self.anchor[row] - bp) / self.bp_per_px
		return self.center + (bp - self.anchor[row]) / self.bp_per_px

	def span(self, row: str) -> tuple[float, float]:
		genes = self.data.genes[row]
		xs = [self.x(row, g.start) for g in genes] + [self.x(row, g.end) for g in genes]
		return min(xs), max(xs)

	def drawn_strand(self, row: str, gene: Gene) -> str:
		if not self.reversed[row]:
			return gene.strand
		return "-" if gene.strand == "+" else "+"

	def arrow(self, row: str, gene: Gene) -> Arrow:
		xa, xb = self.x(row, gene.start), self.x(row, gene.end)
		x0, x1 = min(xa, xb), max(xa, xb)
		cy = self.y(row)
		strand = self.drawn_strand(row, gene)
		h = self.gene_h
		if self.spec.classic:
			length = max(x1 - x0, self.spec.value("min_gene_width"))
			head = min(self.spec.value("arrow_head"), length * 0.35)
		else:
			length = max(x1 - x0, 6)
			head = min(h, length * 0.5)
		top, bot = cy - h / 2, cy + h / 2
		if strand == "-":
			pts = [(x0, cy), (x0 + head, top), (x0 + length, top), (x0 + length, bot), (x0 + head, bot)]
		else:
			pts = [(x0, top), (x0 + length - head, top), (x0 + length, cy), (x0 + length - head, bot), (x0, bot)]
		key = (row, gene.accession)
		if key not in self._clips:
			self._clips[key] = "clip{}".format(len(self._clips))
		return Arrow(x0, x0 + length, cy, strand, pts, self._clips[key])

	def x_in_gene(self, arrow: Arrow, gene: Gene, residue: int) -> float:
		length = max((gene.end - gene.start) // 3, 1)
		frac = min(max(residue / length, 0.0), 1.0)
		return arrow.x1 - frac * arrow.width if arrow.strand == "-" else arrow.x0 + frac * arrow.width

	def genes_in(self, row: str) -> list[Gene]:
		return self.data.genes[row]

	def clip_defs(self) -> str:
		parts = []
		for row in self.rows:
			for gene in self.data.genes[row]:
				arrow = self.arrow(row, gene)
				parts.append('<clipPath id="{}"><polygon points="{}"/></clipPath>'.format(arrow.clip, arrow.polygon()))
		return "<defs>" + "".join(parts) + "</defs>"
