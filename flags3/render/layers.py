from typing import Optional

from flags3.render import palettes, text
from flags3.render.layout import Layout
from flags3.schema import MISSING, Annotation, split_bp_subject

PASTEL = 0.62


def _by_subject(annotations: list[Annotation]) -> dict[str, list[Annotation]]:
	out: dict[str, list[Annotation]] = {}
	for a in annotations:
		out.setdefault(a.subject, []).append(a)
	return out


class GeneStyle:
	def __init__(self, layout: Layout, fills: dict[str, str], has_overlay: set[str], classic: bool):
		self.layout = layout
		self.fills = fills
		self.has_overlay = has_overlay
		self.classic = classic

	def accent(self, gene) -> Optional[str]:
		if gene.accession == "pseudogene*":
			return palettes.PSEUDO[1]
		if gene.is_rna or (gene.accession.endswith("*") and "rna" in gene.accession.lower()):
			return palettes.RNA[1]
		if gene.offset == 0:
			return palettes.QUERY
		return None

	def fill(self, gene) -> str:
		if self.classic and gene.offset == 0:
			return "#000000"
		if gene.accession == "pseudogene*":
			return palettes.PSEUDO[0]
		default = palettes.RNA[0] if gene.is_rna else (palettes.WHITE if self.classic else palettes.GREY)
		colour = self.fills.get(gene.accession, default)
		if gene.accession in self.has_overlay:
			return colour if palettes.is_pale(colour) else palettes.lighten(colour, PASTEL)
		return colour

	def outline(self, gene, fill: str) -> str:
		accent = self.accent(gene)
		if accent:
			return accent
		return palettes.MIDGREY if fill in (palettes.GREY, palettes.WHITE) else fill


def draw_genes(layout: Layout, style: GeneStyle) -> str:
	parts = []
	for row in layout.rows:
		for gene in layout.genes_in(row):
			arrow = layout.arrow(row, gene)
			fill = style.fill(gene)
			parts.append('<polygon points="{}" fill="{}" stroke="{}" stroke-width="2"/>'.format(
				arrow.polygon(), fill, style.outline(gene, fill)))
	return '<g id="layer-genes">{}</g>'.format("".join(parts))


HATCH_ANGLE = {"sismis": 45, "genomad": -45, "defence": 0}
HATCH_STROKE = 5
HATCH_GAP = 5


def hatch_id(stage: str, category_index: int) -> str:
	return "hatch-{}-{}".format(stage, category_index)


def hatch_defs(stage: str, colours: dict[str, str]) -> str:
	angle = HATCH_ANGLE.get(stage, 45)
	period = HATCH_STROKE + HATCH_GAP
	parts = []
	for i, colour in enumerate(colours.values()):
		parts.append('<pattern id="{}" patternUnits="userSpaceOnUse" width="{}" height="{}" patternTransform="rotate({})">'
			'<line x1="{}" y1="0" x2="{}" y2="{}" stroke="{}" stroke-width="{}"/></pattern>'.format(
				hatch_id(stage, i), period, period, angle, HATCH_STROKE / 2, HATCH_STROKE / 2, period, colour, HATCH_STROKE))
	return "".join(parts)


def draw_bands(layout: Layout, stage: str, annotations: list[Annotation], colours: dict[str, str],
		codes: dict[str, str], opacity: float) -> tuple[str, dict[str, list], set[str]]:
	by_subject = _by_subject(annotations)
	index = {c: i for i, c in enumerate(colours)}
	parts, side, drawn = [hatch_defs(stage, colours)], {}, set()
	for row in layout.rows:
		window = layout.data.windows[row]
		subject = "{}|{}".format(window.assembly, window.contig)
		hits = by_subject.get(subject, [])
		if not hits:
			continue
		row_lo, row_hi = layout.span(row)
		top, height = layout.band_box(row)
		for a in hits:
			xa, xb = layout.x(row, a.start), layout.x(row, a.end)
			x0, x1 = max(min(xa, xb), row_lo - 4), min(max(xa, xb), row_hi + 4)
			if x1 <= x0:
				continue
			colour = colours[a.category]
			drawn.add(a.category)
			parts.append('<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="url(#{})" fill-opacity="{}"/>'.format(
				x0, top, x1 - x0, height, hatch_id(stage, index[a.category]), opacity))
			code = codes.get(a.category, "")
			if code and (code, colour) not in side.setdefault(row, []):
				side[row].append((code, colour))
	return '<g id="layer-{}">{}</g>'.format(stage, "".join(parts)), side, drawn


def draw_side_codes(layout: Layout, per_row: dict[str, list]) -> str:
	parts = []
	size = layout.font - 2
	for row, items in per_row.items():
		x = layout.span(row)[1] + 8
		y = layout.y(row) + size / 3
		for code, colour in items:
			parts.append(text.text(x, y, code, size, None, colour, "bold"))
			x += text.width(code, size) + 6
	return '<g id="layer-band-codes">{}</g>'.format("".join(parts))


def side_width(layout: Layout, per_row: dict[str, list]) -> float:
	size = layout.font - 2
	return max((sum(text.width(code, size) + 6 for code, _ in items) + 8 + layout.span(row)[1] - layout.genes_right
		for row, items in per_row.items()), default=0.0)


def draw_wedges(layout: Layout, stage: str, annotations: list[Annotation], colours: dict[str, str],
		codes: dict[str, str], domain_h: float) -> tuple[str, dict[str, list]]:
	by_subject = _by_subject(annotations)
	parts, labels = [], {}
	u = domain_h / 4.0
	for row in layout.rows:
		for gene in layout.genes_in(row):
			hits = by_subject.get(gene.accession)
			if not hits:
				continue
			arrow = layout.arrow(row, gene)
			cy = arrow.cy
			inner = []
			for a in sorted(hits, key=lambda a: (a.start, a.end)):
				s, e = layout.x_in_gene(arrow, gene, a.start), layout.x_in_gene(arrow, gene, a.end)
				colour = colours[a.category]
				pts = [(s, cy + 2 * u), (s, cy + u), (e, cy - 2 * u), (e, cy + 2 * u)]
				inner.append('<polygon points="{}" fill="{}"/>'.format(
					" ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts), colour))
				labels.setdefault(row, []).append(((s + e) / 2, codes.get(a.category, ""), colour))
			parts.append('<g clip-path="url(#{})">{}</g>'.format(arrow.clip, "".join(inner)))
	return '<g id="layer-{}">{}</g>'.format(stage, "".join(parts)), labels


def draw_features(layout: Layout, stage: str, annotations: list[Annotation]) -> str:
	by_subject = _by_subject(annotations)
	parts, on_top = [], []
	h = layout.gene_h
	for row in layout.rows:
		for gene in layout.genes_in(row):
			hits = by_subject.get(gene.accession)
			if not hits:
				continue
			arrow = layout.arrow(row, gene)
			cy = arrow.cy
			inner = []
			for a in hits:
				xs, xe = layout.x_in_gene(arrow, gene, a.start), layout.x_in_gene(arrow, gene, a.end)
				lo, hi = min(xs, xe), max(xs, xe)
				if a.kind == "segment":
					top, bot = cy - h / 2, cy + h / 2
					for k in range(1, 5):
						hy = top + (bot - top) * k / 5
						inner.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="#d40000" stroke-width="0.8"/>'.format(lo, hy, hi, hy))
					inner.append('<line x1="{0:.1f}" y1="{1:.1f}" x2="{0:.1f}" y2="{2:.1f}" stroke="#fff" stroke-width="1.2"/>'.format((lo + hi) / 2, top, bot))
				elif a.kind == "triangle":
					tip = lo if arrow.strand == "+" else hi
					w = max(hi - lo, 3) * 0.6
					pts = [(tip, cy - h / 2), (tip - w / 2, cy + h / 2), (tip + w / 2, cy + h / 2)]
					on_top.append('<polygon points="{}" fill="#000"/>'.format(" ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)))
			if inner:
				parts.append('<g clip-path="url(#{})">{}</g>'.format(arrow.clip, "".join(inner)))
	return '<g id="layer-{}">{}{}</g>'.format(stage, "".join(parts), "".join(on_top))


def draw_outlines(layout: Layout, style: GeneStyle, overlaid: set[str]) -> str:
	parts = []
	for row in layout.rows:
		for gene in layout.genes_in(row):
			if gene.accession not in overlaid:
				continue
			arrow = layout.arrow(row, gene)
			parts.append('<polygon points="{}" fill="none" stroke="{}" stroke-width="2"/>'.format(
				arrow.polygon(), style.outline(gene, style.fill(gene))))
	return '<g id="layer-outlines">{}</g>'.format("".join(parts))


def draw_row_labels(layout: Layout) -> str:
	parts = []
	for row in layout.rows:
		parts.append(text.text(layout.pad + layout.tree_w, layout.y(row) + layout.font / 3, layout.labels[row]))
	return '<g id="layer-row-labels">{}</g>'.format("".join(parts))


def place_labels(layout: Layout, per_row: dict[str, list], step: int) -> str:
	parts = []
	size = layout.font - 4
	for row, items in per_row.items():
		y = layout.y(row) - layout.gene_h / 2 - 4
		right = None
		for x, label, colour in sorted(items, key=lambda t: t[0]):
			if not label:
				continue
			half = text.width(label, size) / 2 + 2
			if right is not None and x - half < right:
				x = right + half
			right = x + half
			parts.append(text.text(x, y, label, size, "middle", colour))
	return '<g id="layer-labels">{}</g>'.format("".join(parts))


def draw_classic_numbers(layout: Layout, numbers: dict[str, str], fills: dict[str, str], no_overlaps: bool) -> str:
	parts = []
	size = min(layout.font - 1, layout.gene_h - 4)
	for row in layout.rows:
		for gene in layout.genes_in(row):
			label = numbers.get(gene.accession)
			if not label:
				continue
			arrow = layout.arrow(row, gene)
			if no_overlaps and text.width(label, size) > arrow.width - 2:
				continue
			fill = "#000000" if gene.offset == 0 else fills.get(gene.accession, palettes.WHITE)
			colour = palettes.readable_on(fill) if gene.offset == 0 else "#000"
			parts.append(text.text((arrow.x0 + arrow.x1) / 2, arrow.cy + size / 3, label, size, "middle", colour))
	return '<g id="layer-labels">{}</g>'.format("".join(parts))
