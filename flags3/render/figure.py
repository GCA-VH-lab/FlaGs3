import math
from io import StringIO
from typing import Optional

from flags3.render import layers, palettes, text
from flags3.render.data import RunData
from flags3.render.layout import Layout
from flags3.render.style import BAND_CODES, BAND_ORDER, STAGE_TITLES, Colours, FigureSpec

CLUSTER_STAGES = ("cluster", "cluster_rna")


class Figure:
	def __init__(self, data: RunData, spec: FigureSpec, colours: Colours, rows: Optional[list[str]] = None,
			no_overlaps: bool = False):
		self.data = data
		self.spec = spec
		self.colours = colours
		self.no_overlaps = no_overlaps
		self.active = {s: a for s, a in data.annotations.items() if spec.wants(s)}
		order = data.tree_order if (spec.tree and data.tree_order) else data.order
		self.rows = [r for r in (rows or order) if r in data.genes]
		self.layout = Layout(data, spec, self.rows, [s for s in self.active if s in BAND_ORDER])
		self.legend_panels: list[tuple[str, list[tuple[str, str, str]], Optional[str]]] = []

	def has_content(self) -> bool:
		return bool(self.rows) and any(s not in CLUSTER_STAGES for s in self.active) or (
			bool(self.rows) and any(s in self.active for s in CLUSTER_STAGES)) or (bool(self.rows) and not self.spec.layers)

	def render(self) -> str:
		L = self.layout
		fills, numbers = self._families()
		domain_stage = next((s for s in self.active if s == "domains"), None)
		feature_stage = next((s for s in self.active if s == "features"), None)
		overlaid = set()
		for s in (domain_stage, feature_stage):
			if s:
				overlaid |= {a.subject for a in self.active[s]}
		style = layers.GeneStyle(L, fills, overlaid, self.spec.classic)
		body = [L.clip_defs()]
		side: dict[str, list] = {}
		for stage in L.band_stages:
			table = self.colours.assign(stage, self.active[stage])
			codes = {c: "{}{}".format(BAND_CODES[stage], i + 1) for i, c in enumerate(table)}
			svg, per_row, drawn = layers.draw_bands(L, stage, self.active[stage], table, codes, self.spec.value("band_opacity"))
			body.append(svg)
			for row, items in per_row.items():
				side.setdefault(row, []).extend(items)
			self.legend_panels.append((STAGE_TITLES[stage], [(codes[c], self._band_label(stage, c), table[c]) for c in table if c in drawn],
				layers.hatch_id(stage, 0).rsplit("-", 1)[0]))
		L.extra_right = layers.side_width(L, side)
		body.append(layers.draw_genes(L, style))
		above: dict[str, list] = {}
		if domain_stage:
			table = self.colours.assign(domain_stage, self.active[domain_stage])
			codes = {c: str(i + 1) for i, c in enumerate(table)}
			svg, labels = layers.draw_wedges(L, domain_stage, self.active[domain_stage], table, codes, self.spec.value("domain_height"))
			body.append(svg)
			for row, items in labels.items():
				above.setdefault(row, []).extend(items)
			self.legend_panels.append((STAGE_TITLES["domains"], [(codes[c], self._domain_label(c), table[c]) for c in table], None))
		if overlaid:
			body.append(layers.draw_outlines(L, style, overlaid))
		if feature_stage:
			body.append(layers.draw_features(L, feature_stage, self.active[feature_stage]))
		if side:
			body.append(layers.draw_side_codes(L, side))
		if self.spec.numbers and numbers:
			if self.spec.classic:
				body.append(layers.draw_classic_numbers(L, numbers, fills, self.no_overlaps))
			else:
				for row in self.rows:
					for gene in L.genes_in(row):
						label = numbers.get(gene.accession)
						if label:
							arrow = L.arrow(row, gene)
							if domain_stage and label[:1].isdigit():
								label = "G" + label
							above.setdefault(row, []).append(((arrow.x0 + arrow.x1) / 2, label, "#000"))
		if above:
			body.append(layers.place_labels(L, above, self.spec.value("label_step")))
		body.append(layers.draw_row_labels(L))
		tree_svg, scale_h = self._tree()
		if tree_svg:
			body.insert(1, tree_svg)
		legend_svg, legend_h = self._legend(L.rows_height + 14 + scale_h, L.width, feature_stage is not None)
		height = int(L.rows_height + scale_h + legend_h + L.pad)
		head = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" font-family="{}" font-size="{}">'.format(
			L.width, height, text.FONT_FAMILY, L.font), '<rect width="{}" height="{}" fill="white"/>'.format(L.width, height)]
		return "\n".join(head + body + [legend_svg, "</svg>"])

	def _families(self) -> tuple[dict[str, str], dict[str, str]]:
		fills, numbers = {}, {}
		for stage in CLUSTER_STAGES:
			if stage not in self.active:
				continue
			table = self.colours.assign(stage, self.active[stage])
			for a in self.active[stage]:
				fills[a.subject] = table[a.category]
				numbers[a.subject] = a.label
		return fills, numbers

	def _band_label(self, stage: str, category: str) -> str:
		labels = {a.label for a in self.active[stage] if a.category == category}
		return category if len(labels) != 1 else labels.pop()

	def _domain_label(self, category: str) -> str:
		labels = []
		for a in self.active["domains"]:
			if a.category == category and a.label not in labels:
				labels.append(a.label)
		return category if category in labels or len(labels) > 3 else "{} ({})".format(category, ", ".join(labels))

	def _tree(self) -> tuple[str, float]:
		L = self.layout
		if not L.tree_w:
			return "", 0
		from Bio import Phylo
		tree = Phylo.read(StringIO(self.data.newick), "newick")
		try:
			tree.root_at_midpoint()
		except Exception:
			pass
		tree.ladderize()
		depths = tree.depths()
		if not any(depths.values()):
			depths = tree.depths(unit_branch_lengths=True)
		maxd = max(depths.values()) or 1
		xscale = (L.tree_w - 10) / maxd
		x0 = L.pad
		y_of = {r: L.y(r) for r in self.rows}
		yc = {}

		def assign(clade):
			if clade.is_terminal():
				yc[clade] = y_of.get(clade.name, L.pad)
			else:
				ys = [assign(c) for c in clade.clades]
				yc[clade] = sum(ys) / len(ys)
			return yc[clade]

		assign(tree.root)
		seg = []

		def walk(clade, px):
			x = x0 + depths[clade] * xscale
			y = yc[clade]
			seg.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="#555"/>'.format(px, y, x, y))
			if clade.is_terminal():
				seg.append('<circle cx="{:.1f}" cy="{:.1f}" r="2.2" fill="#555"/>'.format(x, y))
			else:
				cys = [yc[c] for c in clade.clades]
				seg.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="#555"/>'.format(x, min(cys), x, max(cys)))
				if clade.confidence is not None:
					seg.append(text.text(x - 4, y + 9, "{:g}".format(clade.confidence), 8, "end", "#8b0000"))
				for c in clade.clades:
					walk(c, x)

		walk(tree.root, x0)
		target = maxd / 5
		mag = 10 ** math.floor(math.log10(target))
		nice = min((1, 2, 5, 10), key=lambda m: abs(m * mag - target)) * mag
		w = nice * xscale
		y = L.rows_height + 16
		seg.append('<line x1="{0:.1f}" y1="{1}" x2="{2:.1f}" y2="{1}" stroke="#333" stroke-width="1.5"/>'.format(x0, y, x0 + w))
		seg.append(text.text(x0 + w / 2, y + L.font + 2, "{:g}".format(nice), L.font - 1, "middle"))
		return '<g id="layer-tree">{}</g>'.format("".join(seg)), 34

	def _legend(self, y: float, width: float, features: bool) -> tuple[str, float]:
		L = self.layout
		panels = list(self.legend_panels)
		genes = []
		present = {(g.offset == 0, g.is_rna, g.accession == "pseudogene*") for r in self.rows for g in L.genes_in(r)}
		if any(q for q, _, _ in present):
			genes.append(("", "Query protein", "outline"))
		if any(r for _, r, _ in present):
			genes.append(("", "RNA gene", "rna"))
		if any(p for _, _, p in present):
			genes.append(("", "Pseudogene", "pseudo"))
		if features:
			kinds = {a.kind for a in self.active["features"]}
			if "segment" in kinds:
				genes.append(("", "Transmembrane region", "hatch"))
			if "triangle" in kinds:
				genes.append(("", "Signal peptide", "triangle"))
		parts, top = [], y
		for title, items, hatch in panels:
			longest = max((text.width("{}. {}".format(code, label), L.font) for code, label, _ in items), default=0)
			col_w = int(longest) + 26
			cols = max(1, min(len(items), (width - 2 * L.pad) // col_w))
			parts.append(text.text(L.pad, y, title, weight="bold"))
			for i, (code, label, colour) in enumerate(items):
				ex = L.pad + (i % cols) * col_w
				ey = y + 16 + (i // cols) * 16
				fill = "url(#{}-{})".format(hatch, self._hatch_index(hatch, colour)) if hatch else colour
				parts.append('<rect x="{}" y="{}" width="12" height="10" fill="{}" stroke="{}" stroke-width="1"/>'.format(ex + 2, ey - 8, fill, colour))
				parts.append(text.text(ex + 20, ey, "{}. {}".format(code, label)))
			y += (-(-len(items) // cols) + 1) * 16 + 8
		if genes:
			parts.append(text.text(L.pad, y, "Genes", weight="bold"))
			ex, ey = L.pad, y + 16
			for _, label, shape in genes:
				parts.append(self._mark(ex + 2, ey - 8, shape))
				parts.append(text.text(ex + 20, ey, label))
				ex += int(text.width(label, L.font)) + 46
			y += 36
		return '<g id="layer-legend">{}</g>'.format("".join(parts)), y - top

	def _hatch_index(self, hatch: str, colour: str) -> int:
		stage = hatch.split("-", 1)[1]
		table = self.colours.assigned.get(stage, {})
		return list(table.values()).index(colour) if colour in table.values() else 0

	@staticmethod
	def _mark(x, y, shape) -> str:
		if shape == "hatch":
			return "".join('<line x1="{}" y1="{:.1f}" x2="{}" y2="{:.1f}" stroke="#d40000" stroke-width="0.8"/>'.format(
				x, y + 2 + k * 3, x + 12, y + 2 + k * 3) for k in range(3))
		if shape == "triangle":
			return '<polygon points="{},{} {},{} {},{}" fill="#000"/>'.format(x + 6, y, x, y + 10, x + 12, y + 10)
		stroke = {"outline": "#000", "rna": palettes.RNA[1], "pseudo": palettes.PSEUDO[1]}[shape]
		fill = {"outline": palettes.WHITE, "rna": palettes.RNA[0], "pseudo": palettes.PSEUDO[0]}[shape]
		return '<rect x="{}" y="{}" width="12" height="10" fill="{}" stroke="{}" stroke-width="2"/>'.format(x, y, fill, stroke)


def parts_for(data: RunData, spec: FigureSpec, colours: Colours, max_height: int, no_overlaps: bool) -> list[Figure]:
	whole = Figure(data, spec, colours, no_overlaps=no_overlaps)
	if not whole.rows:
		return []
	if whole.layout.tree_w or len(whole.rows) < 2:
		return [whole]
	per_row = whole.layout.row_pitch
	capacity = max(1, (max_height - 200) // per_row)
	if len(whole.rows) <= capacity:
		return [whole]
	rows = whole.rows
	return [Figure(data, spec, Colours(spec.palette, colours.overrides, spec.monochrome), rows[i:i + capacity], no_overlaps)
		for i in range(0, len(rows), capacity)]
