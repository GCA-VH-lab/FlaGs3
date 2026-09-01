import colorsys
from io import StringIO
from typing import Dict, List, Optional

from Bio import Phylo


class _FlaGsBase: 
	GREY = "#d9d9d9"
	PSEUDO = ("#f2f2f3", "#000080")   
	RNA = ("#f2f2f2", "#008000")      
	OTHER = ("#ffffff", "#bebebe")    

	FONT_FAMILY = "Arial, 'Liberation Sans', Helvetica, sans-serif"

	_CHAR_WIDTH = {
		' ':278,'!':278,'"':355,'#':556,'$':556,'%':889,'&':667,"'":191,
		'(':333,')':333,'*':389,'+':584,',':278,'-':333,'.':278,'/':278,
		'0':556,'1':556,'2':556,'3':556,'4':556,'5':556,'6':556,'7':556,'8':556,'9':556,
		':':278,';':278,'<':584,'=':584,'>':584,'?':556,'@':1015,
		'A':667,'B':667,'C':722,'D':722,'E':667,'F':611,'G':778,'H':722,'I':278,
		'J':500,'K':667,'L':556,'M':833,'N':722,'O':778,'P':667,'Q':778,'R':722,
		'S':667,'T':611,'U':722,'V':667,'W':944,'X':667,'Y':667,'Z':611,
		'[':278,'\\':278,']':278,'^':469,'_':556,'`':333,
		'a':556,'b':556,'c':500,'d':556,'e':556,'f':278,'g':556,'h':556,'i':222,
		'j':222,'k':500,'l':222,'m':833,'n':556,'o':556,'p':556,'q':556,'r':333,
		's':500,'t':278,'u':556,'v':500,'w':722,'x':500,'y':500,'z':500,
		'{':334,'|':260,'}':334,'~':584,
	}
	_CHAR_WIDTH_DEFAULT = 556

	@classmethod
	def _text_width(cls, text: str, font_size: float) -> float:
		units = sum(cls._CHAR_WIDTH.get(ch, cls._CHAR_WIDTH_DEFAULT) for ch in text)
		return units / 1000.0 * font_size

	@staticmethod
	def _special_type(accession: str) -> Optional[str]:
		if not accession.endswith("*"):
			return None
		low = accession.lower()
		if low.startswith("pseudo") or low.startswith("ps"):
			return "pseudo"
		if "rna" in low:
			return "rna"
		return "other"

	def _style_for(self, accession: str, color: Dict[str, str], is_rna: bool = False):
		special = self._special_type(accession)
		if special == "pseudo":
			return self.PSEUDO
		if special == "rna":
			return self.RNA
		if special == "other":
			return self.OTHER
		if is_rna:
			return color.get(accession, self.RNA[0]), self.RNA[1]
		return color.get(accession, self.GREY), "#333"

	@staticmethod
	def _palette(n: int) -> List[str]:
		colors = []
		for i in range(n):
			r, g, b = colorsys.hsv_to_rgb(i / n if n else 0, 0.55, 0.85)
			colors.append("#{:02x}{:02x}{:02x}".format(
				int(r * 255), int(g * 255), int(b * 255)))
		return colors

	@staticmethod
	def _classic_palette(n: int) -> List[str]:
		"""FlaGs2's random_color(): 20 hues on a 5-step grid, L=0.5, S=0.5."""
		hues = [int(h * 3.6) / 100.0 for h in range(0, 100, 5)]
		return ["#%02x%02x%02x" % tuple(int(f * 255) for f in
									    colorsys.hls_to_rgb(hues[i % len(hues)], 0.5, 0.5))
				for i in range(n)]

	def _family_colors(self, families: List[List[str]]) -> Dict[str, str]:
		if getattr(self, "monochrome", False):
			return {acc: self.GREY for fam in families for acc in fam}
		numbers = getattr(self, "numbers", None) or {}
		multi = ([fam for fam in families if any(a in numbers for a in fam)]
				 if numbers else [fam for fam in families if len(fam) > 1])
		palette = (self._classic_palette(len(multi))
				   if getattr(self, "classic", False) else self._palette(len(multi)))
		color = {}
		for fam, c in zip(multi, palette):
			for acc in fam:
				color[acc] = c
		for fam in families:
			for acc in fam:
				color.setdefault(acc, self.GREY)
		return color

	@staticmethod
	def _escape(text: str) -> str:
		return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


	def _draw_tree(self, tree, y_of: Dict[str, float]) -> str:
	  depths = tree.depths()
	  if not any(depths.values()):
	    depths = tree.depths(unit_branch_lengths=True)
	  maxd = max(depths.values()) or 1
	  xscale = (self.tree_w - 10) / maxd
	  x0 = self.pad
	  yc: Dict = {}

	  def assign(clade):
	    if clade.is_terminal():
	      yc[clade] = y_of.get(clade.name, self.pad)
	    else:
	      ys = [assign(c) for c in clade.clades]
	      yc[clade] = sum(ys) / len(ys)
	    return yc[clade]
	  assign(tree.root)

	  seg = []
	  def walk(clade, px):
	    x = x0 + depths[clade] * xscale
	    y = yc[clade]
	    seg.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#555"/>'.format(px, y, x, y))
	    if clade.is_terminal():
	      seg.append('<circle cx="{}" cy="{}" r="2.2" fill="#555"/>'.format(x, y))
	    else:
	      cys = [yc[c] for c in clade.clades]
	      seg.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#555"/>'.format(
	        x, min(cys), x, max(cys)))
	      if clade.confidence is not None:
	        seg.append('<text x="{}" y="{}" font-size="8" fill="#8b0000" '
	          'text-anchor="end">{:g}</text>'.format(x - 4, y + 9, clade.confidence))
	      for c in clade.clades:
	        walk(c, x)
	  walk(tree.root, x0)
	  return "\n".join(seg), xscale, maxd


	def _scale_bar(self, x: float, y: float, xscale: float, maxd: float) -> str:
	  import math
	  target = maxd / 5 if maxd else 0.1
	  if target <= 0:
	    return ""
	  mag = 10 ** math.floor(math.log10(target))
	  nice = min((1, 2, 5, 10), key=lambda m: abs(m * mag - target)) * mag
	  w = nice * xscale
	  return ('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#333" stroke-width="1.5"/>'
	      '<line x1="{0}" y1="{3}" x2="{0}" y2="{4}" stroke="#333"/>'
	      '<line x1="{2}" y1="{3}" x2="{2}" y2="{4}" stroke="#333"/>'
	      '<text x="{5}" y="{6}" text-anchor="middle" font-size="{7}">{8:g}</text>'
	      ).format(x, y, x + w, y - 3, y + 3, x + w / 2, y + self.font + 2,
	        self.font - 1, nice)

	@property
	def classic(self) -> bool:
		return getattr(self, "style", "") == "classic"

	def _gene_style(self, gene, gene_color):
		if self.classic and gene.offset == 0:
			return "#000000", "#000000"
		special = self._special_type(gene.accession)
		default = {"pseudo": self.PSEUDO[0], "rna": self.RNA[0],
				   "other": self.OTHER[0]}.get(
					   special, self.RNA[0] if gene.is_rna else self.GREY)
		fill = gene_color.get(gene.accession, default)
		if self.classic and fill == self.GREY:
			fill = "#ffffff"
		pale = (self.GREY, "#ffffff", self.PSEUDO[0], self.RNA[0], self.OTHER[0])
		outline = self.MIDGREY if fill in pale else fill
		return fill, (self._accent(gene) or outline)

	QUERY_ACCENT = "#000000"
	PASTEL_FILL = 0.62

	def _is_pale(self, fill):
		return fill in (self.GREY, "#ffffff", self.PSEUDO[0], self.RNA[0],
						self.OTHER[0])

	def _pastel_outline(self, fill):
		return self.MIDGREY if self._is_pale(fill) else self._pastel(fill)

	def _pastel(self, fill):
		return fill if self._is_pale(fill) else self._lighten(fill, self.PASTEL_FILL)

	@staticmethod
	def _lighten(colour, amount):
		"""Blend a hex colour towards white. amount 0 = unchanged, 1 = white."""
		try:
			r, g, b = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
		except (ValueError, IndexError):
			return colour
		mix = lambda c: int(round(c + (255 - c) * amount))
		return "#{:02x}{:02x}{:02x}".format(mix(r), mix(g), mix(b))

	def _stroke_width(self, gene):
		if self._special_type(gene.accession) in ("rna", "pseudo") or gene.is_rna:
			return 2
		return 2 if gene.offset == 0 else 1

	def _accent(self, gene):
		"""Outer ring colour for a gene that is more than just its family."""
		special = self._special_type(gene.accession)
		if special == "pseudo":
			return self.PSEUDO[1]
		if special == "rna" or gene.is_rna:
			return self.RNA[1]
		if gene.offset == 0:
			return self.QUERY_ACCENT
		return None

	def _numbers_for(self, families, rna_accessions):
		return dict(getattr(self, "numbers", None) or {})


class OperonView(_FlaGsBase):
	LABEL_STEP = 12   
	SECRETION_BAND_OPACITY = 0.28   
	SECRETION_BAND_PAD = 2          
	SECRETION_LABEL_SCALE = 0.85

	CLASSIC = {"row_h": 20, "gene_h": 15, "font": 12, "head": 7, "min_w": 13}

	def __init__(self, mode: str = "families", row_h: int = 26, gene_h: int = 8,
				 bp_per_px: float = 10.4, pad: int = 16, font: int = 13,
				 domain_h: int = 6, features=frozenset(), monochrome: bool = False,
				 show_numbers: bool = True, style: str = "versatile"):
		self.style = style
		if style == "classic":
			row_h = self.CLASSIC["row_h"] if row_h == 26 else row_h
			gene_h = self.CLASSIC["gene_h"] if gene_h == 8 else gene_h
			font = self.CLASSIC["font"] if font == 13 else font
		self.mode = mode
		self.row_h = row_h
		self.gene_h = gene_h
		self.bp_per_px = bp_per_px
		self.pad = pad
		self.font = font
		self.domain_h = domain_h
		self.features = frozenset(features)
		self.monochrome = monochrome
		self.show_numbers = show_numbers
		self.tree_w = 0          # set per figure; 0 means no tree gutter
		self.newick = ""

	def _row_spans(self, rows, by_query):
		spans = {}
		row_reversed = {}
		max_left = max_right = 0
		for q in rows:
			genes = by_query[q]
			qg = next((g for g in genes if g.offset == 0), genes[len(genes) // 2])
			spans[q] = (qg.start + qg.end) / 2
			reversed_row = self._row_reversed(genes)
			row_reversed[q] = reversed_row
			if reversed_row:
				left_reach = (spans[q] - max(g.end for g in genes)) / self.bp_per_px
				right_reach = (spans[q] - min(g.start for g in genes)) / self.bp_per_px
			else:
				left_reach = (min(g.start for g in genes) - spans[q]) / self.bp_per_px
				right_reach = (max(g.end for g in genes) - spans[q]) / self.bp_per_px
			max_left = min(max_left, left_reach)
			max_right = max(max_right, right_reach)
		return spans, row_reversed, max_left, max_right

	def _drawn_elements(self, by_query, features):
		drawn = set()
		for gs in by_query.values():
				for gene in gs:
					if gene.offset == 0:
						drawn.add("query")
					if gene.is_rna or self._special_type(gene.accession) == "rna":
						drawn.add("rna")
					if self._special_type(gene.accession) == "pseudo":
						drawn.add("pseudo")
					for kind, _, _ in (features or {}).get(gene.accession, []):
						drawn.add("tm" if kind == "tm" else
								  "signal" if kind == "signal" else kind)
		return drawn

	def _tree_gutter(self, rows, gutter):
		out = []
		y_of = {q: self.pad + i * self.row_h + self.row_h / 2
				for i, q in enumerate(rows)}
		try:
			tree = Phylo.read(StringIO(self.newick), "newick")
			try:
				tree.root_at_midpoint()
			except Exception:
				pass
			tree.ladderize()
			saved, self.tree_w = self.tree_w, gutter
			tree_svg, xscale, maxd = self._draw_tree(tree, y_of)
			out.append(tree_svg)
			self.tree_w = saved
			if xscale and maxd:
				out.append(self._scale_bar(
					self.pad, self.pad + len(rows) * self.row_h + 16,
					xscale, maxd))
		except Exception:
			pass
		return out

	def _draw_row(self, q, y, rows, by_query, spans, row_reversed, center,
				  gutter, labels_out, overlay, gene_color, domains, features):
		svg = []
		q_mid = spans[q]
		reversed_row = row_reversed[q]

		svg.append('<text x="{}" y="{}">{}</text>'.format(
			self.pad + gutter, y + self.font / 3, self._escape(labels_out[q])))
		row_labels = []

		if self._sec:
			gxs = []
			for g in by_query[q]:
				gxs.append(self._x_for(g.start, center, q_mid, reversed_row))
				gxs.append(self._x_for(g.end, center, q_mid, reversed_row))
			row_lo, row_hi = (min(gxs), max(gxs)) if gxs else (0, 0)
			for h in self._sec["row_hits"].get(q, []):
				band, mid = self._secretion_band(
					h, self._sec["color"][h.type], center, q_mid, reversed_row,
					y, row_lo, row_hi)
				if band:
					svg.append(band)
					self._sec_labels.setdefault(q, []).append(
						(h.type, self._sec["color"][h.type]))

		for g in by_query[q]:
			if reversed_row:
				gx0 = center + (q_mid - g.end) / self.bp_per_px
				gx1 = center + (q_mid - g.start) / self.bp_per_px
				drawn = "+" if g.strand == "-" else "-"
			else:
				gx0 = center + (g.start - q_mid) / self.bp_per_px
				gx1 = center + (g.end - q_mid) / self.bp_per_px
				drawn = g.strand
			drawn = g.strand if not reversed_row else drawn

			fill, outline = self._gene_style(g, gene_color)
			stroke = outline
			sw = self._stroke_width(g)

			has_overlay = bool(self._dom or features.get(g.accession))
			if has_overlay and not self.show_numbers:
				pass
			elif has_overlay:
				pastel = self._pastel(fill)
				svg.append(self._gene_fill(
					gx0, gx1, y, g.strand, pastel,
					self._accent(g) or self._pastel_outline(fill), sw))
			else:
				svg.append(self._gene_fill(gx0, gx1, y, g.strand, fill, stroke, sw))
			if has_overlay:
				wedges, wlabels = self._domains_on_gene(
					g, gx0, gx1, y, g.strand, (domains or {}).get(g.accession, []),
					self._dom or {"color": {}, "number": {}})
				feats = self._features_on_gene(
					g, gx0, gx1, y, g.strand, features.get(g.accession, []))
				overlay_svg = wedges + feats
				if overlay_svg:
					clip_id = "clip{}".format(self._clip_n); self._clip_n += 1
					clip_pts = " ".join("{:.1f},{:.1f}".format(px, py)
										for px, py in self._arrow_points(gx0, gx1, y, g.strand))
					svg.append('<clipPath id="{}"><polygon points="{}"/></clipPath>'.format(
						clip_id, clip_pts))
					svg.append('<g clip-path="url(#{})">{}</g>'.format(clip_id, overlay_svg))
				svg.append(self._gene_outline(
					gx0, gx1, y, g.strand,
					self._accent(g) or (self._pastel_outline(fill)
										if self.show_numbers else self.MIDGREY),
					sw))
				row_labels.extend(wlabels)
			num = overlay["number"].get(g.accession) if self.show_numbers else None
			if num is not None and self._dom and num[:1].isdigit():
				num = "G{}".format(num)
			if num is not None:
				if self.classic:
					svg.append(self._number_in_gene(gx0, gx1, y, num, fill))
				else:
					row_labels.append(((gx0 + gx1) / 2, num, "#000"))

		svg.append(self._place_labels(row_labels, y))
		return svg

	def render(self, neighborhoods, families, species=None, order=None,
			   domains=None, clans=None, features=None, labels=None, secretion=None):
		species = species or {}
		features = features or {}
		row_labels_map = labels or {}
		by_query = {}
		for g in neighborhoods:
			by_query.setdefault(g.query, []).append(g)
		for q in by_query:
			by_query[q].sort(key=lambda x: x.start)

		rows = [q for q in (order or list(by_query)) if q in by_query] or list(by_query)
		gene_color = self._family_colors(families)
		rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
		overlay = self._family_overlay(families, rna_accessions)
		self._dom = (self._domain_overlay(rows, by_query, domains or {}, clans or {})
					 if domains else None)
		self._sec = (self._secretion_overlay(rows, by_query, secretion or [])
					 if secretion else None)
		def _row_label(q):
			if q in row_labels_map:
				return row_labels_map[q]
			return q if not species.get(q) else "{}  {}".format(q, species[q])
		labels_out = {q: _row_label(q) for q in rows}
		label_w = int(max((self._text_width(l, self.font) for l in labels_out.values()),
						   default=0)) + 2

		spans, row_reversed, max_left, max_right = self._row_spans(rows, by_query)

		gutter = self.tree_w if (self.tree_w and self.newick) else 0
		center = self.pad + gutter + label_w - max_left
		genes_right = center + max_right 
		sec_w = 0
		if secretion:
			label_size = max(7, int(self.font * self.SECRETION_LABEL_SCALE))
			widest = max((self._text_width(h.type, label_size)
						  for h in secretion), default=0)
			sec_w = int(widest) + 14
		W = int(genes_right + sec_w + self.pad)

		drawn = self._drawn_elements(by_query, features)
		panels = []
		for title, ov in (("Domains", self._dom), ("Secretion systems", self._sec)):
			if ov:
				items, cols, col_w, h = self._legend_layout(ov, W)
				if items:
					panels.append((title, items, cols, col_w, h))
		legend_h = sum(p[4] for p in panels)
		extra_blocks = self._extra_legend(drawn)
		extra_h = sum(36 + 16 * (len(items) // 4) for _, items in extra_blocks)
		scale_h = 26 if (self.tree_w and self.newick) else 0
		H = self.pad * 2 + len(rows) * self.row_h + legend_h + extra_h + scale_h

		svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
			   'font-family="{}" font-size="{}">'.format(W, H, self.FONT_FAMILY, self.font),
			   '<rect width="{}" height="{}" fill="white"/>'.format(W, H)]

		if gutter:
			svg.extend(self._tree_gutter(rows, gutter))

		self._clip_n = 0
		self._sec_labels = {}
		for i, q in enumerate(rows):
			svg.extend(self._draw_row(
				q, self.pad + i * self.row_h + self.row_h / 2, rows, by_query,
				spans, row_reversed, center, gutter, labels_out, overlay,
				gene_color, domains, features))

		for i, q in enumerate(rows):
			seen, x = set(), genes_right + 8
			y = self.pad + i * self.row_h + self.row_h / 2
			for kind, colour in self._sec_labels.get(q, []):
				if kind in seen:
					continue
				seen.add(kind)
				size = max(7, int(self.font * self.SECRETION_LABEL_SCALE))
				svg.append('<text x="{:.1f}" y="{:.1f}" font-size="{}" '
						   'fill="#000000">{}</text>'.format(
							   x, y + size / 3, size, self._escape(kind)))
				x += self._text_width(kind, size) + 6

		ly = self.pad + len(rows) * self.row_h + 14 + scale_h
		for title, items, cols, col_w, h in panels:
			svg.append(self._legend(items, self.pad, ly, cols, col_w, title=title))
			ly += h
		if extra_blocks:
			block_svg, used = self._extra_legend_svg(extra_blocks, self.pad, ly, W)
			svg.append(block_svg)
		svg.append('</svg>')
		return "\n".join(svg)

	def _family_overlay(self, families, rna_accessions):
		return {
			"color": self._family_colors(families),
			"number": self._numbers_for(families, rna_accessions),
			"legend": None,
		}

	def _domain_overlay(self, rows, by_query, domains, clans):
		group_index, group_label, domain_group = {}, {}, {}
		for q in rows:
			for g in by_query[q]:
				for d in domains.get(g.accession, []):
					group = clans.get(d.name) or getattr(d, "group", "") or d.name
					if group not in group_index:
						group_index[group] = len(group_index)
						group_label[group] = group
					domain_group[d.name] = group
		return {
			"number": {name: group_index[grp] + 1 for name, grp in domain_group.items()},
			"color": {name: self._contrast_color(group_index[grp])
					  for name, grp in domain_group.items()},
			"legend": sorted(
				((group_index[grp] + 1, group_label[grp], self._contrast_color(group_index[grp]))
				 for grp in group_index), key=lambda t: t[0]),
		}

	def _secretion_overlay(self, rows, by_query, hits):
		by_loc: Dict[tuple, list] = {}
		for h in hits:
			by_loc.setdefault((h.assembly, h.contig), []).append(h)

		gene_hits: Dict[tuple, list] = {}   
		row_hits: Dict[str, list] = {}      
		seen_types = []
		for q in rows:
			assembly = q.rsplit("|", 1)[-1]
			for g in by_query[q]:
				matches = [h for h in by_loc.get((assembly, g.contig), [])
						   if h.start <= g.end and h.end >= g.start]
				if matches:
					gene_hits[(q, g.offset)] = matches
					for h in matches:
						if h not in row_hits.setdefault(q, []):
							row_hits[q].append(h)
						if h.type not in seen_types:
							seen_types.append(h.type)
		seen_types.sort()
		color = {t: self._contrast_color(i) for i, t in enumerate(seen_types)}
		return {
			"gene_hits": gene_hits,
			"row_hits": row_hits,
			"color": color,
			"legend": [(i + 1, t, color[t]) for i, t in enumerate(seen_types)],
		}

	def _x_for(self, coord, center, q_mid, reversed_row):
		if reversed_row:
			return center + (q_mid - coord) / self.bp_per_px
		return center + (coord - q_mid) / self.bp_per_px

	def _secretion_band(self, hit, color, center, q_mid, reversed_row, y,
						row_lo, row_hi):
		xa = self._x_for(hit.start, center, q_mid, reversed_row)
		xb = self._x_for(hit.end, center, q_mid, reversed_row)
		x0, x1 = (xa, xb) if xa <= xb else (xb, xa)
		x0 = max(x0, row_lo)
		x1 = min(x1, row_hi)
		if x1 <= x0:
			return "", None
		h = self.gene_h + self.SECRETION_BAND_PAD * 2
		rect = ('<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" '
				'fill="{}" fill-opacity="{}" stroke="{}" stroke-opacity="{}" '
				'stroke-width="1" rx="2"/>').format(
					x0, y - h / 2, x1 - x0, h, color, self.SECRETION_BAND_OPACITY,
					color, min(1.0, self.SECRETION_BAND_OPACITY * 2))
		return rect, (x0 + x1) / 2

	def _arrow_points(self, x0, x1, cy, strand):
		h = self.gene_h
		if self.classic:
			length = max(x1 - x0, self.CLASSIC["min_w"])
			head = min(self.CLASSIC["head"], length * 0.35)
		else:
			length = max(x1 - x0, 6)
			head = min(h, length * 0.5)
		top, bot = cy - h / 2, cy + h / 2
		if strand == "-":
			body_l = x0 + head
			return [(x0, cy), (body_l, top), (x0 + length, top),
					(x0 + length, bot), (body_l, bot)]
		body_r = x0 + length - head
		return [(x0, top), (body_r, top), (x0 + length, cy),
				(body_r, bot), (x0, bot)]

	def _gene_outline(self, x0, x1, cy, strand, stroke, sw):
		points = " ".join("{:.1f},{:.1f}".format(px, py)
						  for px, py in self._arrow_points(x0, x1, cy, strand))
		return ('<polygon points="{}" fill="none" stroke="{}" '
				'stroke-width="{}"/>'.format(points, stroke, sw))

	def _gene_fill(self, x0, x1, cy, strand, fill, stroke, sw, opacity=None):
		points = " ".join("{:.1f},{:.1f}".format(px, py)
						  for px, py in self._arrow_points(x0, x1, cy, strand))
		fade = '' if opacity is None else ' fill-opacity="{}"'.format(opacity)
		return ('<polygon points="{}" fill="{}"{} stroke="{}" '
				'stroke-width="{}"/>').format(points, fill, fade, stroke, sw)

	def _domains_on_gene(self, gene, gx0, gx1, cy, drawn_strand, hits, overlay):
		if not hits:
			return "", []
		minus = (drawn_strand == "-")
		u = self.domain_h / 4.0
		prot_len = max((gene.end - gene.start) // 3, 1)
		span = gx1 - gx0

		def res_to_x(res):
			frac = min(max(res / prot_len, 0.0), 1.0)
			return gx1 - frac * span if minus else gx0 + frac * span

		wedges, labels = [], []
		for d in sorted(hits, key=lambda h: (h.start, h.end)):
			s, e = res_to_x(d.start), res_to_x(d.end)
			color = overlay["color"][d.name]
			pts = [(s, cy + 2 * u), (s, cy + 1 * u), (e, cy - 2 * u),
				   (e, cy + 2 * u), (s, cy + 2 * u)]
			points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)
			wedges.append('<polygon points="{}" fill="{}"/>'.format(points, color))
			labels.append(((s + e) / 2, overlay["number"][d.name], color))
		return "".join(wedges), labels

	def _features_on_gene(self, gene, gx0, gx1, cy, drawn_strand, regions):
		if not regions:
			return ""
		minus = (drawn_strand == "-")
		prot_len = max((gene.end - gene.start) // 3, 1)
		span = gx1 - gx0

		def res_to_x(res):
			frac = min(max(res / prot_len, 0.0), 1.0)
			return gx1 - frac * span if minus else gx0 + frac * span

		h = self.gene_h
		out = []
		for kind, start, end in regions:
			x_s, x_e = res_to_x(start), res_to_x(end)
			lo, hi = (x_e, x_s) if x_e < x_s else (x_s, x_e)
			if kind == "tm":
				top, bot = cy - h / 2, cy + h / 2
				n = 4                                   # number of hatch lines
				for k in range(1, n + 1):
					hy = top + (bot - top) * k / (n + 1)
					out.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" '
							   'stroke="#d40000" stroke-width="0.8"/>'.format(lo, hy, hi, hy))
				mid = (lo + hi) / 2
				out.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" '
						   'stroke="#fff" stroke-width="1.2"/>'.format(mid, top, mid, bot))
			elif kind == "signal":
				tip = lo
				w = max(hi - lo, 3) * 0.6
				pts = [(tip, cy - h / 2), (tip - w / 2, cy + h / 2),
					   (tip + w / 2, cy + h / 2)]
				points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)
				out.append('<polygon points="{}" fill="#000"/>'.format(points))
		return "".join(out)

	def _number_in_gene(self, x0, x1, cy, num, fill):
		"""Family number centred inside the arrow, as the original FlaGs3 drew it."""
		width = max(x1 - x0, self.CLASSIC["min_w"])
		size = min(self.font - 1, self.gene_h - 4)
		if self._text_width(str(num), size) > width - 2:
			return ""
		return ('<text x="{:.1f}" y="{:.1f}" font-size="{}" fill="{}" '
				'text-anchor="middle">{}</text>').format(
					(x0 + x1) / 2, cy + size / 3, size,
					self._readable_on(fill), self._escape(str(num)))

	@staticmethod
	def _readable_on(fill: str) -> str:
		"""Black or white, whichever stays legible on the arrow's fill."""
		try:
			r, g, b = (int(fill[i:i + 2], 16) for i in (1, 3, 5))
		except (ValueError, IndexError):
			return "#000"
		return "#000" if (0.299 * r + 0.587 * g + 0.114 * b) > 140 else "#fff"

	def _place_labels(self, row_labels, y):
		step = self.LABEL_STEP
		ly = y - self.gene_h / 2 - 4
		out, placed = [], []
		for lx, num, col in sorted(row_labels, key=lambda t: t[0]):
			if placed and lx - placed[-1] < step:
				lx = placed[-1] + step
			placed.append(lx)
			out.append('<text x="{:.1f}" y="{:.1f}" font-size="{}" fill="{}" '
					   'text-anchor="middle">{}</text>'.format(
						   lx, ly, self.font - 4, col, num))
		return "".join(out)

	MIDGREY = "#bebebe"

	def _extra_legend(self, drawn):
		blocks = []
		feats = []
		if "tm" in drawn:
			feats.append(("Transmembrane helix", "#d40000", "hatch"))
		if "signal" in drawn:
			feats.append(("Signal peptide", "#000000", "triangle"))
		if feats:
			blocks.append(("Protein features", feats))
		marks = []
		if "query" in drawn:
			marks.append(("Query protein", "#000000" if self.classic else "#ffffff",
						  "outline"))
		if "rna" in drawn:
			marks.append(("RNA gene", self.RNA[0], "rna"))
		if "pseudo" in drawn:
			marks.append(("Pseudogene", self.PSEUDO[0], "pseudo"))
		if marks:
			blocks.append(("Genes", marks))
		return blocks

	def _extra_legend_svg(self, blocks, x, y, W):
		parts, ey = [], y
		for title, items in blocks:
			parts.append('<text x="{}" y="{}" font-weight="bold">{}</text>'.format(
				x, ey, title))
			ey += 16
			ex = x
			for label, color, shape in items:
				parts.append(self._legend_mark(ex + 2, ey - 8, color, shape))
				parts.append('<text x="{}" y="{}">{}</text>'.format(
					ex + 20, ey, self._escape(label)))
				ex += int(self._text_width(label, self.font)) + 46
				if ex > W - self.pad - 120:
					ex = x
					ey += 16
			ey += 20
		return "".join(parts), ey - y

	def _legend_mark(self, x, y, color, shape):
		if shape == "hatch":
			return "".join(
				'<line x1="{}" y1="{:.1f}" x2="{}" y2="{:.1f}" stroke="{}" '
				'stroke-width="0.8"/>'.format(x, y + 2 + k * 3, x + 12,
											  y + 2 + k * 3, color)
				for k in range(3))
		if shape == "triangle":
			return ('<polygon points="{},{} {},{} {},{}" fill="{}"/>'.format(
				x + 6, y, x, y + 10, x + 12, y + 10, color))
		stroke = {"outline": "#000", "rna": self.RNA[1],
				  "pseudo": self.PSEUDO[1]}.get(shape, "#333")
		width = 1 if shape == "swatch" else 2
		return ('<rect x="{}" y="{}" width="12" height="10" fill="{}" stroke="{}" '
				'stroke-width="{}"/>'.format(x, y, color, stroke, width))

	def _legend_layout(self, overlay, W):
		items = overlay.get("legend")
		if not items:
			return None, 0, 0, 0
		longest = max(self._text_width("{}. {}".format(num, lbl), self.font)
					  for num, lbl, _ in items)
		col_w = int(longest) + 26
		avail = W - 2 * self.pad
		cols = max(1, min(len(items), avail // col_w))
		legend_rows = -(-len(items) // cols)
		legend_h = (legend_rows + 1) * 16 + 8
		return items, cols, col_w, legend_h

	def _legend(self, items, x, y, cols, col_w, title="Domains"):
		parts = ['<text x="{}" y="{}" font-weight="bold">{}</text>'.format(x, y, title)]
		for idx, (num, label, color) in enumerate(items):
			ex = x + (idx % cols) * col_w
			ey = y + 16 + (idx // cols) * 16
			parts.append('<rect x="{}" y="{}" width="12" height="10" fill="{}"/>'.format(
				ex + 2, ey - 8, color))
			parts.append('<text x="{}" y="{}">{}. {}</text>'.format(
				ex + 20, ey, num, self._escape(label)))
		return "".join(parts)

	@staticmethod
	def _contrast_color(n):
		hue = (n * 0.61803398875) % 1.0
		r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.75)
		return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

	@staticmethod
	def _row_reversed(genes):
		q = next((g for g in genes if g.offset == 0), None)
		if q is None:
			return False
		pos = [g for g in genes if g.offset > 0]
		if not pos:
			neg = [g for g in genes if g.offset < 0]
			if not neg:
				return False
			ref = max(neg, key=lambda g: g.offset)
			return ref.start > q.start
		ref = max(pos, key=lambda g: g.offset)
		return ref.start < q.start


class NeighborhoodVisualizer(_FlaGsBase): 

  def __init__(self, gene_h: int = 20, gene_gap: int = 1,
      row_h: int = 24, tree_w: int = 320, pad: int = 16, font: int = 13):
    self.gene_h = gene_h
    self.gene_w = gene_h * 0.95
    self.gene_gap = gene_gap
    self.row_h = row_h; self.tree_w = tree_w; self.pad = pad
    self.font = font

  def render(self, newick: str,
      neighborhoods: list,
      families: List[List[str]],
      species: Optional[Dict[str, str]] = None,
      labels: Optional[Dict[str, str]] = None) -> str:
    species = species or {}
    row_labels_map = labels or {}
    by_query: Dict[str, list] = {}
    for g in neighborhoods:
      by_query.setdefault(g.query, []).append(g)
    for q in by_query:
      by_query[q].sort(key=lambda x: x.offset)

    color = self._family_colors(families)
    rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
    number = self._numbers_for(families, rna_accessions)
    tree = Phylo.read(StringIO(newick), "newick")
    try:
      tree.root_at_midpoint()
    except Exception:
      pass    
    tree.ladderize()
    tip_order = [t.name for t in tree.get_terminals()]
    rows = [q for q in tip_order if q in by_query] or list(by_query)
    def _lbl(q):
      if q in row_labels_map:
        return row_labels_map[q]
      return q if not species.get(q) else "{}  {}".format(q, species[q])
    labels = {q: _lbl(q) for q in rows}
    label_w = int(max((self._text_width(l, self.font) for l in labels.values()), default=0)) + 14

    cell = self.gene_w + self.gene_gap
    max_off = max((abs(g.offset) for gs in by_query.values() for g in gs), default=0)
    track_w = (2 * max_off + 1) * cell

    W = self.pad + self.tree_w + label_w + track_w + self.pad
    H = self.pad * 2 + len(rows) * self.row_h + 30

    y_of = {q: self.pad + i * self.row_h + self.row_h / 2 for i, q in enumerate(rows)}

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
      'font-family="{}" font-size="{}">'.format(W, H, self.FONT_FAMILY, self.font),
      '<rect width="{}" height="{}" fill="white"/>'.format(W, H)]

    tree_svg, xscale, maxd = self._draw_tree(tree, y_of)
    svg.append(tree_svg)

    track_x0 = self.pad + self.tree_w + label_w
    center = track_x0 + track_w / 2  

    for q in rows:
      y = y_of[q]
      svg.append('<text x="{}" y="{}">{}</text>'.format(
        self.pad + self.tree_w + 6, y + self.font / 3, self._escape(labels[q])))
      for g in by_query[q]:
        cx = center + g.offset * cell
        fill, outline = self._style_for(g.accession, color, g.is_rna)
        stroke = outline
        sw = self._stroke_width(g)
        svg.append(self._arrow(cx, y, g.strand, fill, stroke, sw,
          number.get(g.accession)))

    base_y = self.pad + len(rows) * self.row_h + 16
    if xscale and maxd:
      svg.append(self._scale_bar(self.pad, base_y, xscale, maxd))
    svg.append('</svg>')
    return "\n".join(svg)

  def _arrow(self, cx, cy, strand, fill, stroke, sw, number=None) -> str:
    hw, hh = self.gene_w / 2, self.gene_h / 2
    if strand == "+":
      pts = [(cx - hw, cy - hh), (cx + hw, cy), (cx - hw, cy + hh)]
    else:
      pts = [(cx + hw, cy - hh), (cx - hw, cy), (cx + hw, cy + hh)]
    points = " ".join("{},{}".format(x, y) for x, y in pts)
    out = ['<polygon points="{}" fill="{}" stroke="{}" stroke-width="{}"/>'.format(
      points, fill, stroke, sw)]
    if number is not None:
      tx = cx - hw * 0.2 if strand == "+" else cx + hw * 0.2
      out.append('<text x="{}" y="{}" font-size="{}" fill="black" '
        'text-anchor="middle" dominant-baseline="central">{}</text>'.format(
          tx, cy, self.font - 5, number))
    return "".join(out)