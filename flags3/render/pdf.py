import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

from flags3.render.text import CHAR_WIDTH, DEFAULT_WIDTH

NS = "{http://www.w3.org/2000/svg}"
NAMED = {"white": "#ffffff", "black": "#000000"}
INHERITED = ("stroke", "stroke-width", "stroke-opacity")
KAPPA = 0.5523


def colour(value: str) -> tuple | None:
	value = NAMED.get(value, value)
	if value == "none":
		return None
	hexes = value[1:]
	if len(hexes) == 3:
		hexes = "".join(ch * 2 for ch in hexes)
	return tuple(int(hexes[i:i + 2], 16) / 255 for i in (0, 2, 4))


def rgb(c: tuple, op: str) -> str:
	return "{:.3f} {:.3f} {:.3f} {}".format(c[0], c[1], c[2], op)


def escape(content: str) -> str:
	return content.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


class Writer:
	def __init__(self, svg: Path):
		self.name = svg.name
		self.root = ET.parse(svg).getroot()
		self.width = float(self.root.get("width"))
		self.height = float(self.root.get("height"))
		self.font = float(self.root.get("font-size", "13"))
		self.clips = {clip.get("id"): clip for clip in self.root.iter(NS + "clipPath")}
		self.alphas = {}
		self.ops = ["1 0 0 -1 0 {:.2f} cm".format(self.height)]
		self.replaced = set()

	def run(self, pdf: Path) -> None:
		self.group(self.root, {"stroke": "none", "stroke-width": "1", "stroke-opacity": "1"})
		if self.replaced:
			print("Warning: {} has characters outside the PDF text encoding ({}); they are written as ? in the PDF.".format(
				self.name, " ".join(sorted(self.replaced))))
		self.write(pdf)

	def group(self, el: ET.Element, style: dict) -> None:
		for child in el:
			tag = child.tag[len(NS):]
			if tag in ("defs", "clipPath"):
				continue
			local = dict(style)
			for key in INHERITED:
				if child.get(key) is not None:
					local[key] = child.get(key)
			clip = child.get("clip-path")
			if clip:
				self.ops.append("q")
				self.clip(self.clips[clip[5:-1]])
			if tag == "g":
				self.group(child, local)
			else:
				draw = getattr(self, "draw_" + tag, None)
				if draw is None:
					raise ValueError("{}: the PDF writer does not draw <{}>".format(self.name, tag))
				draw(child, local)
			if clip:
				self.ops.append("Q")

	def clip(self, clip: ET.Element) -> None:
		for shape in clip:
			self.ops.append(self.outline(shape) + " W n")

	def outline(self, el: ET.Element) -> str:
		tag = el.tag[len(NS):]
		if tag == "rect":
			return "{} {} {} {} re".format(el.get("x", "0"), el.get("y", "0"), el.get("width"), el.get("height"))
		if tag == "polygon":
			pts = [p.split(",") for p in el.get("points").split()]
			return "{} {} m ".format(*pts[0]) + " ".join("{} {} l".format(*p) for p in pts[1:]) + " h"
		cx, cy, r = (float(el.get(k)) for k in ("cx", "cy", "r"))
		k = KAPPA * r
		arcs = (
			(cx + r, cy + k, cx + k, cy + r, cx, cy + r), (cx - k, cy + r, cx - r, cy + k, cx - r, cy),
			(cx - r, cy - k, cx - k, cy - r, cx, cy - r), (cx + k, cy - r, cx + r, cy - k, cx + r, cy))
		return "{} {} m ".format(cx + r, cy) + " ".join("{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} c".format(*a) for a in arcs) + " h"

	def paint(self, el: ET.Element, style: dict) -> None:
		fill = colour(el.get("fill", "black"))
		stroke = colour(style["stroke"])
		if fill is None and stroke is None:
			return
		ops = ["q"]
		if fill:
			ops.append(rgb(fill, "rg"))
		if stroke:
			ops.append(rgb(stroke, "RG") + " {} w /{} gs".format(style["stroke-width"], self.alpha(style["stroke-opacity"])))
		ops.append("{} {} Q".format(self.outline(el), "B" if fill and stroke else ("f" if fill else "S")))
		self.ops.append(" ".join(ops))

	draw_rect = paint
	draw_polygon = paint
	draw_circle = paint

	def draw_line(self, el: ET.Element, style: dict) -> None:
		stroke = colour(style["stroke"])
		if stroke is None:
			return
		self.ops.append("q {} {} w /{} gs {} {} m {} {} l S Q".format(rgb(stroke, "RG"), style["stroke-width"],
			self.alpha(style["stroke-opacity"]), el.get("x1"), el.get("y1"), el.get("x2"), el.get("y2")))

	def draw_text(self, el: ET.Element, style: dict) -> None:
		content = " ".join((el.text or "").split())
		self.replaced.update(ch for ch in content if not ch.isascii() and ch.encode("cp1252", "replace") == b"?")
		content = content.encode("cp1252", "replace").decode("cp1252")
		size = float(el.get("font-size", self.font))
		x, y = float(el.get("x")), float(el.get("y"))
		if el.get("text-anchor") == "middle":
			x -= sum(CHAR_WIDTH.get(ch, DEFAULT_WIDTH) for ch in content) / 2000.0 * size
		font = "F2" if el.get("font-weight") == "bold" else "F1"
		self.ops.append("BT {} /{} {:g} Tf 1 0 0 -1 {:.2f} {:.2f} Tm ({}) Tj ET".format(
			rgb(colour(el.get("fill", "black")), "rg"), font, size, x, y, escape(content)))

	def alpha(self, value: str) -> str:
		return self.alphas.setdefault("{:.3f}".format(float(value)), "GS{}".format(len(self.alphas)))

	def write(self, pdf: Path) -> None:
		stream = zlib.compress("\n".join(self.ops).encode("latin-1"))
		states = " ".join("/{} {} 0 R".format(name, 7 + i) for i, name in enumerate(self.alphas.values()))
		objects = [
			"<< /Type /Catalog /Pages 2 0 R >>",
			"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
			"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {:.2f} {:.2f}] /Contents 4 0 R "
			"/Resources << /Font << /F1 5 0 R /F2 6 0 R >> /ExtGState << {} >> >> >>".format(self.width, self.height, states),
			b"<< /Length %d /Filter /FlateDecode >>\nstream\n" % len(stream) + stream + b"\nendstream",
			"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>",
			"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>",
		] + ["<< /Type /ExtGState /CA {} >>".format(value) for value in self.alphas]
		out = bytearray(b"%PDF-1.4\n")
		offsets = []
		for i, obj in enumerate(objects, 1):
			offsets.append(len(out))
			out += b"%d 0 obj\n" % i + (obj if isinstance(obj, bytes) else obj.encode("latin-1")) + b"\nendobj\n"
		xref = len(out)
		out += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objects) + 1)
		out += b"".join(b"%010d 00000 n \n" % o for o in offsets)
		out += b"trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objects) + 1, xref)
		pdf.write_bytes(out)


def convert(svg: Path, pdf: Path) -> None:
	Writer(svg).run(pdf)
