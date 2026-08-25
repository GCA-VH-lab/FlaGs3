import importlib.util
import os
import shutil
import subprocess

from flags_log import debug

BACKENDS = ("cairosvg", "svglib", "rsvg-convert", "inkscape")


def _cairosvg(svg_path: str, pdf_path: str):
	import cairosvg
	cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)


def _svglib(svg_path: str, pdf_path: str):
	from svglib.svglib import svg2rlg
	from reportlab.graphics import renderPDF
	drawing = svg2rlg(svg_path)
	if drawing is None:
		raise RuntimeError("svglib could not parse the SVG")
	renderPDF.drawToFile(drawing, pdf_path)


def _binary(name):
	def run(svg_path: str, pdf_path: str):
		exe = shutil.which(name)
		if exe is None:
			raise FileNotFoundError("{} is not on PATH".format(name))
		cmd = ([exe, "-f", "pdf", "-o", pdf_path, svg_path] if name == "rsvg-convert"
			   else [exe, "--export-type=pdf", "--export-filename=" + pdf_path, svg_path])
		result = subprocess.run(cmd, capture_output=True, text=True)
		if result.returncode != 0:
			raise RuntimeError("{} exited {}: {}".format(
				name, result.returncode, (result.stderr or "").strip()[:200]))
	return run


_RUNNERS = {"cairosvg": _cairosvg, "svglib": _svglib,
			"rsvg-convert": _binary("rsvg-convert"), "inkscape": _binary("inkscape")}


def available() -> str:
	for name in BACKENDS:
		if name in ("cairosvg", "svglib"):
			needed = ("cairosvg",) if name == "cairosvg" else ("svglib", "reportlab")
			if all(importlib.util.find_spec(m) for m in needed):
				return name
		elif shutil.which(name):
			return name
	return ""


def convert(svg_path: str, pdf_path: str = None, backend: str = None) -> str:
	pdf_path = pdf_path or os.path.splitext(svg_path)[0] + ".pdf"
	chosen = backend or available()
	if not chosen:
		raise RuntimeError(
			"no SVG to PDF backend found. Install one of: cairosvg, svglib, "
			"librsvg (rsvg-convert), inkscape.")
	if chosen not in _RUNNERS:
		raise ValueError("unknown --pdf backend {!r}; choose from {}".format(
			chosen, ", ".join(BACKENDS)))
	_RUNNERS[chosen](svg_path, pdf_path)
	debug("pdf: {} -> {} via {}".format(os.path.basename(svg_path),
										os.path.basename(pdf_path), chosen))
	return pdf_path


def convert_all(svg_paths, backend: str = None):
	chosen = backend or available()
	if not chosen:
		print("Warning: --pdf needs one of cairosvg, svglib, rsvg-convert or "
			  "inkscape; writing SVG only.")
		return []
	written = []
	for svg_path in svg_paths:
		try:
			written.append(convert(svg_path, backend=chosen))
		except Exception as e:
			print("Warning: could not convert {} to PDF ({}).".format(
				os.path.basename(svg_path), e))
	return written
