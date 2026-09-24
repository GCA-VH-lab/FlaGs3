from pathlib import Path


def available() -> bool:
	try:
		import cairosvg
		return True
	except ImportError:
		return False


def convert(svg: Path, pdf: Path) -> None:
	import cairosvg
	cairosvg.svg2pdf(url=str(svg), write_to=str(pdf))
