import colorsys

GREY = "#d9d9d9"
MIDGREY = "#bebebe"
WHITE = "#ffffff"
PSEUDO = ("#f2f2f3", "#000080")
RNA = ("#f2f2f2", "#008000")
QUERY = "#000000"
OKABE_ITO = ("#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7", "#999999")
TOL_LIGHT = ("#77aadd", "#ee8866", "#eedd88", "#ffaabb", "#99ddff", "#44bb99", "#bbcc33", "#aaaa00", "#dddddd")


def _hex(r, g, b) -> str:
	return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


GOLDEN = 0.618033988749895


def bright(n: int) -> list[str]:
	out = []
	for i in range(n):
		hue = (i * GOLDEN) % 1.0
		value = 0.85 if (i // 7) % 2 == 0 else 0.68
		out.append(_hex(*colorsys.hsv_to_rgb(hue, 0.55, value)))
	return out


def pastel(n: int) -> list[str]:
	return [lighten(c, 0.45) for c in bright(n)]


def classic(n: int) -> list[str]:
	hues = [int(h * 3.6) / 100.0 for h in range(0, 100, 5)]
	return [_hex(*colorsys.hls_to_rgb(hues[i % len(hues)], 0.5, 0.5)) for i in range(n)]


def colourblind(n: int) -> list[str]:
	base = OKABE_ITO + TOL_LIGHT
	return [base[i % len(base)] for i in range(n)]


def monochrome(n: int) -> list[str]:
	return [GREY for _ in range(n)]


PALETTES = {"bright": bright, "pastel": pastel, "classic": classic, "colourblind": colourblind,
	"colorblind": colourblind, "monochrome": monochrome}


def lighten(colour: str, amount: float) -> str:
	try:
		r, g, b = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
	except (ValueError, IndexError):
		return colour
	mix = lambda c: int(round(c + (255 - c) * amount))
	return "#{:02x}{:02x}{:02x}".format(mix(r), mix(g), mix(b))


def readable_on(fill: str) -> str:
	try:
		r, g, b = (int(fill[i:i + 2], 16) for i in (1, 3, 5))
	except (ValueError, IndexError):
		return "#000"
	return "#000" if (0.299 * r + 0.587 * g + 0.114 * b) > 140 else "#fff"


def is_pale(fill: str) -> bool:
	return fill in (GREY, WHITE, PSEUDO[0], RNA[0], MIDGREY)
