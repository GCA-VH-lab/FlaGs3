FONT_FAMILY = "Arial, 'Liberation Sans', Helvetica, sans-serif"

CHAR_WIDTH = {
	' ': 278, '!': 278, '"': 355, '#': 556, '$': 556, '%': 889, '&': 667, "'": 191,
	'(': 333, ')': 333, '*': 389, '+': 584, ',': 278, '-': 333, '.': 278, '/': 278,
	'0': 556, '1': 556, '2': 556, '3': 556, '4': 556, '5': 556, '6': 556, '7': 556, '8': 556, '9': 556,
	':': 278, ';': 278, '<': 584, '=': 584, '>': 584, '?': 556, '@': 1015,
	'A': 667, 'B': 667, 'C': 722, 'D': 722, 'E': 667, 'F': 611, 'G': 778, 'H': 722, 'I': 278,
	'J': 500, 'K': 667, 'L': 556, 'M': 833, 'N': 722, 'O': 778, 'P': 667, 'Q': 778, 'R': 722,
	'S': 667, 'T': 611, 'U': 722, 'V': 667, 'W': 944, 'X': 667, 'Y': 667, 'Z': 611,
	'[': 278, '\\': 278, ']': 278, '^': 469, '_': 556, '`': 333,
	'a': 556, 'b': 556, 'c': 500, 'd': 556, 'e': 556, 'f': 278, 'g': 556, 'h': 556, 'i': 222,
	'j': 222, 'k': 500, 'l': 222, 'm': 833, 'n': 556, 'o': 556, 'p': 556, 'q': 556, 'r': 333,
	's': 500, 't': 278, 'u': 556, 'v': 500, 'w': 722, 'x': 500, 'y': 500, 'z': 500,
	'{': 334, '|': 260, '}': 334, '~': 584,
}
DEFAULT_WIDTH = 556


def width(text: str, size: float) -> float:
	return sum(CHAR_WIDTH.get(ch, DEFAULT_WIDTH) for ch in text) / 1000.0 * size


def escape(text: str) -> str:
	return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x: float, y: float, content: str, size=None, anchor=None, fill=None, weight=None) -> str:
	attrs = ['x="{:.1f}"'.format(x), 'y="{:.1f}"'.format(y)]
	if size is not None:
		attrs.append('font-size="{}"'.format(size))
	if anchor:
		attrs.append('text-anchor="{}"'.format(anchor))
	if fill:
		attrs.append('fill="{}"'.format(fill))
	if weight:
		attrs.append('font-weight="{}"'.format(weight))
	return "<text {}>{}</text>".format(" ".join(attrs), escape(content))
