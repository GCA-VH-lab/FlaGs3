import csv
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Optional

from flags3.render import palettes
from flags3.schema import MISSING, Annotation

TABLE_NAME = "visualisation_table.tsv"
STAGE_TITLES = {
	"cluster": "Gene families", "cluster_rna": "RNA families", "domains": "Domains",
	"features": "Protein features", "sismis": "Secretion systems", "genomad": "Mobile elements",
	"defence": "Defence systems",
}
BAND_CODES = {"sismis": "S", "genomad": "M", "defence": "D"}
BAND_ORDER = ("sismis", "genomad", "defence")
GEOMETRY = {
	"font_size": 13, "row_height": 26, "gene_height": 8, "pad": 16, "domain_height": 6,
	"bases_per_pixel": 10.4, "label_step": 12, "arrow_head": 7, "min_gene_width": 13,
	"band_opacity": 0.85, "tree_width": 320,
}
CLASSIC_GEOMETRY = {"font_size": 12, "row_height": 20, "gene_height": 15, "bases_per_pixel": 16.0}


class StyleError(ValueError):
	pass


@dataclass
class FigureSpec:
	name: str
	layers: list[str]
	mode: str = "versatile"
	tree: bool = False
	numbers: bool = True
	palette: str = "bright"
	monochrome: bool = False
	geometry: dict = field(default_factory=dict)

	@property
	def classic(self) -> bool:
		return self.mode == "classic"

	def value(self, key: str):
		if key in self.geometry:
			return self.geometry[key]
		if self.classic and key in CLASSIC_GEOMETRY:
			return CLASSIC_GEOMETRY[key]
		return GEOMETRY[key]

	def wants(self, stage: str) -> bool:
		return stage in self.layers


def _bool(raw: str) -> bool:
	return raw.strip().lower() in ("1", "true", "yes")


def parse_row(row: dict[str, str], origin: str) -> FigureSpec:
	name = row.get("name", "").strip()
	if not name:
		raise StyleError("{}: a figure row has no name".format(origin))
	layers = [l.strip() for l in row.get("layers", "").split(",") if l.strip()]
	mode = (row.get("mode") or "versatile").strip().lower()
	if mode not in ("versatile", "classic"):
		raise StyleError("{}: figure {} has mode {!r}; versatile or classic".format(origin, name, mode))
	palette = (row.get("palette") or "bright").strip().lower()
	if palette not in palettes.PALETTES:
		raise StyleError("{}: figure {} has palette {!r}; one of {}".format(origin, name, palette, ", ".join(palettes.PALETTES)))
	geometry = {}
	for key in GEOMETRY:
		raw = (row.get(key) or "").strip()
		if raw and raw.lower() != "default":
			try:
				geometry[key] = float(raw) if "." in raw else int(raw)
			except ValueError:
				raise StyleError("{}: figure {} has {}={!r}, not a number".format(origin, name, key, raw))
	return FigureSpec(name, layers, mode, _bool(row.get("tree", "")), _bool(row.get("numbers", "true")), palette,
		_bool(row.get("monochrome", "")), geometry)


def read_table(path: Optional[Path]) -> list[FigureSpec]:
	if path is None:
		with resources.as_file(resources.files("flags3.data") / TABLE_NAME) as shipped:
			return read_table(shipped)
	specs = []
	with open(path, newline="", encoding="utf-8") as handle:
		reader = csv.DictReader((ln for ln in handle if ln.strip()), delimiter="\t")
		for raw in reader:
			row = {(k or "").lstrip("#").strip(): (v or "").strip() for k, v in raw.items()}
			specs.append(parse_row(row, str(path)))
	return specs


class Colours:
	def __init__(self, palette: str, overrides: Optional[dict[tuple[str, str], str]] = None, monochrome: bool = False):
		self.palette = palettes.PALETTES[palette]
		self.overrides = overrides or {}
		self.monochrome = monochrome
		self.assigned: dict[str, dict[str, str]] = {}

	def assign(self, stage: str, annotations: list[Annotation]) -> dict[str, str]:
		categories = []
		for a in annotations:
			if a.category not in categories:
				categories.append(a.category)
		chooser = palettes.monochrome if (self.monochrome and stage.startswith("cluster")) else self.palette
		colours = chooser(len(categories))
		table = {c: self.overrides.get((stage, c), colours[i]) for i, c in enumerate(categories)}
		self.assigned[stage] = table
		return table

	def of(self, stage: str, category: str) -> str:
		return self.assigned.get(stage, {}).get(category, palettes.GREY)


def read_overrides(path: Path) -> dict[tuple[str, str], str]:
	out = {}
	if not path.is_file():
		return out
	with open(path, newline="", encoding="utf-8") as handle:
		for raw in csv.DictReader((ln for ln in handle if ln.strip()), delimiter="\t"):
			row = {(k or "").lstrip("#").strip(): (v or "").strip() for k, v in raw.items()}
			if row.get("stage") and row.get("category") and row.get("colour"):
				out[(row["stage"], row["category"])] = row["colour"]
	return out
