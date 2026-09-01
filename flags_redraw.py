#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
from typing import Dict, List, NamedTuple, Optional, Set

TABLE_NAME = "visualisation_table.tsv"

COLUMNS = ("name", "mode", "tree_width", "features_allowed", "family_numbers",
		   "font_size", "row_height", "gene_height", "gene_gap",
		   "bases_per_pixel", "pad", "domain_height", "label_step",
		   "arrow_head", "min_gene_width", "band_opacity")

NUMERIC = {"font_size": int, "row_height": int, "gene_height": int,
		   "gene_gap": int, "bases_per_pixel": float, "pad": int,
		   "domain_height": int, "label_step": int, "arrow_head": float,
		   "min_gene_width": float, "band_opacity": float}

MODES = ("versatile", "triangles", "classic")

FEATURES = ("cluster_rna", "domains", "tmhmm", "signalp", "sismis",
			"monochrome", "none")

FALSEY = ("false", "no", "off", "0", "")
DEFAULTY = ("", "-", "default", "from the code", "auto")

def default_table_path() -> str:
	return os.path.join(os.path.dirname(os.path.abspath(__file__)), TABLE_NAME)


def read_default_table() -> str:
	path = default_table_path()
	if not os.path.isfile(path):
		raise FileNotFoundError(
			"no {} next to {}. It defines the standard figures; restore it or "
			"pass --format.".format(TABLE_NAME, os.path.basename(__file__)))
	with open(path, encoding="utf-8") as fh:
		return fh.read()


class FigureSpec(NamedTuple):
	name: str
	mode: str = "versatile"
	tree_width: Optional[float] = None      # None = no tree panel
	features: Set[str] = frozenset()
	family_numbers: bool = True
	geometry: Dict[str, float] = {}         # every numeric knob, absent = default

	@property
	def monochrome(self) -> bool:
		return "monochrome" in self.features

	def wants(self, feature: str) -> bool:
		return feature in self.features


def _is_default(value: str) -> bool:
	return value.strip().lower() in DEFAULTY


def _as_bool(value: str) -> bool:
	return value.strip().lower() not in FALSEY


def _as_number(value: str, cast, field: str, name: str):
	if _is_default(value):
		return None
	if value.strip().lower() in FALSEY:
		return False
	try:
		return cast(value)
	except ValueError:
		raise ValueError("figure {!r}: {} must be a number, False, or 'default' "
						 "(got {!r})".format(name, field, value))


def parse_row(row: Dict[str, str]) -> FigureSpec:
	name = (row.get("name") or "").strip()
	if not name:
		raise ValueError("a row has no name")
	if not re.match(r"^[A-Za-z0-9._-]+$", name):
		raise ValueError("figure {!r}: name is used as a filename, so keep it to "
						 "letters, digits, dot, dash and underscore".format(name))

	mode = (row.get("mode") or "versatile").strip().lower() or "versatile"
	if mode not in MODES:
		raise ValueError("figure {!r}: mode must be one of {} (got {!r})".format(
			name, ", ".join(MODES), mode))

	raw = (row.get("features_allowed") or "").strip()
	features = set()
	for token in re.split(r"[,;\s]+", raw):
		token = token.strip().lower()
		if not token or token == "none":
			continue
		if token not in FEATURES:
			raise ValueError("figure {!r}: unknown feature {!r}; allowed: {}".format(
				name, token, ", ".join(FEATURES)))
		features.add(token)

	tree_width = _as_number(row.get("tree_width", ""), float, "tree_width", name)
	if tree_width is False:
		tree_width = None          # False = no tree panel
	elif tree_width is None:
		tree_width = 1.0 if mode in ("triangles", "classic") else None
	if mode == "triangles" and tree_width is None:
		tree_width = 1.0           # a triangles figure without a tree is just rows

	bpp = _as_number(row.get("bases_per_pixel", ""), float, "bases_per_pixel", name)
	if bpp is not None and bpp is not False and bpp <= 0:
		raise ValueError("figure {!r}: bases_per_pixel must be positive".format(name))

	geometry = {}
	for field, cast in NUMERIC.items():
		if field == "bases_per_pixel":
			continue
		value = _as_number(row.get(field, ""), cast, field, name)
		if value is not None and value is not False:
			geometry[field] = value
	if bpp is not None:
		geometry["bases_per_pixel"] = bpp
	return FigureSpec(
		name=name, mode=mode, tree_width=tree_width, features=frozenset(features),
		family_numbers=_as_bool(row.get("family_numbers", "TRUE")),
		geometry=geometry)


def read_table(path: str) -> List[FigureSpec]:
	with open(path, newline="", encoding="utf-8") as fh:
		rows = list(csv.DictReader(
			(ln for ln in fh if ln.strip() and not ln.startswith("##")),
			delimiter="\t"))
	if not rows:
		raise ValueError("{} has no figure rows".format(path))
	fields = [(f or "").lstrip("#").strip() for f in (rows[0].keys())]
	if "name" not in fields:
		raise ValueError("{}: first column must be 'name' (or '#name')".format(path))
	specs, seen = [], set()
	for row in rows:
		clean = {(k or "").lstrip("#").strip(): (v or "") for k, v in row.items()}
		spec = parse_row(clean)
		if spec.name in seen:
			raise ValueError("figure {!r} is listed twice".format(spec.name))
		seen.add(spec.name)
		specs.append(spec)
	return specs


def default_specs() -> List[FigureSpec]:
	path = default_table_path()
	if not os.path.isfile(path):
		raise FileNotFoundError(
			"no {} next to flags_redraw.py. It defines the standard figures; "
			"restore it or pass --format.".format(TABLE_NAME))
	return read_table(path)


def write_default_table(path: str):
	with open(path, "w") as out:
		out.write(read_default_table())



class Gene(NamedTuple):
	"""Stands in for FlaGs3.FlankingGene, rebuilt from _operon.tsv."""
	accession: str
	strand: str
	start: int
	end: int
	product: str
	offset: int
	query: str
	is_rna: bool = False
	contig: str = ""


class Hit(NamedTuple):
	assembly: str
	contig: str
	start: int
	end: int
	type: str
	probability: float = 0.0


class RunData(NamedTuple):
	genes: List[Gene]
	families: List[List[str]]
	species: Dict[str, str]
	labels: Dict[str, str]
	order: Optional[List[str]]
	newick: str
	numbers: Dict[str, str]
	domains: Dict[str, list]
	clans: Dict[str, str]
	features: Dict[str, list]
	secretion: List[Hit]


def _rows(path: str):
	with open(path, newline="", encoding="utf-8") as fh:
		for row in csv.DictReader(
				(ln.lstrip("#") if i == 0 else ln
				 for i, ln in enumerate(fh) if ln.strip()), delimiter="\t"):
			yield {(k or "").strip(): (v or "").strip() for k, v in row.items()}


def load_run(directory: str, prefix: str = None) -> RunData:
	directory = os.path.abspath(directory)
	if prefix is None:
		prefix = os.path.basename(os.path.normpath(directory))
	path = lambda suffix: os.path.join(directory, prefix + suffix)

	operon = path("_operon.tsv")
	if not os.path.isfile(operon):
		raise FileNotFoundError(
			"{} does not look like a FlaGs3 output directory: no {}".format(
				directory, os.path.basename(operon)))

	genes, species, labels, seen_rows, numbers = [], {}, {}, [], {}
	label_assembly = False
	for row in _rows(operon):
		query, assembly = row.get("query", ""), row.get("assembly", "")
		row_id = "{}|{}".format(query, assembly) if assembly and assembly != "-" else query
		if row_id not in labels:
			seen_rows.append(row_id)
		genes.append(Gene(
			accession=row.get("accession", ""), strand=row.get("strand", "+") or "+",
			start=int(row.get("start") or 0), end=int(row.get("end") or 0),
			product=row.get("product", ""), offset=int(row.get("offset") or 0),
			query=row_id,
			is_rna=(row.get("is_rna", "") or "").lower() in ("true", "1", "yes"),
			contig="" if row.get("contig") in (None, "-", "") else row["contig"]))
		if row.get("species") and row["species"] != "-":
			species[row_id] = row["species"]
		labels[row_id] = query
		fam = row.get("family", "-")
		if fam and fam != "-":
			numbers[row.get("accession", "")] = fam
	base = {}
	for row_id in seen_rows:
		base.setdefault(labels[row_id], []).append(row_id)
	label_assembly = any(len(v) > 1 for v in base.values())
	for row_id in seen_rows:
		name = row_id if label_assembly else labels[row_id]
		labels[row_id] = ("{}  {}".format(name, species[row_id])
						  if species.get(row_id) else name)

	families = []
	clusters = path("_clusters.tsv")
	if os.path.isfile(clusters):
		for row in _rows(clusters):
			members = [m for m in re.split(r"[,;\s]+", row.get("members", "")) if m]
			if members:
				families.append(members)

	domains, clans = {}, {}
	dom_path = path("_domains.tsv")
	if os.path.isfile(dom_path):
		from flags_domains import DomainHit
		for row in _rows(dom_path):
			protein = row.get("protein", "")
			if not protein:
				continue
			domains.setdefault(protein, []).append(DomainHit(
				protein=protein, name=row.get("domain", ""),
				start=int(row.get("start") or 0), end=int(row.get("end") or 0),
				evalue=float(row.get("evalue") or 0.0),
				accession=row.get("pfam", ""),
				database=row.get("database", ""),
				group=row.get("group", "")))
			if row.get("clan") and row["clan"] != "-":
				clans[row.get("domain", "")] = row["clan"]

	features = {}
	feat_path = path("_features.tsv")
	if os.path.isfile(feat_path):
		for row in _rows(feat_path):
			protein = row.get("protein", "")
			if protein:
				features.setdefault(protein, []).append(
					(row.get("kind", ""), int(row.get("start") or 0),
					 int(row.get("end") or 0)))

	secretion = []
	sec_path = path("_secretion.tsv")
	if os.path.isfile(sec_path):
		for row in _rows(sec_path):
			if not row.get("contig") or not (row.get("start") or "").isdigit():
				continue      # older files without the normalised columns
			secretion.append(Hit(
				assembly=row.get("assembly", ""), contig=row["contig"],
				start=int(row["start"]), end=int(row.get("end") or 0),
				type=row.get("type", "") or "?",
				probability=float(row.get("probability") or 0.0)))

	newick, order = "", None
	tree_path = os.path.join(directory, "tree", prefix + "_tree.nwk")
	if not os.path.isfile(tree_path):
		tree_path = path("_tree.nwk")
	if os.path.isfile(tree_path):
		with open(tree_path) as fh:
			newick = fh.read().strip()
		if newick:
			try:
				from flags_tree import ladderized_leaf_order
				order = [r for r in ladderized_leaf_order(newick) if r in labels]
			except Exception:
				order = None

	return RunData(genes=genes, families=families, species=species, labels=labels,
				   numbers=numbers,
				   order=order or seen_rows, newick=newick, domains=domains,
				   clans=clans, features=features, secretion=secretion)


def render_figure(spec: FigureSpec, data: RunData) -> str:
	"""Draw one figure from a spec. Returns SVG text, or '' if nothing applies."""
	from flags_view import OperonView, NeighborhoodVisualizer

	g = dict(spec.geometry)
	classic = spec.mode == "classic"
	pick = lambda key, default: g.get(key, default)

	domains = data.domains if spec.wants("domains") else None
	secretion = data.secretion if spec.wants("sismis") else None
	features = {}
	if spec.wants("tmhmm") or spec.wants("signalp"):
		want = {k for k in ("tmhmm", "signalp") if spec.wants(k)}
		for protein, regions in data.features.items():
			kept = [r for r in regions if _feature_kind(r[0]) in want]
			if kept:
				features[protein] = kept
	families = data.families
	if not spec.wants("cluster_rna"):
		rna = {g_.accession for g_ in data.genes if g_.is_rna}
		families = [fam for fam in families if not set(fam) <= rna]

	if spec.mode == "triangles":
		if not data.newick:
			return ""
		viz = NeighborhoodVisualizer(
			gene_h=int(pick("gene_height", 20)), gene_gap=int(pick("gene_gap", 1)),
			row_h=int(pick("row_height", 24)), pad=int(pick("pad", 16)),
			font=int(pick("font_size", 13)),
			tree_w=int(320 * (spec.tree_width or 1.0)))
		viz.numbers = data.numbers
		viz.show_numbers = spec.family_numbers
		viz.monochrome = spec.monochrome
		return viz.render(data.newick, data.genes, families,
						  data.species, labels=data.labels)

	view = OperonView(
		style=spec.mode,
		row_h=int(pick("row_height", 20 if classic else 26)),
		gene_h=int(pick("gene_height", 15 if classic else 8)),
		font=int(pick("font_size", 12 if classic else 13)),
		pad=int(pick("pad", 16)),
		domain_h=int(pick("domain_height", 6)),
		bp_per_px=pick("bases_per_pixel", 10.4) or 10.4,
		features=spec.features, monochrome=spec.monochrome,
		show_numbers=spec.family_numbers)
	view.LABEL_STEP = int(pick("label_step", OperonView.LABEL_STEP))
	view.SECRETION_BAND_OPACITY = pick("band_opacity",
									   OperonView.SECRETION_BAND_OPACITY)
	if classic:
		view.CLASSIC = dict(OperonView.CLASSIC,
							head=pick("arrow_head", OperonView.CLASSIC["head"]),
							min_w=pick("min_gene_width", OperonView.CLASSIC["min_w"]))
	view.numbers = data.numbers
	if spec.tree_width and data.newick:
		view.tree_w = int(320 * spec.tree_width)
		view.newick = data.newick
	order = data.order
	return view.render(data.genes, families, data.species, order,
					   domains=domains, clans=data.clans, features=features,
					   labels=data.labels, secretion=secretion)


def _feature_kind(kind: str) -> str:
	low = (kind or "").lower()
	if "signal" in low:
		return "signalp"
	return "tmhmm"


def render_all(specs: List[FigureSpec], data: RunData, out_path, verbose=False,
			   pdf=False):
	"""Render every spec. Returns the list of figure names actually written."""
	written = []
	for spec in specs:
		try:
			svg = render_figure(spec, data)
		except Exception as e:
			print("Warning: could not draw figure {!r} ({}).".format(spec.name, e))
			continue
		if not svg:
			if verbose:
				print(">> figure {!r} skipped: nothing to draw".format(spec.name))
			continue
		with open(out_path("_{}.svg".format(spec.name)), "w") as out:
			out.write(svg)
		written.append(spec.name)
	if pdf and written:
		import flags_pdf
		flags_pdf.convert_all([out_path("_{}.svg".format(n)) for n in written])
	return written


def main():
    parser = argparse.ArgumentParser(
        description="Redraw the figures of a finished FlaGs3 run.")
    parser.add_argument("-d", "--data", required=True,
                        help="Output directory of a finished FlaGs3 run.")
    parser.add_argument("-f", "--format",
                        help="Figure table (TSV). Default: the "
                             "{} inside --data, else the built-in "
                             "set.".format(TABLE_NAME))
    parser.add_argument("-o", "--output",
                        help="Where to write the  Default: --data, "
                             "which overwrites the existing ones.")
    parser.add_argument("--prefix",
                        help="File prefix inside --data. Default: the "
                             "directory's own name, which is what FlaGs3 uses.")
    parser.add_argument("--write_table", action="store_true",
                        help="Write a starting {} into --data and "
                             "exit, so you have something to edit.".format(
                                 TABLE_NAME))
    parser.add_argument("-pdf", "--pdf", action="store_true",
                        help="Also write a PDF beside each SVG. Needs one of "
                             "cairosvg, svglib, rsvg-convert or inkscape.")
    parser.add_argument("-vb", "--verbose", action="store_true")
    args = parser.parse_args()

    if not os.path.isdir(args.data):
        sys.exit("Error: --data directory not found: {}".format(args.data))

    if args.write_table:
        target = os.path.join(args.data, TABLE_NAME)
        write_default_table(target)
        print("wrote {}".format(target))
        return

    table = args.format
    if table is None:
        table = os.path.join(args.data, TABLE_NAME)
        if not os.path.isfile(table):
            write_default_table(table)
    if not os.path.isfile(table):
        sys.exit("Error: --format table not found: {}".format(table))

    try:
        specs = read_table(table)
    except ValueError as e:
        sys.exit("Error in {}: {}".format(table, e))

    try:
        data = load_run(args.data, args.prefix)
    except (OSError, ValueError) as e:
        sys.exit("Error: {}".format(e))

    out_dir = os.path.abspath(args.output or args.data)
    os.makedirs(out_dir, exist_ok=True)
    prefix = args.prefix or os.path.basename(os.path.normpath(
        os.path.abspath(args.data)))
    out_path = lambda suffix: os.path.join(out_dir, prefix + suffix)

    if args.verbose:
        print(">> {} rows, {} families, {} figures from {}".format(
            len({g.query for g in data.genes}), len(data.families), len(specs),
            os.path.basename(table)))

    written = render_all(specs, data, out_path, verbose=args.verbose,
                         pdf=args.pdf)
    if not written:
        sys.exit("Error: nothing was drawn. Check the figure table and that "
                 "--data holds the tables the figures need.")
    print("figures in {}/".format(
        os.path.relpath(out_dir) if out_dir.startswith(os.getcwd() + os.sep)
        else out_dir))
    for name in written:
        print("  {}_{}.svg".format(prefix, name))
        if args.pdf and os.path.isfile(out_path("_{}.pdf".format(name))):
            print("  {}_{}.pdf".format(prefix, name))


if __name__ == "__main__":
    main()