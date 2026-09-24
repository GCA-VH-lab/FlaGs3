from pathlib import Path

from flags3.log import note
from flags3.render import pdf
from flags3.render.data import RunData
from flags3.render.figure import parts_for
from flags3.render.style import Colours, read_overrides, read_table
from flags3.stage import Stage

COLOURS_FILE = "colours.tsv"


class Figures(Stage):
	name = "figures"
	requires = ("extract",)
	optional = True
	barrier = True

	def wanted(self, config) -> bool:
		return not config.flag("no_figures")

	def run(self, run, config, out: Path) -> None:
		specs = read_table(config.path("figures"))
		data = RunData.load(run)
		overrides = read_overrides(run.path / COLOURS_FILE)
		max_height = config.integer("figure_height", 16383)
		written = []
		for spec in specs:
			extra = [s for s in spec.layers if s not in ("cluster", "cluster_rna")]
			if extra and not any(s in data.annotations for s in extra):
				note("{}: skipped, no {} annotations".format(spec.name, " or ".join(extra)))
				continue
			if spec.tree and not data.newick:
				note("{}: skipped, no tree".format(spec.name))
				continue
			colours = Colours(spec.palette, overrides, spec.monochrome)
			figures = parts_for(data, spec, colours, max_height, config.flag("no_overlaps"))
			for i, figure in enumerate(figures):
				name = spec.name if len(figures) == 1 else "{}_part{}".format(spec.name, i + 1)
				path = out / (name + ".svg")
				path.write_text(figure.render())
				written.append(path)
				if config.flag("pdf"):
					if pdf.available():
						pdf.convert(path, path.with_suffix(".pdf"))
					else:
						print("Warning: --pdf needs cairosvg (pip install cairosvg); SVG only.")
		note("{} figures written".format(len(written)))
		if not written:
			raise RuntimeError("no figure had anything to draw")
