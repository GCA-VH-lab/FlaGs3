import re
from pathlib import Path

import pytest

from flags3 import fasta
from flags3.render import palettes
from flags3.render.data import RunData
from flags3.render.figure import Figure, parts_for
from flags3.render.style import Colours, FigureSpec, StyleError, parse_row, read_overrides, read_table
from flags3.run import RunDir
from flags3.schema import Annotation
from flags3.stage import Runner
from flags3.stages.cluster import Cluster
from flags3.stages.extract import Extract
from flags3.stages.fetch import Fetch
from flags3.stages.figures import Figures
from tests.synth import write_genome


def _run(tmp_path, queries=("WP_004", "WP_009", "WP_001"), **config):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("\n".join(queries) + "\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "gene": 2,
		"cluster_method": "jackhmmer", "cpu": 1, **config})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	for stage in (Fetch(), Extract(), Cluster()):
		assert runner.execute(stage)
	return run, runner


def _fake_stage(run, name, annotations):
	out = run.reset(name)
	Annotation.write(out / Annotation.FILE, annotations)
	status = run.status(name)
	status.set("status", "ok")
	status.save()


def test_palettes_and_table():
	for name, fn in palettes.PALETTES.items():
		assert len(fn(5)) == 5 and all(c.startswith("#") for c in fn(5))
	assert palettes.colourblind(3) == list(palettes.OKABE_ITO[:3])
	assert palettes.monochrome(2) == [palettes.GREY, palettes.GREY]
	specs = {s.name: s for s in read_table(None)}
	assert specs["sismis"].monochrome and specs["sismis"].palette == "bright"
	assert specs["classic"].classic and specs["tree"].tree
	assert specs["classic"].value("gene_height") == 15 and specs["neighbors"].value("gene_height") == 8
	with pytest.raises(StyleError):
		parse_row({"name": "x", "layers": "cluster", "palette": "neon"}, "t")
	with pytest.raises(StyleError):
		parse_row({"name": "x", "layers": "cluster", "row_height": "big"}, "t")


def test_colours_and_overrides(tmp_path):
	ann = [Annotation("a", "-", None, None, "fill", "family:1", "1", "t", None),
		Annotation("b", "-", None, None, "fill", "family:2", "2", "t", None)]
	c = Colours("bright", {("cluster", "family:2"): "#123456"})
	table = c.assign("cluster", ann)
	assert table["family:2"] == "#123456" and table["family:1"] == palettes.bright(2)[0]
	assert c.of("cluster", "family:9") == palettes.GREY
	mono = Colours("bright", monochrome=True)
	assert set(mono.assign("cluster", ann).values()) == {palettes.GREY}
	assert mono.assign("sismis", ann)["family:1"] != palettes.GREY
	path = tmp_path / "colours.tsv"
	path.write_text("#stage\tcategory\tcolour\ncluster\tfamily:1\t#ff0000\n")
	assert read_overrides(path) == {("cluster", "family:1"): "#ff0000"}


def test_layout_flips_minus_strand_query(tmp_path):
	run, _ = _run(tmp_path)
	data = RunData.load(run)
	assert data.order == ["WP_004|GCF_TEST", "WP_009|GCF_TEST", "WP_001|GCF_TEST"]
	fig = Figure(data, FigureSpec("n", ["cluster"]), Colours("bright"))
	L = fig.layout
	assert L.reversed["WP_004|GCF_TEST"] and not L.reversed["WP_001|GCF_TEST"]
	query = next(g for g in data.genes["WP_004|GCF_TEST"] if g.offset == 0)
	arrow = L.arrow("WP_004|GCF_TEST", query)
	assert arrow.strand == "+" and arrow.points[2][0] == pytest.approx(arrow.x1)
	assert L.x("WP_004|GCF_TEST", query.end) == pytest.approx(L.center)
	right_gene = next(g for g in data.genes["WP_004|GCF_TEST"] if g.offset == -1)
	assert L.x("WP_004|GCF_TEST", right_gene.start) < L.center
	assert L.x_in_gene(arrow, query, 1) == pytest.approx(arrow.x0, abs=1)
	assert L.x_in_gene(arrow, query, 200) == pytest.approx(arrow.x1)
	assert L.y("WP_009|GCF_TEST") - L.y("WP_004|GCF_TEST") == L.row_pitch


def test_figure_layers_and_legend(tmp_path):
	run, _ = _run(tmp_path)
	_fake_stage(run, "domains", [Annotation("WP_004", "aa", 10, 100, "wedge", "PF1", "Dom_1", "pfam", 1e-5)])
	_fake_stage(run, "features", [Annotation("WP_005", "aa", 1, 20, "triangle", "signal", "SP", "signalp", None),
		Annotation("WP_005", "aa", 30, 50, "segment", "tm", "TM", "deeptmhmm", None)])
	_fake_stage(run, "sismis", [Annotation("GCF_TEST|A", "bp", 1000, 3000, "band", "T3SS", "T3SS", "sismis", 0.9)])
	_fake_stage(run, "defence", [Annotation("GCF_TEST|A", "bp", 2000, 2600, "band", "RM", "RM", "padloc", None)])
	data = RunData.load(run)
	spec = FigureSpec("all", ["cluster", "domains", "features", "sismis", "genomad", "defence"], numbers=True)
	svg = Figure(data, spec, Colours("bright")).render()
	for layer in ("layer-genes", "layer-domains", "layer-features", "layer-sismis", "layer-defence", "layer-labels", "layer-legend", "layer-row-labels", "layer-outlines"):
		assert 'id="{}"'.format(layer) in svg, layer
	assert "layer-genomad" not in svg
	assert svg.index('id="layer-sismis"') < svg.index('id="layer-genes"') < svg.index('id="layer-domains"')
	assert 'id="layer-band-codes"' in svg
	assert ">S1<" in svg and ">D1<" in svg and "S1. T3SS" in svg and "D1. RM" in svg and "1. PF1 (Dom_1)" in svg
	assert "Transmembrane region" in svg and "Signal peptide" in svg and "Query protein" in svg
	assert 'clip-path="url(#clip' in svg
	assert svg.count("<svg") == 1 and svg.strip().endswith("</svg>")


def test_classic_and_parts(tmp_path):
	run, _ = _run(tmp_path)
	data = RunData.load(run)
	classic = Figure(data, FigureSpec("c", ["cluster"], mode="classic", palette="classic"), Colours("classic")).render()
	assert 'fill="#000000"' in classic
	parts = parts_for(data, FigureSpec("n", ["cluster"]), Colours("bright"), max_height=260, no_overlaps=False)
	assert len(parts) == 2 and [len(p.rows) for p in parts] == [2, 1]
	assert len(parts_for(data, FigureSpec("n", ["cluster"]), Colours("bright"), max_height=5000, no_overlaps=False)) == 1


def test_figures_stage(tmp_path):
	run, runner = _run(tmp_path)
	assert runner.execute(Figures())
	names = sorted(p.name for p in run.stage_dir("figures").glob("*.svg"))
	assert names == ["classic.svg", "neighbors.svg"]
	cfg = run.config()
	cfg.set("no_figures", True)
	cfg.save()
	assert not Figures().wanted(cfg)
	table = tmp_path / "figs.tsv"
	table.write_text("#name\tlayers\tmode\ttree\tnumbers\tpalette\nmine\tcluster\tversatile\t0\ttrue\tcolourblind\n")
	cfg.set("no_figures", False)
	cfg.set("figures", str(table))
	cfg.save()
	runner.config = cfg
	assert runner.execute(Figures())
	assert sorted(p.name for p in run.stage_dir("figures").glob("*.svg")) == ["mine.svg"]
