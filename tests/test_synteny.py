import random
import shutil

import pytest

from flags3 import fasta
from flags3.cli import build_parser, _imply
from flags3.run import RunDir
from flags3.schema import GenomeFiles, Window
from flags3.stage import Runner
from flags3.stages import PIPELINE
from flags3.stages.synteny import Link, Pair, Synteny
from flags3.stages.tree import Leaf, Tree
from flags3.tools import Tools

needs_blastn = pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn is not installed")


def _sequence(seed: int, length: int) -> str:
	rng = random.Random(seed)
	return "".join(rng.choice("ACGT") for _ in range(length))


def _ok(run, name):
	out = run.reset(name)
	status = run.status(name)
	status.set("status", "ok")
	status.save()
	return out


def _run(tmp_path, genomes: dict[str, str], windows: list[Window]):
	run = RunDir(tmp_path / "out").create("test", "")
	paths = {}
	for assembly, sequence in genomes.items():
		paths[assembly] = tmp_path / (assembly + ".fna")
		fasta.write(paths[assembly], [("c", sequence)])
	GenomeFiles.write(_ok(run, "fetch") / GenomeFiles.FILE,
		[GenomeFiles(a, "local", "-", "-", "-", str(paths[a])) if a in paths else GenomeFiles(a, "local", "-", "-", "-", "-")
			for a in sorted({w.assembly for w in windows})])
	Window.write(_ok(run, "extract") / Window.FILE, windows)
	Leaf.write(_ok(run, "tree") / Leaf.FILE, [Leaf(i + 1, w.row_id) for i, w in enumerate(windows)])
	cfg = run.config()
	cfg.update({"synteny": True, "cpu": 2})
	cfg.save()
	return run, cfg


def _window(row: str, assembly: str, lo: int, hi: int, length: int) -> Window:
	return Window(row, assembly, "c", length, lo, hi, 1, length, 1, length)


def test_flag_implies_tree_order_and_stage_order():
	args = build_parser().parse_args(["run", "-i", "x", "-sy"])
	_imply(args)
	assert args.synteny and args.tree_order and args.tree
	assert PIPELINE.index(Synteny) == PIPELINE.index(Tree) + 1
	assert Synteny.needs == ("genome",) and Synteny.optional
	tool = Tools.load().tools["synteny"]
	assert tool.argv(**{"in": "a", "subject": "b", "out": "o", "evalue": "1e-5"})[:7] == [
		"blastn", "-task", "blastn", "-query", "a", "-subject", "b"]


def test_pair_filters_and_places_hits(tmp_path):
	table = tmp_path / "1.tsv"
	table.write_text("11\t510\t1\t500\t95.5\t500\n1\t80\t300\t221\t99.0\t80\n600\t900\t900\t600\t50.0\t301\n100\t400\t700\t400\t88.0\t301\n")
	pair = Pair(_window("a", "A", 1001, 3000, 5000), _window("b", "B", 201, 2000, 5000), tmp_path, tmp_path, table)
	links = pair.links(0.6, 100)
	assert links == [Link("a", 1011, 1510, "b", 201, 700, "+", 95.5, 500), Link("a", 1100, 1400, "b", 600, 900, "-", 88.0, 301)]


@needs_blastn
def test_stage_compares_only_neighbours_in_tree_order(tmp_path):
	shared, other = _sequence(1, 1500), _sequence(2, 1200)
	genomes = {
		"A": _sequence(10, 1000) + shared + _sequence(11, 500),
		"B": _sequence(20, 300) + fasta.reverse_complement(shared) + other + _sequence(21, 200),
		"C": _sequence(30, 600) + other + _sequence(31, 700),
		"D": shared + _sequence(40, 400),
	}
	windows = [_window("a|A", "A", 801, 2700, 3000), _window("b|B", "B", 1, 3200, 3200),
		_window("c|C", "C", 401, 2000, 2500), _window("d|D", "D", 1, 1900, 1900), _window("e|E", "E", 1, 100, 100)]
	run, cfg = _run(tmp_path, genomes, windows)
	assert Synteny().wanted(cfg)
	assert Runner(run, cfg, report=lambda m: None).execute(Synteny())
	links = Link.read(run.stage_file("synteny", Link.FILE))
	assert {(l.row_a, l.row_b) for l in links} == {("a|A", "b|B"), ("b|B", "c|C")}
	first = next(l for l in links if l.row_a == "a|A")
	assert (first.start_a, first.end_a, first.start_b, first.end_b) == pytest.approx((1001, 2500, 301, 1800), abs=3)
	assert first.strand == "-" and first.identity > 99.0 and first.length >= 1500
	second = next(l for l in links if l.row_a == "b|B")
	assert (second.start_a, second.end_a, second.start_b, second.end_b) == pytest.approx((1801, 3000, 601, 1800), abs=3)
	assert second.strand == "+"
	assert sorted(p.name for p in run.stage_dir("synteny").glob("raw/*.tsv")) == ["1.tsv", "2.tsv", "3.tsv"]


@needs_blastn
def test_stage_fails_softly_without_the_tool(tmp_path):
	table = tmp_path / "tools.tsv"
	table.write_text("#name\tcommand\tdirectory\tscan_range\tengine\toptions\nsynteny\t/nonexistent/blastn -query {in} -subject {subject} -out {out}\t\t\t\t\n")
	run, cfg = _run(tmp_path, {"A": _sequence(1, 500), "B": _sequence(2, 500)},
		[_window("a|A", "A", 1, 500, 500), _window("b|B", "B", 1, 500, 500)])
	cfg.set("tools", str(table))
	cfg.save()
	assert not Runner(run, cfg, report=lambda m: None).execute(Synteny())
	assert "could not run" in run.status("synteny").get("error")
