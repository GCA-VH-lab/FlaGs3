import os
import stat
from pathlib import Path

import pytest

from flags3 import tree
from flags3.run import RunDir
from flags3.stage import Runner
from flags3.stages.tree import Leaf, Tree
from flags3.tools import Tools

FAKE_MAFFT = """#!/bin/sh
python3 - "$1" <<'PY'
import sys
names, seqs, n = [], [], None
for line in open(sys.argv[1]):
	if line.startswith('>'):
		names.append(line[1:].strip()); seqs.append('')
	else:
		seqs[-1] += line.strip()
width = max(len(s) for s in seqs)
for n, s in zip(names, seqs):
	print('>' + n); print(s + '-' * (width - len(s)))
PY
"""
FAKE_TRIMAL = """#!/bin/sh
cp "$2" "$4"
"""
FAKE_TREE = """#!/bin/sh
names=$(grep '^>' "$1" | sed 's/^>//')
set -- $names
echo "(($1:0.1,$2:0.2):0.3,$3:0.5);"
"""


def _fake_tools(tmp_path, with_trimal=True):
	bin_dir = tmp_path / "bin"
	bin_dir.mkdir()
	scripts = {"mafft": FAKE_MAFFT, "VeryFastTree": FAKE_TREE}
	if with_trimal:
		scripts["trimal"] = FAKE_TRIMAL
	for name, body in scripts.items():
		path = bin_dir / name
		path.write_text(body)
		path.chmod(path.stat().st_mode | stat.S_IEXEC)
	table = tmp_path / "tools.tsv"
	lines = ["#name\tcommand",
		"mafft\t{} {{in}}".format(bin_dir / "mafft"),
		"veryfasttree\t{} {{in}}".format(bin_dir / "VeryFastTree"),
		"trimal\t{} -in {{in}} -out {{out}} -fasta {{mode}}".format(bin_dir / ("trimal" if with_trimal else "nope"))]
	table.write_text("\n".join(lines) + "\n")
	return table


def test_builder_runs_the_three_tools(tmp_path):
	tools = Tools.load(_fake_tools(tmp_path))
	builder = tree.TreeBuilder(tools, tmp_path / "work")
	newick = builder.build({"a|x": "MKV", "b|y": "MKVLL", "c|z": "MK"})
	assert newick == "((a|x:0.1,b|y:0.2):0.3,c|z:0.5);"
	assert builder.raw_alignment == {"a|x": "MKV--", "b|y": "MKVLL", "c|z": "MK---"}
	assert builder.alignment == builder.raw_alignment
	assert len(builder.commands) == 3 and "-gt 0.1" in builder.commands[1]
	with pytest.raises(tree.TreeError):
		builder.build({"a": "M", "b": "M"})


def test_internal_trim_when_trimal_missing(tmp_path, capsys):
	tools = Tools.load(_fake_tools(tmp_path, with_trimal=False))
	builder = tree.TreeBuilder(tools, tmp_path / "work", trimal_value=0.5)
	builder.build({"a": "MKV", "b": "MKVLL", "c": "MK"})
	assert builder.alignment == {"a": "MKV", "b": "MKV", "c": "MK-"}
	assert "trimal not found" in capsys.readouterr().out
	assert len(builder.commands) == 2


def test_leaf_order_ladderizes():
	newick = "((a:0.1,b:0.2):0.3,(c:0.1,(d:0.1,e:0.1):0.1):0.4);"
	assert tree.leaf_order(newick, ladderize=False) == ["a", "b", "c", "d", "e"]
	order = tree.leaf_order(newick)
	assert sorted(order) == ["a", "b", "c", "d", "e"] and order != ["a", "b", "c", "d", "e"]


def test_tree_stage(tmp_path):
	table = _fake_tools(tmp_path)
	run = RunDir(tmp_path / "out").create("test", "")
	ext = run.reset("extract")
	(ext / "queries.faa").write_text(">q1|g1\nMKV\n>q2|g2\nMKVLL\n>q3|g3\nMK\n")
	status = run.status("extract")
	status.set("status", "ok")
	status.save()
	cfg = run.config()
	cfg.update({"tree": True, "tools": str(table)})
	cfg.save()
	stage = Tree()
	assert stage.wanted(cfg)
	assert Runner(run, cfg, report=lambda m: None).execute(stage)
	assert run.stage_file("tree", "tree.nwk").read_text().startswith("((q1|g1")
	leaves = Leaf.read(run.stage_file("tree", Leaf.FILE))
	assert {l.row_id for l in leaves} == {"q1|g1", "q2|g2", "q3|g3"}
	assert [l.position for l in leaves] == [1, 2, 3]
	assert run.stage_file("tree", "commands.txt").read_text().count("\n") == 3
	assert run.stage_file("tree", "raw/queries.fasta").is_file()


def test_tree_stage_with_two_queries_fails_softly(tmp_path):
	table = _fake_tools(tmp_path)
	run = RunDir(tmp_path / "out").create("test", "")
	ext = run.reset("extract")
	(ext / "queries.faa").write_text(">q1|g1\nMKV\n>q2|g2\nMKVLL\n")
	status = run.status("extract")
	status.set("status", "ok")
	status.save()
	cfg = run.config()
	cfg.update({"tree_order": True, "tools": str(table)})
	cfg.save()
	messages = []
	assert Runner(run, cfg, report=messages.append).execute(Tree()) is False
	assert "at least 3" in run.status("tree").get("error")
	assert not Tree().wanted(run.config()) is False
