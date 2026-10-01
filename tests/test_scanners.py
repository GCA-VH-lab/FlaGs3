import stat
from pathlib import Path

import pytest

from flags3 import features
from flags3.run import RunDir
from flags3.schema import Annotation
from flags3.stage import Runner
from flags3.stages.defence import Defence, DefenceCall, ToolStatus
from flags3.stages.extract import Extract
from flags3.stages.features import Feature, Features
from flags3.stages.fetch import Fetch
from flags3.stages.genomad import Genomad, database, plasmid_name, virus_name
from tests.synth import write_genome

FAKE_GENOMAD = """#!/bin/sh
out="$3"; mkdir -p "$out/x_summary"
printf 'seq_name\\tlength\\ttopology\\tcoordinates\\tn_genes\\tvirus_score\\ttaxonomy\\n' > "$out/x_summary/x_virus_summary.tsv"
printf 'w0|provirus_1\\t300\\tProvirus\\t101-400\\t3\\t0.97\\tViruses;Duplodnaviria;Caudoviricetes\\n' >> "$out/x_summary/x_virus_summary.tsv"
printf 'seq_name\\tlength\\tplasmid_score\\tconjugation_genes\\tamr_genes\\n' > "$out/x_summary/x_plasmid_summary.tsv"
printf 'w1\\t2300\\t0.8\\tNA\\tblaTEM\\n' >> "$out/x_summary/x_plasmid_summary.tsv"
"""
FAKE_DF = """#!/bin/sh
out=""; while [ $# -gt 0 ]; do case "$1" in -o) out="$2"; shift;; esac; shift; done
mkdir -p "$out"
printf 'sys_id\\ttype\\tsubtype\\tprotein_in_syst\\n' > "$out/neighbourhoods_defense_finder_systems.tsv"
printf 'r0_1\\tRM\\tRM_Type_I\\tr0_1,r0_2\\n' >> "$out/neighbourhoods_defense_finder_systems.tsv"
"""
FAKE_PADLOC = """#!/bin/sh
out=""; while [ $# -gt 0 ]; do case "$1" in --outdir) out="$2"; shift;; esac; shift; done
mkdir -p "$out"
printf 'system.number,seqid,system,target.name\\n' > "$out/neighbourhoods_padloc.csv"
printf '1,r0,RM_Type_I,r0_1\\n1,r0,RM_Type_I,r0_2\\n2,r1,CBASS,r1_1\\n' >> "$out/neighbourhoods_padloc.csv"
"""
FAKE_TMHMM = """#!/sh
"""
FAKE_TM = """#!/bin/sh
out=""; while [ $# -gt 0 ]; do case "$1" in --output-dir) out="$2"; shift;; esac; shift; done
mkdir -p "$out"
printf '>WP_004 | TM\\nMKV\\nSSSSSMMMMMOOOOOMMMMMIII\\n>WP_005 | GLOB\\nMK\\nOOOO\\n' > "$out/predicted_topologies.3line"
"""
FAKE_SP = """#!/bin/sh
out=""; while [ $# -gt 0 ]; do case "$1" in --output_dir) out="$2"; shift;; esac; shift; done
mkdir -p "$out"
printf '# SignalP-6.0\\nWP_004\\tSP\\t0.01\\t0.99\\tCS pos: 22-23. Pr: 0.9\\nWP_005\\tOTHER\\t0.99\\t0.01\\n' > "$out/prediction_results.txt"
"""


def _script(tmp_path, name, body):
	path = tmp_path / name
	path.write_text(body)
	path.chmod(path.stat().st_mode | stat.S_IEXEC)
	return path


def _tools(tmp_path):
	table = tmp_path / "tools.tsv"
	table.write_text("#name\tcommand\tdirectory\tscan_range\tengine\toptions\n"
		"genomad\t{} end-to-end {{in}} {{out}} {{db}}\t\t\t\tdb={}\n"
		"defensefinder\t{} run -o {{out}} {{faa}}\n"
		"padloc\t{} --faa {{faa}} --gff {{gff}} --outdir {{out}}\n"
		"deeptmhmm\t{} --fasta {{fasta}} --output-dir {{out}}\n"
		"signalp\t{} --fastafile {{fasta}} --output_dir {{out}}\n".format(
			_script(tmp_path, "genomad", FAKE_GENOMAD), tmp_path / "gdb",
			_script(tmp_path, "defense-finder", FAKE_DF), _script(tmp_path, "padloc", FAKE_PADLOC),
			_script(tmp_path, "deeptmhmm", FAKE_TM), _script(tmp_path, "signalp6", FAKE_SP)))
	(tmp_path / "gdb").mkdir()
	(tmp_path / "gdb" / "version.txt").write_text("1")
	return table


def _run(tmp_path, queries, **config):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("\n".join(queries) + "\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "fetch.slots": "genome",
		"tools": str(_tools(tmp_path)), **config})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	runner.execute(Fetch())
	runner.execute(Extract())
	return run, runner


def test_genomad_names_and_database(tmp_path):
	assert virus_name("Viruses;Duplodnaviria;Caudoviricetes") == "Caudoviricetes"
	assert virus_name("Viruses;unclassified") == "virus"
	assert plasmid_name("NA", "blaTEM") == "plasmid (AMR)"
	assert plasmid_name("3", "0") == "conjugative plasmid"
	from flags3.tools import Tool
	t = Tool("genomad")
	with pytest.raises(Exception):
		database(t, "")
	(tmp_path / "genomad_db").mkdir()
	(tmp_path / "genomad_db" / "version.txt").write_text("1")
	assert database(t, str(tmp_path)) == tmp_path / "genomad_db"


def test_genomad_stage(tmp_path):
	run, runner = _run(tmp_path, ["WP_004", "WP_009"], scan_range=1000, scan_margin=100, gene=1, genomad=True)
	assert runner.execute(Genomad())
	ann = sorted(Annotation.read(run.stage_file("genomad", Annotation.FILE)), key=lambda a: a.category)
	assert [(a.subject, a.category, a.label, a.start, a.end, a.kind) for a in ann] == [
		("GCF_TEST|B", "plasmid", "plasmid (AMR)", 1, 2300, "band"),
		("GCF_TEST|A", "virus", "Caudoviricetes", 1000, 1299, "band")]
	assert ann[1].score == 0.97
	lines = run.stage_file("genomad", "mobile_elements.tsv").read_text().splitlines()
	assert "genomad_taxonomy" in lines[0]


def test_defence_stage_merges_tools(tmp_path):
	run, runner = _run(tmp_path, ["WP_004", "WP_009"], gene=1, defensefinder=True, padloc=True)
	assert runner.execute(Defence())
	gff = run.stage_file("defence", "raw/neighbourhoods.gff").read_text()
	assert "##sequence-region r0 2000 3500" in gff and "pseudogene" not in gff
	calls = DefenceCall.read(run.stage_file("defence", DefenceCall.FILE))
	assert len(calls) == 2
	rm = next(c for c in calls if c.type == "RM_Type_I")
	assert rm.called_by == "defensefinder,padloc"
	assert rm.genes == "WP_004,WP_005" and (rm.start, rm.end) == (2000, 3500)
	assert rm.rows == "WP_004|GCF_TEST"
	cbass = next(c for c in calls if c.type == "CBASS")
	assert cbass.called_by == "padloc" and cbass.genes == "WP_008"
	ann = Annotation.read(run.stage_file("defence", Annotation.FILE))
	assert {(a.category, a.tool) for a in ann} == {("RM_Type_I", "defensefinder,padloc"), ("CBASS", "padloc")}
	assert all(a.kind == "band" and a.space == "bp" for a in ann)
	status = {s.tool: s for s in ToolStatus.read(run.stage_file("defence", ToolStatus.FILE))}
	assert status["defensefinder"].systems == 1 and status["padloc"].systems == 2


def test_defence_missing_tool_is_reported(tmp_path):
	run, runner = _run(tmp_path, ["WP_004"], gene=1, defensefinder=True)
	cfg = run.config()
	cfg.set("tools", "")
	cfg.save()
	runner.config = cfg
	assert runner.execute(Defence()) is False
	assert "not run" in ToolStatus.read(run.stage_file("defence", ToolStatus.FILE))[0].status


def test_feature_parsers():
	regions = features.parse_3line(">p1 | TM\nMKV\nSSSMMMOOMMII\n>p2\nMK\nOOOO\n", want_signal=True)
	assert [(r.protein, r.kind, r.start, r.end) for r in regions] == [("p1", "tm", 4, 6), ("p1", "tm", 9, 10), ("p1", "signal", 1, 3)]
	assert features.parse_3line(">p1\nM\nSSMM\n", want_signal=False) == [features.Region("p1", "tm", 3, 4)]
	sp = features.parse_signalp("# hdr\np1\tSP\t0.1\t0.9\tCS pos: 22-23. Pr: 0.9\np2\tOTHER\t1\t0\n")
	assert sp == [features.Region("p1", "signal", 1, 22)]


def test_features_stage_local(tmp_path):
	run, runner = _run(tmp_path, ["WP_004"], gene=1, tmhmm=True, local_tmhmm=True, signalp=True, local_signalp=True)
	assert runner.execute(Features())
	feats = Feature.read(run.stage_file("features", Feature.FILE))
	assert [(f.protein, f.kind, f.start, f.end, f.tool) for f in feats] == [
		("WP_004", "signal", 1, 22, "signalp"), ("WP_004", "tm", 6, 10, "deeptmhmm"), ("WP_004", "tm", 16, 20, "deeptmhmm")]
	ann = Annotation.read(run.stage_file("features", Annotation.FILE))
	assert {(a.kind, a.category, a.label) for a in ann} == {("triangle", "signal", "SP"), ("segment", "tm", "TM")}
	assert all(a.space == "aa" for a in ann)


def test_features_without_biolib_fails_softly(tmp_path, monkeypatch):
	run, runner = _run(tmp_path, ["WP_004"], gene=1, tmhmm=True)
	import builtins
	real_import = builtins.__import__
	def no_biolib(name, *a, **k):
		if name == "biolib":
			raise ImportError
		return real_import(name, *a, **k)
	monkeypatch.setattr(builtins, "__import__", no_biolib)
	assert runner.execute(Features()) is False
	assert "pybiolib" in run.status("features").get("error")
