import stat
from pathlib import Path

import pytest

from flags3 import fasta, windows
from flags3.run import RunDir
from flags3.scan import Diagnostic
from flags3.schema import Annotation, Window
from flags3.stage import Runner
from flags3.stages.extract import Extract
from flags3.stages.fetch import Fetch
from flags3.stages.sismis import Sismis, parse_clusters
from tests.synth import write_genome

FAKE_SISMIS = """#!/bin/sh
in=""; out=""
while [ $# -gt 0 ]; do case "$1" in -g) in="$2"; shift;; -o) out="$2"; shift;; esac; shift; done
mkdir -p "$out"
name=$(basename "$in" .fna)
printf 'sequence_id\\tcluster_id\\tstart\\tend\\tmax_p\\ttype\\tproteins\\n' > "$out/$name.clusters.tsv"
printf 'w0\\tw0_cluster_1\\t101\\t400\\t0.95\\tT3SS\\tw0_1;w0_2\\n' >> "$out/$name.clusters.tsv"
printf 'w0\\tw0_cluster_2\\t2\\t50\\t0.5\\tFlagellum\\tw0_3\\n' >> "$out/$name.clusters.tsv"
"""


def _tools(tmp_path):
	script = tmp_path / "sismis"
	script.write_text(FAKE_SISMIS)
	script.chmod(script.stat().st_mode | stat.S_IEXEC)
	table = tmp_path / "tools.tsv"
	table.write_text("#name\tcommand\nsismis\t{} run -g {{in}} -o {{out}}\n".format(script))
	return table


def _run(tmp_path, queries, **config):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("\n".join(queries) + "\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "fetch.slots": "genome",
		"sismis": True, "tools": str(_tools(tmp_path)), **config})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	runner.execute(Fetch())
	runner.execute(Extract())
	return run, runner


def test_merge_and_cut(tmp_path):
	run, _ = _run(tmp_path, ["WP_001", "WP_004", "WP_009"], scan_range=1000, scan_margin=200)
	rows = Window.read(run.stage_file("extract", Window.FILE))
	merged = windows.merge(rows)
	assert [(c.contig, c.cut_lo, c.cut_hi, c.scan_lo, c.scan_hi) for c in merged["GCF_TEST"]] == [
		("A", 1, 3800, 1, 3600), ("B", 1, 2400, 1, 2200)]
	assert merged["GCF_TEST"][0].rows == "WP_001|GCF_TEST,WP_004|GCF_TEST"
	batches = windows.Batches(tmp_path / "batches", max_bases=3000)
	paths = batches.write({"GCF_TEST": tmp_path / "genomes" / "GCF_TEST_genomic.fna"}, merged)
	assert [p.name for p in paths] == ["batch000.fna", "batch001.fna"]
	records = list(fasta.read(paths[0]))
	assert records[0][0] == "w0" and len(records[0][2]) == 3800
	assert "query=WP_001,WP_004" in records[0][1] and "genomic=1-3800" in records[0][1]
	assert batches.per_assembly() == {"GCF_TEST": (2, 3800 + 2400)}
	cut = batches.by_record()["w0"]
	assert cut.place(101, 400) == (101, 400, "full")
	assert cut.place(3500, 3700) == (3500, 3700, "partial")
	assert windows.Cut.read(tmp_path / "batches" / windows.Cut.FILE)[1].record == "w1"


def test_parse_clusters_requires_columns(tmp_path):
	bad = tmp_path / "x.clusters.tsv"
	bad.write_text("sequence_id\tstart\n")
	with pytest.raises(Exception):
		list(parse_clusters(bad))


def test_sismis_stage(tmp_path):
	run, runner = _run(tmp_path, ["WP_004", "WP_009"], scan_range=1000, scan_margin=100, gene=1)
	stage = Sismis()
	assert stage.needs == ("genome",) and stage.wanted(run.config())
	assert runner.execute(stage)
	ann = Annotation.read(run.stage_file("sismis", Annotation.FILE))
	assert len(ann) == 2
	first = next(a for a in ann if a.category == "T3SS")
	assert first.subject == "GCF_TEST|A" and first.space == "bp" and first.kind == "band"
	assert (first.start, first.end) == (1000, 1299) and first.score == 0.95 and first.tool == "sismis"
	lines = run.stage_file("sismis", "secretion.tsv").read_text().splitlines()
	assert lines[0].split("\t")[:8] == ["assembly", "contig", "start", "end", "type", "score", "coverage", "rows"]
	assert "sismis_cluster_id" in lines[0] and "sismis_proteins" in lines[0]
	t3ss = next(l for l in lines if "\tT3SS\t" in l).split("\t")
	assert t3ss[6] == "full" and t3ss[7] == "WP_004|GCF_TEST"
	flag = next(l for l in lines if "\tFlagellum\t" in l).split("\t")
	assert flag[6] == "partial" and flag[7] == "-"
	diag = {d.assembly: d for d in Diagnostic.read(run.stage_file("sismis", Diagnostic.FILE))}
	assert diag["GCF_TEST"].windows == 2 and diag["GCF_TEST"].hits == 2 and diag["GCF_TEST"].status == "ok"
	assert run.stage_file("sismis", "raw/batch000_sismis/batch000.clusters.tsv").is_file()


def test_missing_sismis_fails_softly(tmp_path):
	run, runner = _run(tmp_path, ["WP_004"], scan_range=1000)
	cfg = run.config()
	cfg.set("tools", "")
	cfg.save()
	runner.config = cfg
	assert runner.execute(Sismis()) is False
	assert "flags3 install sismis" in run.status("sismis").get("error")
