from pathlib import Path

import pytest

from flags3 import cluster
from flags3.run import RunDir
from flags3.schema import Annotation, ClusterHit, Family, Gene, MISSING, RowInfo
from flags3.stage import Runner
from flags3.stages.cluster import Cluster, ClusterRna, FamilyTables
from flags3.stages.extract import Extract
from flags3.stages.fetch import Fetch
from flags3.tools import Tool, Tools
from tests.synth import write_genome


def test_connected_components_are_symmetric_and_sorted():
	adjacency = {"a": {"b"}, "b": set(), "c": {"a"}, "d": set(), "e": {"f"}, "f": {"e"}}
	assert cluster.connected_components(adjacency) == [["a", "b", "c"], ["e", "f"], ["d"]]


def test_by_name_groups_products():
	fams = cluster.by_name({"r1": "tRNA-Ala", "r2": " trna-ala ", "r3": "16S rRNA", "r4": ""})
	assert fams == [["r1", "r2"], ["r3"], ["r4"]]


def test_method_lookup_and_options():
	tools = Tools.load()
	assert cluster.methods(tools) == ["jackhmmer", "mmseqs_cluster", "mmseqs_cluster_exhaustive", "nhmmer"]
	assert cluster.option(tools["jackhmmer"], "iterations", 1) == 3
	assert cluster.option(tools["jackhmmer"], "incE", 1.0) == 1e-3
	with pytest.raises(cluster.ClusterError):
		cluster.method(tools, "mafft")
	bad = Tool("x")
	bad.options = {"chunk": "many"}
	with pytest.raises(cluster.ClusterError):
		cluster.option(bad, "chunk", 100)


def test_mmseqs_command_and_pairs(tmp_path):
	tools = Tools.load()
	c = cluster.build(tools, "mmseqs_cluster", 2, tmp_path / "raw")
	argv = c.command_line(Path("in.fa"), Path("hits.tsv"), 10)
	assert argv[:3] == ["mmseqs", "easy-search", "in.fa"]
	assert "--threads" in argv and argv[argv.index("--threads") + 1] == "2"
	assert argv[argv.index("--max-seqs") + 1] == "300"
	assert argv[argv.index("-s") + 1] == "7.5"
	hits = tmp_path / "hits.tsv"
	hits.write_text("a\tb\nb\tzzz\nc\tc\n")
	assert cluster.MmseqsClusterer.read_pairs(hits, {"a": "", "b": "", "c": ""}) == {"a": {"b"}, "b": set(), "c": {"c"}}
	broken = Tool("m")
	broken.engine = "mmseqs"
	broken.command = "mmseqs {in} {undefined}"
	with pytest.raises(cluster.ClusterError):
		cluster.MmseqsClusterer(broken, 1, tmp_path).command_line(Path("a"), Path("b"), 1)


def test_family_labels(tmp_path):
	families = [["a", "b"], ["q", "x"], ["s"], ["t"]]
	occurrences = {"a": 1, "b": 1, "q": 1, "x": 1, "s": 3, "t": 1}
	tables = FamilyTables(families, {"a": {"b"}}, occurrences, {"q"}, "", "jackhmmer")
	assert tables.labels == {2: "1", 0: "2", 1: "Q1"}
	tables.write(tmp_path)
	fams = Family.read(tmp_path / Family.FILE)
	assert [(f.family, f.label, f.size, f.occurrences) for f in fams] == [(1, "2", 2, 2), (2, "Q1", 2, 2), (3, "1", 1, 3), (4, MISSING, 1, 1)]
	ann = Annotation.read(tmp_path / Annotation.FILE)
	assert {a.subject for a in ann} == {"a", "b", "q", "x", "s"}
	assert all(a.kind == "fill" and a.whole for a in ann)
	assert {a.category for a in ann if a.subject == "s"} == {"family:3"}
	hits = {h.accession: h for h in ClusterHit.read(tmp_path / ClusterHit.FILE)}
	assert hits["a"].hits == "b" and hits["b"].hits == MISSING and hits["t"].family == 4
	rna = FamilyTables([["r1", "r2"]], {}, {"r1": 1, "r2": 1}, set(), "R", "nhmmer")
	assert rna.labels == {0: "R1"}


def test_jackhmmer_stage_on_synthetic_genome(tmp_path):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("WP_004\nWP_009\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "fetch.slots": "rna,genome",
		"gene": 4, "cluster_method": "jackhmmer", "rna_method": "nhmmer", "cluster_rna": True, "cpu": 1})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	for stage in (Fetch(), Extract(), Cluster(), ClusterRna()):
		assert runner.execute(stage)
	fams = Family.read(run.stage_file("cluster", Family.FILE))
	members = {a for f in fams for a in f.accessions}
	assert members == {"WP_001", "WP_002", "WP_004", "WP_005", "WP_006", "WP_007", "WP_008", "WP_009"}
	assert all(f.label == MISSING for f in fams)
	assert Annotation.read(run.stage_file("cluster", Annotation.FILE)) == []
	rna = Family.read(run.stage_file("cluster_rna", Family.FILE))
	assert [f.members for f in rna] == ["rna_001"]
