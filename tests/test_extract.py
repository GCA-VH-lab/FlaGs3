import pytest

from flags3 import fasta
from flags3.run import RunDir
from flags3.schema import Gene, MISSING, RangeReport, RowInfo, Unmatched, Window
from flags3.stage import Runner
from flags3.stages.extract import Extract, SCAN_GENES, SCAN_PROTEINS, PROTEINS, QUERIES, RNA
from flags3.stages.fetch import Fetch
from tests.synth import write_genome


def _run(tmp_path, queries, **config):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("\n".join(queries) + "\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "fetch.slots": "rna,genome", **config})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	runner.execute(Fetch())
	runner.execute(Extract())
	return run


def _genes(run):
	return Gene.read(run.stage_file("extract", Gene.FILE))


def test_flank_window_on_minus_strand_query(tmp_path):
	run = _run(tmp_path, ["WP_004"], gene=2)
	genes = _genes(run)
	assert [g.accession for g in genes] == ["rna_001", "pseudogene*", "WP_004", "WP_005", "WP_006"]
	assert [g.offset for g in genes] == [2, 1, 0, -1, -2]
	assert [g.strand for g in genes] == ["+", "+", "-", "+", "+"]
	assert genes[0].is_rna and not genes[1].is_rna
	assert all(g.row_id == "WP_004|GCF_TEST" for g in genes)
	row = RowInfo.read(run.stage_file("extract", RowInfo.FILE))[0]
	assert row.species == "Escherichia coli K-12"
	assert row.strand == "-" and row.contig == "A"
	window = Window.read(run.stage_file("extract", Window.FILE))[0]
	assert (window.lo, window.hi) == (1000, 4800)
	assert (window.scan_lo, window.scan_hi) == (1, 20000)
	assert (window.cut_lo, window.cut_hi) == (1, 20000)
	assert window.contig_length == 20000


def test_range_window_and_scan_span(tmp_path):
	run = _run(tmp_path, ["WP_004"], range=1200, scan_range=500, scan_margin=100)
	genes = _genes(run)
	assert [g.accession for g in genes] == ["WP_002", "rna_001", "pseudogene*", "WP_004", "WP_005"]
	window = Window.read(run.stage_file("extract", Window.FILE))[0]
	assert (window.scan_lo, window.scan_hi) == (1500, 3100)
	assert (window.cut_lo, window.cut_hi) == (1400, 3200)
	scan = Gene.read(run.stage_file("extract", SCAN_GENES))
	assert [g.accession for g in scan] == ["WP_004", "WP_005"]
	assert {n for n, _, _ in fasta.read(run.stage_file("extract", SCAN_PROTEINS))} == {"WP_004", "WP_005"}
	report = RangeReport.read(run.stage_file("extract", RangeReport.FILE))[0]
	assert report.contig_length == 20000
	assert (report.up_available, report.down_available) == (17400, 1999)
	assert (report.up_reached, report.down_reached) == (900, 1500)
	assert (report.genes_up, report.genes_down) == (1, 3)


def test_contig_edge_and_sequences(tmp_path):
	run = _run(tmp_path, ["WP_009", "WP_001"], gene=4)
	genes = _genes(run)
	by_row = {}
	for g in genes:
		by_row.setdefault(g.row_id, []).append(g.accession)
	assert by_row["WP_009|GCF_TEST"] == ["WP_008", "WP_009"]
	assert by_row["WP_001|GCF_TEST"] == ["WP_001", "WP_002", "rna_001", "pseudogene*", "WP_004"]
	proteins = {n: s for n, _, s in fasta.read(run.stage_file("extract", PROTEINS))}
	assert set(proteins) == {"WP_001", "WP_002", "WP_004", "WP_008", "WP_009"}
	assert "pseudogene*" not in proteins
	queries = dict((n, s) for n, _, s in fasta.read(run.stage_file("extract", QUERIES)))
	assert set(queries) == {"WP_009|GCF_TEST", "WP_001|GCF_TEST"}
	rna = {n: s for n, _, s in fasta.read(run.stage_file("extract", RNA))}
	assert rna["rna_001"].startswith("GC")
	assert not run.stage_file("extract", SCAN_GENES).exists()


def test_unmatched_queries_are_reported(tmp_path):
	run = _run(tmp_path, ["WP_004", "WP_999", "WP_001\tGCF_NOPE"])
	unmatched = Unmatched.read(run.stage_file("extract", Unmatched.FILE))
	assert {(u.query, u.assembly) for u in unmatched} == {("WP_999", MISSING), ("WP_001", "GCF_NOPE")}
	assert all(u.reason.startswith("unresolved") for u in unmatched)
	assert len(RowInfo.read(run.stage_file("extract", RowInfo.FILE))) == 1


def test_paired_assembly_matches_local_basename(tmp_path):
	run = _run(tmp_path, ["WP_005\tGCF_TEST"], gene=1)
	assert [g.accession for g in _genes(run)] == ["WP_004", "WP_005", "WP_006"]


def test_rerun_overwrites(tmp_path):
	run = _run(tmp_path, ["WP_004"], gene=1)
	assert len(_genes(run)) == 3
	cfg = run.config()
	cfg.set("gene", 3)
	cfg.save()
	Runner(run, cfg, report=lambda m: None).execute(Extract())
	assert len(_genes(run)) == 7
