from flags3 import fasta
from flags3.run import RunDir
from flags3.schema import MISSING
from flags3.stage import Runner
from flags3.stages.cluster import Cluster
from flags3.stages.extract import Extract
from flags3.stages.fetch import Fetch
from flags3.stages.report import NeighbourhoodRow, QueryStatus, Report
from tests.synth import write_genome


def _run(tmp_path):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "myrun").create("test", "flags3 run -i list.txt")
	(run.input_dir / "list.txt").write_text("WP_004\nWP_999\nWP_009\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "gene": 2,
		"cluster_method": "jackhmmer", "cpu": 1})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	for stage in (Fetch(), Extract(), Cluster(), Report()):
		assert runner.execute(stage)
	return run


def test_report_tables(tmp_path):
	run = _run(tmp_path)
	rows = NeighbourhoodRow.read(run.stage_file("report", NeighbourhoodRow.FILE))
	assert [r.accession for r in rows if r.row_id == "WP_004|GCF_TEST"] == ["WP_006", "WP_005", "WP_004", "pseudogene*", "rna_001"]
	assert rows[0].row_id == "WP_004|GCF_TEST" and rows[0].strand == "+" and rows[0].offset == -2 and rows[0].domains == MISSING
	status = {q.query: q for q in QueryStatus.read(run.stage_file("report", QueryStatus.FILE))}
	assert list(status) == ["WP_004", "WP_999", "WP_009"]
	assert status["WP_004"].rows == 1 and status["WP_004"].status == "ok"
	assert status["WP_999"].rows == 0 and status["WP_999"].status.startswith("unresolved")
	summary = run.stage_file("report", "run_summary.txt").read_text()
	assert "3 queries, 2 rows" in summary and "cluster      ok" in summary and "flags3 run -i list.txt" in summary
	assert "elapsed: " in summary and "stages total" in summary
	families = run.stage_file("report", "families.tsv").read_text().splitlines()
	assert families[0].startswith("family\tlabel\tsize")


def test_report_fastas_and_legacy(tmp_path):
	run = _run(tmp_path)
	flanking = list(fasta.read(run.stage_file("report", "flanking.faa")))
	assert ("WP_005", "protein five") in [(n, d) for n, d, _ in flanking]
	assert "WP_004" not in {n for n, _, _ in flanking}
	queries = list(fasta.read(run.stage_file("report", "queries.faa")))
	assert ("WP_004|GCF_TEST", "Escherichia coli K-12") in [(n, d) for n, d, _ in queries]
	legacy = run.stage_file("report", "legacy")
	names = sorted(p.name for p in legacy.iterdir())
	assert "myrun_operon.tsv" in names and "myrun_outdesc.txt" in names and "myrun_all.fasta" in names
	operon = (legacy / "myrun_operon.tsv").read_text().splitlines()
	assert operon[0].startswith("#query\tassembly\tspecies\tfamily\tstrand")
	query_line = next(l for l in operon if "\tWP_004\t" in l).split("\t")
	assert query_line[4] == "+" and query_line[5] == "0"
	first = next(l for l in operon if "\t-2\t" in l).split("\t")
	assert first[11] == "WP_006" and first[4] == "-"
	status = (legacy / "myrun_QueryStatus.txt").read_text().splitlines()
	assert status[1] == "WP_004\tGCF_TEST\tYes" and status[2] == "WP_999\t-\tNo"
	all_fasta = (legacy / "myrun_all.fasta").read_text().splitlines()
	assert all_fasta[0] == ">WP_004|GCF_TEST" and len(all_fasta[1]) > 60
	assert ">WP_005|protein five" in all_fasta
