import stat

import pytest

from flags3 import blast
from flags3.run import RunDir
from flags3.schema import QueryTarget
from flags3.stage import Runner
from flags3.stages.blast import Blast, HitRow
from flags3.stages.fetch import Fetch
from tests.synth import write_genome

FAKE_BLASTP = """#!/bin/sh
printf 'WP_004\\t1e-50\\t200.5\\tquery protein [Escherichia coli]\\n'
printf 'WP_008.1\\t1e-20\\t90.0\\tprotein eight\\n'
printf 'WP_008\\t1e-19\\t80.0\\tduplicate\\n'
printf 'garbage line\\n'
"""


def test_parse_query():
	q = blast.parse_query(["WP_000001.1"], "x")
	assert q.accession == "WP_000001.1" and q.sequence == ""
	q = blast.parse_query([">myprot desc", "MKV", "LLA"], "x")
	assert q.name == "myprot" and q.sequence == "MKVLLA" and q.accession is None
	assert blast.parse_query(["mkvlla"], "x").sequence == "MKVLLA"
	with pytest.raises(blast.BlastError):
		blast.parse_query([">h"], "x")
	with pytest.raises(blast.BlastError):
		blast.parse_query(["not a protein 123"], "x")
	assert blast.accession_in("ref|WP_000001.1|", "desc") == "WP_000001.1"
	assert blast.accession_in("gb|ABC12345.1|", "") == "ABC12345.1"
	assert blast.accession_in("xyz", "") is None


def test_dedupe_and_aliases():
	hits = [blast.BlastHit("WP_1.1", 1, 1, ""), blast.BlastHit("WP_1.2", 1, 1, ""), blast.BlastHit("WP_2", 1, 1, "")]
	assert [h.accession for h in blast.dedupe(hits)] == ["WP_1.1", "WP_2"]
	assert blast.RemoteBlast("refseq_select", 1e-5, 10, "", 1).database == "refseq_select_prot"


def test_local_blast_stage_feeds_fetch(tmp_path):
	script = tmp_path / "blastp"
	script.write_text(FAKE_BLASTP)
	script.chmod(script.stat().st_mode | stat.S_IEXEC)
	table = tmp_path / "tools.tsv"
	table.write_text("#name\tcommand\nblastp\t{} -query {{in}} -db {{db}} -evalue {{evalue}} -max_target_seqs {{hits}}\n".format(script))
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("WP_009\nMKVLLAKQ\tBLAST\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True, "tools": str(table),
		"blast_mode": "local", "blast_db": "mydb", "blast_inline": True, "blast_hits": 5})
	cfg.save()
	runner = Runner(run, cfg, report=lambda m: None)
	stage = Blast()
	assert stage.wanted(cfg)
	assert runner.execute(stage)
	hits = HitRow.read(run.stage_file("blast", HitRow.FILE))
	assert [h.accession for h in hits] == ["WP_004", "WP_008.1"]
	assert hits[0].query == "query" and hits[0].evalue == 1e-50
	lines = [l for l in run.stage_file("blast", "accessions.txt").read_text().splitlines() if not l.startswith("#")]
	assert lines == ["WP_004", "WP_008.1"]
	assert runner.execute(Fetch())
	targets = QueryTarget.read(run.stage_file("fetch", QueryTarget.FILE))
	assert [t.query for t in targets] == ["WP_009", "WP_004", "WP_008.1"]
	assert targets[1].status == "ok"


def test_blast_stage_without_query_fails(tmp_path):
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("WP_009\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "blast_inline": True, "user_email": "x@y"})
	cfg.save()
	with pytest.raises(blast.BlastError):
		Runner(run, cfg, report=lambda m: None).execute(Blast())
	cfg.set("blast_inline", False)
	assert not Blast().wanted(cfg)
