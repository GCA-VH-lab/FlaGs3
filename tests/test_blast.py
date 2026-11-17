import stat

import pytest

from flags3 import blast
from flags3.run import RunDir
from flags3.schema import QueryTarget
from flags3.stage import Runner
from flags3.stages.blast import Blast, HitRow
from flags3.stages.fetch import Fetch
from flags3.tools import Tools
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


def _local_run(tmp_path, script_body, **config):
	tmp_path.mkdir(exist_ok=True)
	script = tmp_path / "blastp"
	script.write_text(script_body)
	script.chmod(script.stat().st_mode | stat.S_IEXEC)
	table = tmp_path / "tools.tsv"
	table.write_text("#name\tcommand\nblastp\t{} -query {{in}} -db {{db}} -evalue {{evalue}} -max_target_seqs {{hits}}\n".format(script))
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("MKVLLAKQ\tBLAST\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "tools": str(table), "blast_mode": "local", "blast_db": "mydb", "blast_inline": True,
		"blast_hits": 5, "cpu": 1, **config})
	cfg.save()
	return run, Runner(run, cfg, report=lambda m: None)


def test_blast_exclude_reaches_local_blastp(tmp_path):
	record = tmp_path / "argv.txt"
	body = "#!/bin/sh\necho \"$@\" > {}\nprintf 'WP_004\\t1e-50\\t200.5\\tprotein\\n'\n".format(record)
	run, runner = _local_run(tmp_path, body, blast_exclude="2759, 10239")
	assert runner.execute(Blast())
	assert record.read_text().strip().endswith("-negative_taxids 2759,10239")
	run, runner = _local_run(tmp_path / "plain", body.replace(str(record), str(tmp_path / "plain.txt")))
	assert runner.execute(Blast())
	assert "negative_taxids" not in (tmp_path / "plain.txt").read_text()


def test_blast_exclude_is_validated_and_local_only(tmp_path):
	assert blast.parse_taxids("") == [] and blast.parse_taxids(" 2759 ,2 ") == ["2759", "2"]
	for bad in ("eukaryota", "2759;rm", "-2759", "2759,"):
		with pytest.raises(blast.BlastError):
			blast.parse_taxids(bad)
	run, runner = _local_run(tmp_path, "#!/bin/sh\n", blast_exclude="2759", blast_mode="remote", user_email="a@b.c")
	with pytest.raises(blast.BlastError, match="local database"):
		runner.execute(Blast())
	run, runner = _local_run(tmp_path / "bad", "#!/bin/sh\n", blast_exclude="plants")
	with pytest.raises(blast.BlastError, match="taxonomy ids"):
		runner.execute(Blast())


def test_real_blastp_honours_the_exclusion(tmp_path):
	import shutil
	import subprocess
	if not (shutil.which("blastp") and shutil.which("makeblastdb")):
		pytest.skip("BLAST+ is not installed")
	seq_a = "MKVLLAKQRTVDEAGHILPWYSNCFMKVLLAKQRTVDEAGHILPWYSNCF"
	seq_b = "MKVLLAKQRTVDEAGHILPWYSNCFMKVLLAKQRTVDEAGHILPWYSNCA"
	(tmp_path / "db.faa").write_text(">WP_000000001.1 kept protein\n{}\n>XP_000000002.1 dropped protein\n{}\n".format(seq_a, seq_b))
	(tmp_path / "map.txt").write_text("WP_000000001.1 562\nXP_000000002.1 9606\n")
	subprocess.run(["makeblastdb", "-in", str(tmp_path / "db.faa"), "-dbtype", "prot", "-parse_seqids", "-taxid_map", str(tmp_path / "map.txt"),
		"-out", str(tmp_path / "db")], check=True, capture_output=True)
	tools = Tools.load(None)
	hits = blast.LocalBlast(tools["blastp"], str(tmp_path / "db"), 1e-5, 10, 1, tmp_path / "all").search(">q\n{}\n".format(seq_a))
	assert sorted(h.accession.split(".")[0] for h in blast.dedupe(hits)) == ["WP_000000001", "XP_000000002"]
	kept = blast.LocalBlast(tools["blastp"], str(tmp_path / "db"), 1e-5, 10, 1, tmp_path / "some", ["9606"]).search(">q\n{}\n".format(seq_a))
	assert [h.accession.split(".")[0] for h in blast.dedupe(kept)] == ["WP_000000001"]


def test_local_blast_asks_for_versioned_accessions():
	from flags3.tools import Tools
	argv = Tools.load().tools["blastp"].argv(db="d", evalue=1, hits=5, **{"in": "q"})
	assert argv[argv.index("-outfmt") + 1].split()[:3] == ["6", "qseqid", "saccver"]


BATCH_BLASTP = """#!/usr/bin/env python3
import sys
query = sys.argv[sys.argv.index("-query") + 1]
with open(query + ".calls", "a") as calls:
	calls.write("call\\n")
names = [line[1:].strip() for line in open(query) if line.startswith(">")]
lead = (lambda name: name + "\\t") if any("qseqid" in a for a in sys.argv) else (lambda name: "")
for i, name in enumerate(names):
	print("{}WP_00{}.1\\t1e-{}0\\t100.0\\thit".format(lead(name), i + 4, 5 - i))
print(lead(names[0]) + "WP_009.1\\t1e-9\\t50.0\\tshared")
print(lead(names[-1]) + "WP_009.1\\t1e-8\\t40.0\\tshared")
"""


def _batch_run(tmp_path, row):
	script = tmp_path / "blastp"
	script.write_text(BATCH_BLASTP)
	script.chmod(script.stat().st_mode | stat.S_IEXEC)
	table = tmp_path / "tools.tsv"
	table.write_text("#name\tcommand\nblastp\t{} -query {{in}} -db {{db}} {}\n".format(script, row))
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("WP_100.1\tBLAST\nWP_200.1\tBLAST\nWP_300.1\tBLAST\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "tools": str(table), "blast_mode": "local", "blast_db": "mydb", "blast_inline": True,
		"blast_hits": 5, "cpu": 1})
	cfg.save()
	return run, Runner(run, cfg, report=lambda m: None)


def test_several_local_queries_are_one_search(tmp_path, monkeypatch):
	monkeypatch.setattr(blast, "fetch_sequence", lambda accession: "MKVLLAKQ")
	run, runner = _batch_run(tmp_path, '-outfmt "6 qseqid saccver evalue bitscore stitle"')
	assert runner.execute(Blast())
	raw = run.stage_dir("blast") / "raw"
	assert (raw / "query.fasta.calls").read_text().count("call") == 1
	assert [l for l in (raw / "query.fasta").read_text().splitlines() if l.startswith(">")] == [">q1", ">q2", ">q3"]
	hits = HitRow.read(run.stage_file("blast", HitRow.FILE))
	queries = list(dict.fromkeys(h.query for h in hits))
	assert queries == ["WP_100.1", "WP_200.1", "WP_300.1"]
	assert [h.accession for h in hits if h.query == queries[0]] == ["WP_004.1", "WP_009.1"]
	assert [h.accession for h in hits if h.query == queries[1]] == ["WP_005.1"]
	assert [h.accession for h in hits if h.query == queries[2]] == ["WP_006.1", "WP_009.1"]
	lines = [l for l in run.stage_file("blast", "accessions.txt").read_text().splitlines() if not l.startswith("#")]
	assert lines == ["WP_004.1", "WP_009.1", "WP_005.1", "WP_006.1"]


def test_a_blastp_row_without_qseqid_searches_one_by_one(tmp_path, monkeypatch):
	monkeypatch.setattr(blast, "fetch_sequence", lambda accession: "MKVLLAKQ")
	run, runner = _batch_run(tmp_path, "")
	assert runner.execute(Blast())
	raw = run.stage_dir("blast") / "raw"
	assert (raw / "query.fasta.calls").read_text().count("call") == 3
	hits = HitRow.read(run.stage_file("blast", HitRow.FILE))
	assert [h.accession for h in hits if h.query == "WP_200.1"] == ["WP_004.1", "WP_009.1"]
	assert blast.LocalBlast(Tools.load().tools["blastp"], "d", 1, 5, 1, tmp_path).batched
