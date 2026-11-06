import gzip
import os
import time

import pytest

from flags3 import cli
from flags3.cache import CacheReport
from flags3.stages.fetch import GenomeCache

DAY = 86400


def _genome(directory, base, age_days, slots=("gff", "faa", "fna", "rna.fna")):
	stamp = time.time() - age_days * DAY
	for ext in slots:
		path = directory / "{}.{}.gz".format(base, ext)
		with gzip.open(path, "wt") as out:
			out.write(">x\nACGT\n" if "fna" in ext else "##gff-version 3\n")
		os.utime(path, (stamp, stamp))


def test_cache_hit_touches_atime_but_not_mtime(tmp_path):
	_genome(tmp_path, "GCF_000000001.1", 100)
	cache = GenomeCache(tmp_path)
	before = {p: p.stat() for p in cache.index["GCF_000000001.1"].values()}
	got = cache.files("GCF_000000001.1", ["gff", "faa"])
	assert set(got) == {"gff", "faa"}
	for slot, path in cache.index["GCF_000000001.1"].items():
		stat = path.stat()
		assert stat.st_mtime == before[path].st_mtime
		if slot in got:
			assert time.time() - stat.st_atime < 60
		else:
			assert stat.st_atime == before[path].st_atime


def test_report_and_prune(tmp_path):
	_genome(tmp_path, "GCF_000000001.1", 100)
	_genome(tmp_path, "GCA_000000002.1", 40, ("gff", "faa"))
	_genome(tmp_path, "MGYG000000003", 2)
	(tmp_path / "junk.part").write_text("")
	report = CacheReport(tmp_path)
	assert [e.base for e in report.entries] == ["GCA_000000002.1", "GCF_000000001.1", "MGYG000000003"]
	assert "3 genomes" in report.summary() and "between 2 and 100 days ago" in report.summary()
	assert [e.base for e in report.older_than(30)] == ["GCA_000000002.1", "GCF_000000001.1"]
	stale = report.prune(30, dry_run=True)
	assert len(stale) == 2 and all(p.exists() for e in stale for p in e.files.values())
	stale = report.prune(30)
	assert len(stale) == 2 and not any(p.exists() for e in stale for p in e.files.values())
	remaining = sorted(p.name for p in tmp_path.iterdir())
	assert remaining == ["MGYG000000003.faa.gz", "MGYG000000003.fna.gz", "MGYG000000003.gff.gz", "MGYG000000003.rna.fna.gz", "junk.part"]
	assert CacheReport(tmp_path).prune(30) == []


def test_recent_use_protects_an_old_download(tmp_path):
	_genome(tmp_path, "GCF_000000001.1", 100)
	GenomeCache(tmp_path).files("GCF_000000001.1", ["gff", "faa", "genome", "rna"])
	assert CacheReport(tmp_path).older_than(30) == []


def test_cli(tmp_path, capsys):
	_genome(tmp_path, "GCF_000000001.1", 100)
	_genome(tmp_path, "GCF_000000002.1", 1)
	assert cli.main(["cache", str(tmp_path)]) == 0
	assert "2 genomes" in capsys.readouterr().out
	assert cli.main(["cache", str(tmp_path), "--prune", "30", "--dry-run"]) == 0
	out = capsys.readouterr().out
	assert "would delete GCF_000000001.1" in out and "would delete 1 genome," in out and (tmp_path / "GCF_000000001.1.gff.gz").exists()
	assert cli.main(["cache", str(tmp_path), "--prune", "30"]) == 0
	assert "deleted GCF_000000001.1" in capsys.readouterr().out and not (tmp_path / "GCF_000000001.1.gff.gz").exists()
	with pytest.raises(SystemExit):
		cli.main(["cache", str(tmp_path / "missing")])
	assert not (tmp_path / "missing").exists()
