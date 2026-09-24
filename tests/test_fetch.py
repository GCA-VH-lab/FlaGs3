import gzip
import shutil
from pathlib import Path

import pytest

from flags3 import mgnify, ncbi, net
from flags3.run import RunDir
from flags3.schema import MISSING, Failure, GenomeFiles, QueryTarget
from flags3.stage import Runner
from flags3.stages.fetch import Fetch, GenomeCache
from tests.synth import write_genome

IPG_TEXT = """Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tProtein Name\tOrganism\tStrain\tAssembly
1\tRefSeq\tNZ_CP001\t100\t400\t+\tWP_001\tprotein one\tEscherichia coli\tK-12\tGCF_000000001.1
1\tINSDC\tCP001\t100\t400\t+\tAAA001.1\tprotein one\tEscherichia coli\tK-12\tGCA_000000001.1
1\tRefSeq\tNZ_CP002\t100\t400\t+\tWP_001\tprotein one\tEscherichia coli\tB\tGCF_000000002.1
2\tRefSeq\tNZ_CP003\t1\t9\t+\tWP_777\tother\tSomething\t\tGCF_000000003.1
"""


class FakeMapper(ncbi.IpgMapper):
	def __init__(self, *a, **k):
		super().__init__(*a, **k)
		self.calls = []

	def _fetch_ipg(self, proteins):
		self.calls.append(list(proteins))
		return IPG_TEXT


class FakeNcbi(ncbi.NcbiGenomes):
	source = None

	def listing(self, url):
		assembly = url.rstrip("/").split("/")[-1]
		return '<a href="{0}_ASM1v1/">{0}_ASM1v1/</a>'.format("GCF_000000001.1") if url.endswith("/000/001") else ""

	def stream(self, url, local):
		if not url.startswith("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/001/GCF_000000001.1_ASM1v1/GCF_000000001.1_ASM1v1"):
			self.failures[url] = "HTTP 404"
			return False
		slot = next(s for s, suffix in ncbi.SUFFIX.items() if url.endswith(suffix))
		plain = {"gff": "GCF_TEST_genomic.gff", "faa": "GCF_TEST_protein.faa",
			"rna": "GCF_TEST_rna_from_genomic.fna", "genome": "GCF_TEST_genomic.fna"}[slot]
		with open(self.source / plain, "rb") as src, gzip.open(local, "wb") as dst:
			shutil.copyfileobj(src, dst)
		return True


def _setup(tmp_path, monkeypatch, queries, **config):
	source = tmp_path / "source"
	write_genome(source)
	FakeNcbi.source = source
	monkeypatch.setattr(ncbi, "IpgMapper", FakeMapper)
	monkeypatch.setattr(ncbi, "NcbiGenomes", FakeNcbi)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("\n".join(queries) + "\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(tmp_path / "cache"), "user_email": "x@y.z", **config})
	cfg.save()
	messages = []
	runner = Runner(run, cfg, report=messages.append)
	return run, runner


def _read(run):
	genomes = {g.assembly: g for g in GenomeFiles.read(run.stage_file("fetch", GenomeFiles.FILE))}
	targets = QueryTarget.read(run.stage_file("fetch", QueryTarget.FILE))
	failures = Failure.read(run.stage_file("fetch", Failure.FILE))
	return genomes, targets, failures


def test_ipg_resolution_and_download(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_001", "WP_999"])
	assert runner.execute(Fetch())
	genomes, targets, failures = _read(run)
	assert [t.assembly for t in targets] == ["GCF_000000001.1", MISSING]
	assert targets[0].accessions == "WP_001"
	assert targets[1].status.startswith("unresolved")
	g = genomes["GCF_000000001.1"]
	assert g.source == "ncbi" and g.usable
	assert g.rna == MISSING and g.genome == MISSING
	assert Path(g.gff).name == "GCF_000000001.1_genomic.gff.gz"
	assert net.intact(Path(g.faa))


def test_slots_and_cache_reuse(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_001"])
	run.config().update({"fetch.slots": "rna,genome"})
	cfg = run.config()
	cfg.set("fetch.slots", "rna,genome")
	cfg.save()
	runner.config = cfg
	runner.execute(Fetch())
	genomes, _, _ = _read(run)
	assert genomes["GCF_000000001.1"].rna != MISSING and genomes["GCF_000000001.1"].genome != MISSING
	calls = []
	monkeypatch.setattr(FakeNcbi, "stream", lambda self, url, local: calls.append(url) or False)
	runner.execute(Fetch())
	genomes, _, _ = _read(run)
	assert calls == []
	assert genomes["GCF_000000001.1"].source == "cache" and genomes["GCF_000000001.1"].usable


def test_broken_cache_file_is_redownloaded(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_001"])
	runner.execute(Fetch())
	genomes, _, _ = _read(run)
	faa = Path(genomes["GCF_000000001.1"].faa)
	faa.write_bytes(faa.read_bytes()[:40])
	assert not net.intact(faa)
	runner.execute(Fetch())
	genomes, _, _ = _read(run)
	assert net.intact(Path(genomes["GCF_000000001.1"].faa))
	assert genomes["GCF_000000001.1"].source == "ncbi"


def test_paired_query_keeps_assembly_and_remap_fills_aliases(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_001\tGCF_000000001.1", "WP_001\tGCF_000000009.1"])
	runner.execute(Fetch())
	_, targets, _ = _read(run)
	assert targets[0].accessions == MISSING
	assert targets[1].status.startswith("unresolved")
	cfg = run.config()
	cfg.set("remap", True)
	cfg.save()
	runner.config = cfg
	runner.execute(Fetch())
	_, targets, failures = _read(run)
	assert targets[0].accessions == "WP_001"
	assert targets[1].assembly == "GCF_000000001.1"
	assert any("remapped from GCF_000000009.1" in f.reason for f in failures)


def test_cross_db_and_max_assemblies(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_001"], max_assemblies=3)
	runner.execute(Fetch())
	_, targets, _ = _read(run)
	assert [t.assembly for t in targets] == ["GCF_000000001.1", "GCF_000000002.1", "GCA_000000001.1"]
	assert targets[2].accessions == "AAA001.1"
	cfg = run.config()
	cfg.set("no_cross_db", True)
	cfg.save()
	runner.config = cfg
	runner.execute(Fetch())
	_, targets, failures = _read(run)
	assert [t.assembly for t in targets] == ["GCF_000000001.1", "GCF_000000002.1"]
	assert any("cross-database" in f.reason for f in failures)


def test_email_required_unless_offline(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_001"], user_email="")
	with pytest.raises(RuntimeError):
		runner.execute(Fetch())


def test_offline_prefix_match_and_protein_scan(tmp_path):
	genomes = tmp_path / "cache"
	write_genome(genomes)
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("WP_005\tGCF_TES\nWP_008\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(genomes), "offline": True})
	cfg.save()
	Runner(run, cfg, report=lambda m: None).execute(Fetch())
	genomes_, targets, _ = _read(run)
	assert [t.assembly for t in targets] == ["GCF_TEST", "GCF_TEST"]
	assert genomes_["GCF_TEST"].source == "local"


def test_cache_ignores_partial_files(tmp_path):
	write_genome(tmp_path)
	(tmp_path / "GCF_OTHER_genomic.gff.gz.part").write_bytes(b"x")
	cache = GenomeCache(tmp_path)
	assert set(cache.index) == {"GCF_TEST"}
	assert cache.find("GCF_TEST.1") == "GCF_TEST"


def test_mgnify_url_map():
	payload = {"downloads": [{"url": "https://x/y/MGYG1.gff.gz?token=1"}, {"links": {"download": "https://x/MGYG1.faa"}}, "https://x/MGYG1.fna"]}
	assert mgnify.url_map(payload) == {"MGYG1.gff.gz": "https://x/y/MGYG1.gff.gz?token=1", "MGYG1.faa": "https://x/MGYG1.faa", "MGYG1.fna": "https://x/MGYG1.fna"}
	assert mgnify.is_mgnify("MGYG000454827") and not mgnify.is_mgnify("GCF_1")


def test_local_first_mixes_own_genome_with_ipg(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["WP_008", "WP_777"], local_first=True)
	write_genome(tmp_path / "cache")
	assert runner.execute(Fetch())
	_, targets, _ = _read(run)
	assert [(t.query, t.assembly) for t in targets] == [("WP_008", "GCF_TEST"), ("WP_777", "GCF_000000003.1")]
	assert targets[0].status == "ok"
