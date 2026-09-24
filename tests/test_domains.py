import gzip
from pathlib import Path

import pytest

from flags3 import domains
from flags3.run import RunDir
from flags3.schema import Annotation
from flags3.stage import Runner
from flags3.stages.domains import Domain, Domains

MOTIF = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQCGVWPAAVRESVPSLL"
NOISE = "MSSGGGGSSGGGSGSGGGSGGGSSGGGGSGGSSGGGGSSGGGSGSGGGSGGGSSGGGGSGGSSGGGGSSGGGSGSGGGSGGGSSGGGGSGG"


def _hmm(path: Path, name="Tiny__dom", accession="PF99999.1"):
	import pyhmmer
	from pyhmmer.easel import Alphabet, TextSequence
	alphabet = Alphabet.amino()
	seq = TextSequence(name=name.encode(), sequence=MOTIF).digitize(alphabet)
	builder = pyhmmer.plan7.Builder(alphabet)
	background = pyhmmer.plan7.Background(alphabet)
	hmm, _, _ = builder.build(seq, background)
	hmm.accession = accession.encode()
	with open(path, "wb") as out:
		hmm.write(out)
	return path


def test_source_parsing_and_files(tmp_path):
	hmm = _hmm(tmp_path / "Pfam-A.hmm")
	s = domains.HmmSource.parse(str(hmm))
	assert s.name == "Pfam-A" and s.files() == [hmm]
	s = domains.HmmSource.parse("df=" + str(tmp_path), {"df": (0.7, 0.5)})
	assert s.name == "df" and s.query_cov == 0.7 and s.files() == [hmm]
	assert s.group("Cas__Cas9") == "Cas" and s.group("PF00001") == "PF00001"
	with pytest.raises(domains.DomainError):
		domains.HmmSource.parse(str(tmp_path / "missing.hmm")).files()
	assert domains.parse_coverage(["0.5", "df=0.7,0.5"]) == {"": (0.5, 0.0), "df": (0.7, 0.5)}
	with pytest.raises(domains.DomainError):
		domains.parse_coverage(["1.5"])


def test_scanner_finds_motif_and_applies_coverage(tmp_path):
	hmm = _hmm(tmp_path / "tiny.hmm")
	sequences = {"hit": NOISE + MOTIF + NOISE, "miss": NOISE * 3, "half": NOISE + MOTIF[:120] + NOISE}
	hits = domains.DomainScanner([domains.HmmSource.parse(str(hmm))]).scan(sequences)
	by = {}
	for h in hits:
		by.setdefault(h.protein, []).append(h)
	assert set(by) == {"hit", "half"}
	h = by["hit"][0]
	assert h.name == "Tiny__dom" and h.accession == "PF99999.1" and h.group == "Tiny"
	assert len(NOISE) < h.start < h.end <= len(NOISE) + len(MOTIF) + 5
	strict = domains.DomainScanner([domains.HmmSource.parse("t=" + str(hmm), {"t": (0.0, 0.8)})]).scan(sequences)
	assert {h.protein for h in strict} == {"hit"}


def test_clans_and_interpro(tmp_path):
	clans = tmp_path / "clans.tsv.gz"
	with gzip.open(clans, "wt") as out:
		out.write("PF99999\tCL0001\tTinyClan\tTiny__dom\tdesc\nPF00002\t\t\tNoClan\tdesc\n")
	assert domains.load_clans(clans) == {"PF99999": "TinyClan", "Tiny__dom": "TinyClan"}
	table = tmp_path / "interpro.tsv"
	table.write_text("accession\tname\ttype\tpfam_members\tcharacterization_status\tinformativeness_label\tinterpretation\n"
		"IPR000001\tTiny domain\tDomain\tPF99999;PF00002\tchar\tinfo\tinterp\n")
	ip = domains.InterPro()
	assert ip.load(table) == 2
	assert ip.get("PF99999.1")["name"] == "Tiny domain"
	assert ip.get("PF12345") == {}
	assert domains.find_interpro(table) == table
	with pytest.raises(domains.DomainError):
		domains.find_interpro(tmp_path / "nope.tsv")


def test_domains_stage(tmp_path):
	hmm = _hmm(tmp_path / "tiny.hmm")
	clans = tmp_path / "clans.tsv"
	clans.write_text("PF99999\tCL0001\tTinyClan\tTiny__dom\tdesc\n")
	run = RunDir(tmp_path / "out").create("test", "")
	ext = run.reset("extract")
	(ext / "proteins.faa").write_text(">p1\n{}\n>p2\n{}\n".format(NOISE + MOTIF, NOISE * 2))
	status = run.status("extract")
	status.set("status", "ok")
	status.save()
	cfg = run.config()
	cfg.update({"domains": True, "hmmdb": ["tiny=" + str(hmm)], "clans": str(clans)})
	cfg.save()
	assert Runner(run, cfg, report=lambda m: None).execute(Domains())
	rows = Domain.read(run.stage_file("domains", Domain.FILE))
	assert [r.protein for r in rows] == ["p1"]
	assert rows[0].database == "tiny" and rows[0].clan == "TinyClan" and rows[0].pfam == "PF99999.1"
	assert rows[0].interpro == "-"
	ann = Annotation.read(run.stage_file("domains", Annotation.FILE))
	assert len(ann) == 1 and ann[0].kind == "wedge" and ann[0].space == "aa"
	assert ann[0].category == "TinyClan" and ann[0].label == "Tiny__dom" and ann[0].tool == "tiny"
	assert ann[0].score is not None and ann[0].score < 1e-10
	cfg.set("clans", None)
	cfg.save()
	Runner(run, cfg, report=lambda m: None).execute(Domains())
	ann = Annotation.read(run.stage_file("domains", Annotation.FILE))
	assert ann[0].category == "Tiny"


def test_missing_database_fails_softly(tmp_path):
	run = RunDir(tmp_path / "out").create("test", "")
	ext = run.reset("extract")
	(ext / "proteins.faa").write_text(">p1\nMKV\n")
	status = run.status("extract")
	status.set("status", "ok")
	status.save()
	cfg = run.config()
	cfg.update({"domains": True, "hmmdb": [str(tmp_path / "none.hmm")]})
	cfg.save()
	assert Runner(run, cfg, report=lambda m: None).execute(Domains()) is False
	assert "flags3 install pfam" in run.status("domains").get("error")
