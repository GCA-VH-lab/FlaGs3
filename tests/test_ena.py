import gzip
from pathlib import Path

import pytest

from flags3 import ena, fasta, ncbi, net
from flags3.genome import PSEUDOGENE, GeneTable, GenomeFasta, ProteinFasta, RnaFasta
from flags3.run import RunDir
from flags3.schema import MISSING, Failure, GenomeFiles, QueryTarget
from flags3.stage import Runner
from flags3.stages.fetch import Fetch

ORGANISM = "Testus syntheticus"
PROTEIN = "MKKLLPTAAAGLLLLAAQPAMAKKKLLPTAAAGLLLLAAQPAMAKKKLLPTAAAGLLL"


def _feature(kind, location, **qualifiers):
	lines = ["FT   {:<16}{}".format(kind, location)]
	for key, value in qualifiers.items():
		lines.append("FT                   /{}{}".format(key, "" if value is None else '="{}"'.format(value)))
	return "\n".join(lines)


def _record(accession, topology, length, features, description="synthetic record"):
	seq = ("acgt" * (length // 4 + 1))[:length]
	body = ["ID   {}; SV 1; {}; genomic DNA; STD; PRO; {} BP.".format(accession, topology, length), "XX",
		"AC   {};".format(accession), "XX", "DE   {} {}".format(ORGANISM, description), "XX",
		"OS   {}".format(ORGANISM), "OC   Bacteria.", "XX", "FH   Key             Location/Qualifiers", "FH",
		_feature("source", "1..{}".format(length), organism=ORGANISM)]
	body.extend(features)
	body.extend(["XX", "SQ   Sequence {} BP; {} A; {} C; {} G; {} T; 0 other;".format(length, *[seq.count(b) for b in "acgt"])])
	for i in range(0, length, 60):
		chunk = seq[i:i + 60]
		body.append("     " + " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10)).ljust(66) + str(min(i + 60, length)))
	body.append("//")
	return "\n".join(body) + "\n"


CHROMOSOME = _record("CP000001", "circular", 1200, [
	_feature("gene", "1..180", locus_tag="ENA_0001"),
	_feature("CDS", "1..180", locus_tag="ENA_0001", product="protein one", protein_id="AAA00001.1", translation=PROTEIN),
	_feature("gene", "complement(200..380)", locus_tag="ENA_0002"),
	_feature("CDS", "complement(200..380)", locus_tag="ENA_0002", product="broken protein", pseudo=None),
	_feature("gene", "400..480", locus_tag="ENA_0003", pseudo=None),
	_feature("gene", "500..575", locus_tag="ENA_0004"),
	_feature("tRNA", "500..575", locus_tag="ENA_0004", product="tRNA-Ala"),
	_feature("gene", "600..780", locus_tag="ENA_0005"),
	_feature("CDS", "600..780", locus_tag="ENA_0005", product="no id; odd=product", translation=PROTEIN),
	_feature("gene", "complement(800..980)", locus_tag="ENA_0006"),
	_feature("CDS", "complement(800..980)", locus_tag="ENA_0006", product="query protein", protein_id="AAA00006.1", translation=PROTEIN),
], "chromosome, complete genome.")
PLASMID = _record("CP000002", "linear", 400, [
	_feature("CDS", "10..190", locus_tag="ENA_1001", product="plasmid protein", protein_id="AAA01001.1", translation=PROTEIN),
], "plasmid pTest.")


def _convert(tmp_path, text=CHROMOSOME + PLASMID, name="CP000001"):
	flat = tmp_path / "flat.embl"
	flat.write_text(text)
	return ena.EmblRecords(flat).write(tmp_path, name)


def test_accession_patterns():
	assert ena.is_ena("CP003200") and ena.is_ena("LT906474.1") and ena.is_ena("A00145")
	assert ena.is_ena("CABHOY01") and ena.is_ena("CABHOY010000000") and ena.is_ena("AAAK03000000")
	assert ena.is_ena("CABHOY010000012")
	for other in ("GCA_000005845.2", "GCF_000005845.2", "MGYG000454827", "AAC73113.1", "WP_005328829.1", "ENA_0001"):
		assert not ena.is_ena(other), other
	assert ena.cache_name("CABHOY010000000") == "CABHOY01" and ena.cache_name("CP003200.2") == "CP003200.2"
	assert ena.url_of("CABHOY01") == "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/cab/CABHOY01.dat.gz"
	assert ena.url_of("CP003200") == "https://www.ebi.ac.uk/ena/browser/api/embl/CP003200"


def test_embl_becomes_the_four_slots(tmp_path):
	got = _convert(tmp_path)
	assert sorted(got) == ["faa", "genome", "gff", "rna"]
	assert all(net.intact(p) for p in got.values())
	table = GeneTable(got["gff"])
	assert table.contigs["CP000001.1"].circular and table.contigs["CP000001.1"].length == 1200
	assert not table.contigs["CP000002.1"].circular
	first = [g for g in table.genes if g.contig == "CP000001.1"]
	assert [g.accession for g in first] == ["AAA00001.1", PSEUDOGENE, PSEUDOGENE, "ENA_0004", "ENA_0005", "AAA00006.1"]
	assert [g.strand for g in first] == ["+", "-", "+", "+", "+", "-"]
	assert first[3].is_rna and first[3].biotype == "tRNA"
	assert first[4].product == "no id; odd=product"
	assert [g.accession for g in table.genes if g.contig == "CP000002.1"] == ["AAA01001.1"]
	proteins = ProteinFasta(got["faa"])
	assert set(proteins.sequences) == {"AAA00001.1", "ENA_0005", "AAA00006.1", "AAA01001.1"}
	assert proteins.organism_of("AAA00006.1") == ORGANISM and proteins.get("AAA00001.1") == PROTEIN
	rna = RnaFasta(got["rna"])
	assert len(rna.get(first[3])) == 76
	genome = GenomeFasta(got["genome"])
	assert len(genome.contigs["CP000001.1"]) == 1200 and len(genome.slice(first[5])) == 181


def test_record_without_rna_or_cds(tmp_path):
	bare = _record("CP000003", "linear", 300, [_feature("CDS", "1..90", locus_tag="X_1", product="p", translation="MKK")])
	got = _convert(tmp_path, bare, "CP000003")
	assert sorted(got) == ["faa", "genome", "gff", "rna"]
	assert net.intact(got["rna"]) and list(fasta.read(got["rna"])) == [] and not list(tmp_path.glob("*.part"))
	empty = _record("CP000004", "linear", 300, [])
	got = _convert(tmp_path, empty, "CP000004")
	assert sorted(got) == ["genome", "gff", "rna"] and GeneTable(got["gff"]).genes == []


class FakeEna(ena.EnaGenomes):
	calls = []

	def stream(self, url, local):
		FakeEna.calls.append(url)
		if "CP000001" in url:
			local.write_text(CHROMOSOME + PLASMID)
			return True
		if "CABHOY01" in url:
			with gzip.open(local, "wt") as out:
				out.write(_record("CABHOY010000001", "linear", 400, [
					_feature("CDS", "10..190", locus_tag="WGS_0001", product="contig protein", protein_id="BBB00001.1", translation=PROTEIN)]))
			return True
		self.failures[url] = "HTTP 404"
		return False


class NoNcbi(ncbi.NcbiGenomes):
	def listing(self, url):
		raise AssertionError("NCBI must not be contacted for ENA genomes")


def _setup(tmp_path, monkeypatch, lines, **config):
	monkeypatch.setattr(Fetch, "RETRY_PAUSE", 0.0)
	monkeypatch.setattr(ena, "EnaGenomes", FakeEna)
	monkeypatch.setattr(ncbi, "NcbiGenomes", NoNcbi)
	FakeEna.calls = []
	run = RunDir(tmp_path / "out").create("test", "")
	(run.input_dir / "list.txt").write_text("\n".join(lines) + "\n")
	cfg = run.config()
	cfg.update({"inputs": "list.txt", "genomes": str(tmp_path / "cache"), **config})
	cfg.save()
	return run, Runner(run, cfg, report=lambda m: None)


def _read(run):
	genomes = {g.assembly: g for g in GenomeFiles.read(run.stage_file("fetch", GenomeFiles.FILE))}
	targets = QueryTarget.read(run.stage_file("fetch", QueryTarget.FILE))
	failures = Failure.read(run.stage_file("fetch", Failure.FILE))
	return genomes, targets, failures


def test_paired_ena_query_is_fetched_without_ncbi(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["AAA00006.1\tCP000001", "BBB00001.1\tCABHOY010000000"], **{"fetch.slots": "rna,genome"})
	assert runner.execute(Fetch())
	genomes, targets, failures = _read(run)
	assert [(t.assembly, t.status, t.via) for t in targets] == [("CP000001", "ok", "paired"), ("CABHOY010000000", "ok", "paired")]
	assert genomes["CP000001"].source == "ena" and genomes["CP000001"].usable
	assert Path(genomes["CP000001"].gff).name == "CP000001.gff.gz" and genomes["CP000001"].rna != MISSING
	assert genomes["CABHOY010000000"].usable and Path(genomes["CABHOY010000000"].genome).name == "CABHOY01.fna.gz"
	assert failures == [] and len(FakeEna.calls) == 2
	assert not list((tmp_path / "cache").glob("*.embl")) and not list((tmp_path / "cache").glob("*.dat.gz"))
	runner.execute(Fetch())
	genomes, targets, _ = _read(run)
	assert genomes["CP000001"].source == "cache" and targets[1].assembly == "CABHOY01" and len(FakeEna.calls) == 2


def test_bare_ena_accession_is_rejected(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["CP000001", "AAA00001.1\tCP000001"])
	assert runner.execute(Fetch())
	_, targets, _ = _read(run)
	assert targets[0].assembly == MISSING and "PROTEIN<tab>CP000001" in targets[0].status and targets[0].via == "ena"
	assert targets[1].status == "ok"
	run, runner = _setup(tmp_path / "off", monkeypatch, ["CP000001"], offline=True)
	with pytest.raises(RuntimeError):
		runner.execute(Fetch())
	_, targets, _ = _read(run)
	assert targets[0].status.startswith("unresolved: an ENA sequence accession")


def test_ena_failure_is_reported(tmp_path, monkeypatch):
	run, runner = _setup(tmp_path, monkeypatch, ["AAA00001.1\tCP000009"])
	with pytest.raises(RuntimeError):
		runner.execute(Fetch())
	_, targets, failures = _read(run)
	assert targets[0].status == "unresolved: genome files unavailable"
	assert any(f.subject == "CP000009" and "404" in f.reason for f in failures)


def test_extract_reads_the_converted_genome(tmp_path, monkeypatch):
	from flags3.schema import Gene, RowInfo
	from flags3.stages.extract import Extract, PROTEINS, QUERIES, RNA
	run, runner = _setup(tmp_path, monkeypatch, ["AAA00006.1\tCP000001"], gene=3, **{"fetch.slots": "rna,genome"})
	assert runner.execute(Fetch()) and runner.execute(Extract())
	genes = Gene.read(run.stage_file("extract", Gene.FILE))
	assert [g.accession for g in genes] == [PSEUDOGENE, "ENA_0004", "ENA_0005", "AAA00006.1"]
	rows = RowInfo.read(run.stage_file("extract", RowInfo.FILE))
	assert [(r.accession, r.strand, r.species) for r in rows] == [("AAA00006.1", "-", ORGANISM)]
	proteins = {name for name, _, _ in fasta.read(run.stage_file("extract", PROTEINS))}
	assert "AAA00006.1|CP000001" in proteins or "AAA00006.1" in proteins
	assert run.stage_file("extract", RNA).exists()
