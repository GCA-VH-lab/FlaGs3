import pytest

from flags3.schema import (Annotation, Gene, SchemaError, Window, bp_subject,
	row_id, split_bp_subject, split_row_id)


def test_gene_round_trip(tmp_path):
	genes = [
		Gene("WP_1|GCF_1", "GCF_1", "NZ_1", "WP_1", 100, 400, "+", "hypothetical protein", False, 0),
		Gene("WP_1|GCF_1", "GCF_1", "NZ_1", "rna_1", 500, 580, "-", "tRNA-Ala", True, 1),
	]
	path = tmp_path / Gene.FILE
	assert Gene.write(path, genes) == 2
	assert Gene.read(path) == genes
	assert path.read_text().splitlines()[0] == "\t".join(Gene.columns())


def test_gene_validation():
	with pytest.raises(SchemaError):
		Gene("r", "a", "c", "x", 10, 5, "+", "p", False, 0)
	with pytest.raises(SchemaError):
		Gene("r", "a", "c", "x", 1, 5, "?", "p", False, 0)


def test_window_nesting():
	w = Window("q|a", "a", "c", 5000, 100, 900, 50, 950, 1, 1950)
	assert w.subject == "a|c"
	with pytest.raises(SchemaError):
		Window("q|a", "a", "c", 5000, 100, 900, 200, 6000, 1, 5000)


def test_annotation_kinds_and_spaces(tmp_path):
	rows = [
		Annotation("WP_1", "aa", 12, 140, "wedge", "PF00001", "7tm_1", "hmmer", 1e-30),
		Annotation("GCF_1|NZ_1", "bp", 1200, 8900, "band", "CRISPR-Cas", "CAS-I-E", "padloc", 0.9),
		Annotation("WP_2", "-", None, None, "fill", "family:3", "3", "jackhmmer", None),
	]
	path = tmp_path / Annotation.FILE
	Annotation.write(path, rows)
	back = Annotation.read(path)
	assert back[1:] == rows[1:] and back[0].score == 1e-30
	assert back[2].whole and not back[0].whole
	lines = path.read_text().splitlines()
	assert lines[3].split("\t")[1:4] == ["-", "-", "-"]


def test_annotation_validation():
	with pytest.raises(SchemaError):
		Annotation("x", "aa", 1, 2, "sparkle", "c", "l", "t", None)
	with pytest.raises(SchemaError):
		Annotation("x", "nm", 1, 2, "band", "c", "l", "t", None)
	with pytest.raises(SchemaError):
		Annotation("x", "-", 1, 2, "fill", "c", "l", "t", None)
	with pytest.raises(SchemaError):
		Annotation("x", "aa", None, None, "wedge", "c", "l", "t", None)


def test_missing_column_is_reported(tmp_path):
	path = tmp_path / "bad.tsv"
	path.write_text("subject\tspace\n")
	with pytest.raises(SchemaError):
		Annotation.read(path)


def test_identifiers():
	assert split_row_id(row_id("WP_1", "GCF_1")) == ("WP_1", "GCF_1")
	assert split_bp_subject(bp_subject("GCF_1", "NZ_1")) == ("GCF_1", "NZ_1")
	with pytest.raises(SchemaError):
		split_row_id("WP_1")
