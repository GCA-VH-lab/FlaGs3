

def test_until_stops_the_pipeline(tmp_path, monkeypatch):
	from flags3.cli import main
	from tests.synth import write_genome
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	listing = tmp_path / "list.txt"
	listing.write_text("WP_004\n")
	out = tmp_path / "run"
	assert main(["run", "-i", str(listing), "-o", str(out), "-nt", "-gd", str(genomes), "--offline",
		"-g", "1", "-c", "1", "--until", "extract"]) == 0
	assert (out / "extract" / "genes.tsv").is_file()
	assert not (out / "cluster").exists() and not (out / "report").exists()
	import pytest as _pytest
	with _pytest.raises(SystemExit):
		main(["run", "-i", str(listing), "-o", str(tmp_path / "bad"), "-nt", "-gd", str(genomes), "--offline", "--until", "nonsense"])
