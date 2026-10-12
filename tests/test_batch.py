from pathlib import Path

from flags3.batch import Batch
from flags3.run import KeyValueFile
from tests.synth import write_genome


def test_batch_runs_each_list_and_skips_finished(tmp_path, monkeypatch):
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	lists = tmp_path / "lists"
	lists.mkdir()
	(lists / "alpha.txt").write_text("WP_004\n")
	(lists / "beta.txt").write_text("WP_009\n")
	out = tmp_path / "runs"
	options = ["-gd", str(genomes), "--offline", "-g", "1", "-nf", "-c", "1"]
	assert Batch(lists, out, options, retry=False).run() == 0
	assert (out / "alpha" / "report" / "status.tsv").is_file() and (out / "beta" / "report" / "status.tsv").is_file()
	rows = (out / "batch_status.tsv").read_text().splitlines()
	assert rows[0] == "list\tstatus\tseconds\tfinished" and len(rows) == 3 and all("\tok\t" in r for r in rows[1:])
	stamp = (out / "alpha" / "report" / "status.tsv").stat().st_mtime
	assert Batch(lists, out, options, retry=False).run() == 0
	assert (out / "alpha" / "report" / "status.tsv").stat().st_mtime == stamp
	(out / "beta" / "report" / "status.tsv").unlink()
	assert Batch(lists, out, options, retry=False).run() == 1
	assert Batch(lists, out, options, retry=True).run() == 0


def test_batch_through_the_command_line(tmp_path, monkeypatch, capsys):
	from flags3.cli import main
	genomes = tmp_path / "genomes"
	write_genome(genomes)
	lists = tmp_path / "lists"
	lists.mkdir()
	(lists / "one.txt").write_text("WP_004\n")
	out = tmp_path / "runs"
	code = main(["batch", str(lists), "-o", str(out), "--", "-gd", str(genomes), "--offline", "-g", "1", "-nf", "-c", "1"])
	assert code == 0 and (out / "one" / "report" / "status.tsv").is_file()
	code = main(["batch", "-o", str(out), str(lists), "--", "-gd", str(genomes), "--offline", "-g", "1", "-nf", "-c", "1"])
	assert code == 0 and "skipped" in capsys.readouterr().out
