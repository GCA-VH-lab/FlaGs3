import pytest

from flags3.tools import Tool, ToolError, Tools


def test_shipped_table_loads():
	tools = Tools.load()
	assert "jackhmmer" in tools
	assert tools["jackhmmer"].engine == "jackhmmer"
	assert tools["jackhmmer"].options["iterations"] == "3"
	assert tools["mmseqs_cluster"].engine == "mmseqs"
	assert tools.names("clustering") == ["jackhmmer", "mmseqs_cluster", "mmseqs_cluster_exhaustive", "nhmmer"]


def test_user_table_overrides_but_cannot_add(tmp_path):
	user = tmp_path / "tools_table.tsv"
	user.write_text("#name\tcommand\tdirectory\tscan_range\tengine\toptions\n"
		"sismis\t/opt/sismis/bin/sismis run -g {in} -o {out}\t\t50000\t\t\n")
	tools = Tools.load(user)
	assert tools["sismis"].command.startswith("/opt/sismis")
	assert tools["sismis"].scan_range == 50000
	assert tools["mafft"].command.startswith("mafft")
	user.write_text("#name\tcommand\nnewtool\tx\n")
	with pytest.raises(ToolError):
		Tools.load(user)
	with pytest.raises(ToolError):
		Tools.load(tmp_path / "missing.tsv")


def test_argv_fills_and_splits():
	t = Tool("x")
	t.command = 'prog -o {out} --extra {extra} -e {evalue}'
	assert t.argv(out="/tmp/o", extra="-w 3", evalue=1e-5) == [
		"prog", "-o", "/tmp/o", "--extra", "-w", "3", "-e", "1e-05"]
	assert t.argv(out="/tmp/o", extra="", evalue=1) == ["prog", "-o", "/tmp/o", "--extra", "-e", "1"]
	assert t.missing({"out": "", "evalue": 1}) == ["out"]


def test_round_trip_write(tmp_path):
	tools = Tools.load()
	tools.write(tmp_path / "t.tsv")
	again = Tools.load(tmp_path / "t.tsv")
	assert again["genomad"].command == tools["genomad"].command
	assert again["jackhmmer"].options == tools["jackhmmer"].options


def test_shipped_template_survives_user_override(tmp_path):
	user = tmp_path / "tools_table.tsv"
	user.write_text("#name\tcommand\ndeeptmhmm\t/old/python predict.py --fasta {fasta} --output-dir {out}\n")
	tools = Tools.load(user)
	assert tools["deeptmhmm"].command.startswith("/old/python")
	assert tools.template("deeptmhmm") == "{fasta} {out} --simplify-io"
	assert tools.template("mafft").startswith("--auto")


def test_user_table_holds_only_overrides(tmp_path):
	tools = Tools.load()
	tools["mafft"].command = "/opt/bin/mafft --auto {in}"
	tools.write(tmp_path / "t.tsv")
	text = (tmp_path / "t.tsv").read_text()
	assert "mafft" in text and "jackhmmer" not in text and "incE" not in text
	again = Tools.load(tmp_path / "t.tsv")
	assert again["mafft"].command == "/opt/bin/mafft --auto {in}"
	assert again["jackhmmer"].options["incE"] == "1e-10"
	stale = tmp_path / "stale.tsv"
	stale.write_text("#name\tcommand\tdirectory\tscan_range\tengine\toptions\njackhmmer\t\t\t\tjackhmmer\titerations=3;incE=1e-3\ngenomad\t/x/genomad end-to-end\t\t\t\tdb=/db\n")
	tools = Tools.load(stale)
	assert tools["jackhmmer"].options["incE"] == "1e-3"
	tools.write(stale, installer=True)
	tools = Tools.load(stale)
	assert tools["jackhmmer"].options["incE"] == "1e-10" and tools["genomad"].options["db"] == "/db"
	assert tools["genomad"].command == "/x/genomad end-to-end"
