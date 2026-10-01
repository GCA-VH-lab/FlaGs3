import pytest

from flags3.run import Config, RunDir, RunError


def test_create_and_reopen(tmp_path):
	run = RunDir(tmp_path / "out").create("3.0", "flags3 run -i x")
	assert run.exists()
	assert run.input_dir.is_dir()
	info = RunDir(tmp_path / "out").require().info()
	assert info.get("version") == "3.0"
	with pytest.raises(RunError):
		run.create("3.0", "again")
	with pytest.raises(RunError):
		RunDir(tmp_path / "nowhere").require()


def test_config_types(tmp_path):
	run = RunDir(tmp_path / "out").create("3.0", "")
	cfg = run.config()
	cfg.update({"gene": 4, "domains": True, "range": None, "hmmdb": ["a.hmm", "b.hmm"]})
	cfg.set("blast.evalue", 1e-5)
	cfg.save()
	cfg = run.config()
	assert cfg.integer("gene") == 4
	assert cfg.flag("domains") is True
	assert cfg.integer("range") is None
	assert cfg.text("hmmdb") == "a.hmm,b.hmm"
	assert cfg.number("blast.evalue") == 1e-5
	assert cfg.items("blast.") == {"evalue": "1e-05"}
	assert cfg.path("missing") is None


def test_stage_dirs_and_status(tmp_path):
	run = RunDir(tmp_path / "out").create("3.0", "")
	assert run.stages() == []
	assert not run.has("extract")
	out = run.reset("extract")
	(out / "junk").write_text("x")
	status = run.status("extract")
	status.set("status", "ok")
	status.save()
	assert run.has("extract")
	assert run.stages() == ["extract"]
	run.reset("extract")
	assert not (out / "junk").exists()
	assert not run.has("extract")
