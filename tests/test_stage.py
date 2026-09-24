import pytest

from flags3.run import RunDir
from flags3.stage import Runner, Stage, StageError


class Extract(Stage):
	name = "extract"

	def run(self, run, config, out):
		(out / "genes.tsv").write_text("row_id\n")


class Domains(Stage):
	name = "domains"
	requires = ("extract",)
	optional = True

	def wanted(self, config):
		return config.flag("domains")

	def run(self, run, config, out):
		raise RuntimeError("no hmm database")


class Report(Stage):
	name = "report"
	requires = ("extract", "domains")


def _runner(tmp_path, **config):
	run = RunDir(tmp_path / "out").create("3.0", "")
	cfg = run.config()
	cfg.update(config)
	cfg.save()
	messages = []
	return Runner(run, cfg, report=messages.append), messages


def test_requirements_are_checked(tmp_path):
	runner, _ = _runner(tmp_path)
	with pytest.raises(StageError):
		runner.execute(Report())


def test_ok_stage_records_status(tmp_path):
	runner, messages = _runner(tmp_path)
	assert runner.execute(Extract()) is True
	status = runner.run.status("extract")
	assert status.ok
	assert float(status.get("seconds")) >= 0
	assert runner.run.stage_file("extract", "genes.tsv").is_file()
	assert messages[0] == ">> extract"


def test_optional_failure_is_recorded_and_survived(tmp_path):
	runner, messages = _runner(tmp_path, domains=True)
	runner.execute(Extract())
	assert runner.execute(Domains()) is False
	status = runner.run.status("domains")
	assert status.state == "failed"
	assert "no hmm database" in status.get("error")
	assert runner.run.stage_file("domains", "traceback.txt").is_file()
	assert any(m.startswith("Warning: domains failed") for m in messages)


def test_required_failure_raises(tmp_path):
	runner, _ = _runner(tmp_path)
	failing = Extract()
	failing.run = lambda run, config, out: 1 / 0
	with pytest.raises(ZeroDivisionError):
		runner.execute(failing)
	assert runner.run.status("extract").state == "failed"


def test_execute_all_skips_unwanted(tmp_path):
	runner, _ = _runner(tmp_path, domains=False)
	results = runner.execute_all([Extract(), Domains()])
	assert results == {"extract": True}
	assert not runner.run.stage_dir("domains").exists()


def test_background_stage_runs_alongside_and_is_joined_at_barrier(tmp_path):
	import threading
	import time
	runner, messages = _runner(tmp_path)
	order = []
	gate = threading.Event()

	class Cloud(Stage):
		name = "cloud"
		requires = ("extract",)
		background = True
		optional = True

		def run(self, run, config, out):
			order.append("cloud-start")
			gate.wait(5)
			order.append("cloud-end")

	class Local(Stage):
		name = "local"
		requires = ("extract",)

		def run(self, run, config, out):
			order.append("local")
			gate.set()

	class Final(Stage):
		name = "final"
		requires = ("extract",)
		barrier = True

		def run(self, run, config, out):
			order.append("final")

	results = runner.execute_all([Extract(), Cloud(), Local(), Final()])
	assert results == {"extract": True, "cloud": True, "local": True, "final": True}
	assert order.index("local") < order.index("cloud-end") < order.index("final")
	assert order.index("cloud-start") < order.index("cloud-end")
