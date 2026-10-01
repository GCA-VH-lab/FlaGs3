import threading
import time
import traceback
from pathlib import Path
from typing import Iterable

from flags3.run import Config, RunDir, RunError, Status


class StageError(Exception):
	pass


class Stage:
	name = ""
	requires: tuple[str, ...] = ()
	needs: tuple[str, ...] = ()
	optional = False
	background = False
	barrier = False

	def run(self, run: RunDir, config: Config, out: Path) -> None:
		raise NotImplementedError

	def wanted(self, config: Config) -> bool:
		return True


class Runner:
	def __init__(self, run: RunDir, config: Config, report=print):
		self.run = run
		self.config = config
		self.report = report

	def missing(self, stage: Stage) -> list[str]:
		return [name for name in stage.requires if not self.run.has(name)]

	def execute(self, stage: Stage) -> bool:
		missing = self.missing(stage)
		if missing:
			raise StageError("{} needs {} first".format(
				stage.name, ", ".join(missing)))
		out = self.run.reset(stage.name)
		status = self.run.status(stage.name)
		status.set("status", Status.RUNNING)
		status.set("started", time.strftime("%Y-%m-%d %H:%M:%S"))
		status.save()
		self.report(">> {}".format(stage.name))
		clock = time.perf_counter()
		try:
			stage.run(self.run, self.config, out)
		except Exception as error:
			status.set("status", Status.FAILED)
			status.set("seconds", "{:.2f}".format(time.perf_counter() - clock))
			status.set("error", "{}: {}".format(type(error).__name__, error))
			status.save()
			with open(out / "traceback.txt", "w") as handle:
				handle.write(traceback.format_exc())
			if not stage.optional:
				raise
			self.report("Warning: {} failed ({}); the run continues without it.".format(
				stage.name, error))
			return False
		status.set("status", Status.OK)
		status.set("seconds", "{:.2f}".format(time.perf_counter() - clock))
		status.save()
		self.report(">> {} done in {:.1f} s".format(stage.name, time.perf_counter() - clock))
		return True

	def execute_all(self, stages: Iterable[Stage]) -> dict[str, bool]:
		results: dict[str, bool] = {}
		threads: list[tuple[Stage, threading.Thread]] = []
		for stage in stages:
			if not stage.wanted(self.config):
				continue
			if stage.barrier:
				self.join(threads, results)
			if stage.background and not self.missing(stage):
				thread = threading.Thread(target=self._in_background, args=(stage, results), daemon=True)
				thread.start()
				threads.append((stage, thread))
				continue
			results[stage.name] = self.execute(stage)
		self.join(threads, results)
		return results

	def _in_background(self, stage: Stage, results: dict[str, bool]) -> None:
		try:
			results[stage.name] = self.execute(stage)
		except Exception:
			results[stage.name] = False

	def join(self, threads: list, results: dict[str, bool]) -> None:
		for stage, thread in threads:
			if thread.is_alive():
				self.report(">> waiting for {} to finish".format(stage.name))
			thread.join()
		threads.clear()
