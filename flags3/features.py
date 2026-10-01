import os
import re
import subprocess
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from flags3.log import debug, note, record_command
from flags3.tools import Tool, brief

RESULTS_URL = "https://biolib.com/results/{}/"


class FeatureError(RuntimeError):
	pass


@dataclass(frozen=True)
class Region:
	protein: str
	kind: str
	start: int
	end: int


def runs(topology: str, wanted: str, kind: str, protein: str) -> list[Region]:
	out, start = [], None
	for i, code in enumerate(topology):
		if code in wanted:
			if start is None:
				start = i
		elif start is not None:
			out.append(Region(protein, kind, start + 1, i))
			start = None
	if start is not None:
		out.append(Region(protein, kind, start + 1, len(topology)))
	return out


def parse_3line(text: str, want_signal: bool) -> list[Region]:
	lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
	regions = []
	i = 0
	while i < len(lines):
		if not lines[i].startswith(">"):
			i += 1
			continue
		protein = lines[i][1:].split("|")[0].split()[0].strip()
		topology = lines[i + 2] if i + 2 < len(lines) else ""
		regions += runs(topology, "MB", "tm", protein)
		if want_signal:
			regions += runs(topology, "S", "signal", protein)
		i += 3
	return regions


def parse_signalp(text: str) -> list[Region]:
	regions = []
	position = re.compile(r"CS pos:\s*(\d+)")
	for line in text.splitlines():
		if line.startswith("#") or not line.strip():
			continue
		col = line.split("\t")
		if len(col) < 2 or col[1].strip().upper() in ("", "OTHER"):
			continue
		match = position.search(line)
		if match:
			regions.append(Region(col[0].split()[0].strip(), "signal", 1, int(match.group(1))))
	return regions


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
	with open(path, "w") as out:
		for name, seq in sequences.items():
			out.write(">{}\n{}\n".format(name, seq))


class LocalRunner:
	def __init__(self, tool: Tool, work: Path):
		self.tool = tool
		self.work = work
		self.text = ""

	def start(self, sequences: dict[str, str], result_suffix: str) -> None:
		self.work.mkdir(parents=True, exist_ok=True)
		fasta = self.work / "query.fasta"
		out_dir = self.work / "out"
		write_fasta(fasta, sequences)
		argv = self.tool.argv(fasta=fasta, out=out_dir)
		found, where = self.tool.locate()
		if not found and self.tool.directory and (Path(self.tool.directory) / argv[0]).is_file():
			found = True
		if not found:
			raise FeatureError("{}: {}".format(self.tool.name, where))
		debug("features: " + " ".join(argv))
		done = subprocess.run(argv, cwd=self.tool.directory or None, capture_output=True, text=True, env=self.tool.environment())
		record_command(argv, done.returncode, done.stdout, done.stderr)
		if done.returncode != 0:
			raise FeatureError("{} exited {}: {}".format(argv[0], done.returncode, brief(done.stderr or done.stdout)))
		for root in (out_dir, self.work, Path(self.tool.directory) if self.tool.directory else self.work):
			for path in sorted(root.rglob("*" + result_suffix)) if root.is_dir() else []:
				self.text = path.read_text()
				return
		raise FeatureError("{} produced no {}".format(argv[0], result_suffix))

	def finish(self) -> str:
		return self.text


class BioLibRunner:
	LOCK = threading.Lock()

	def __init__(self, app: str, args_template: str, batch_size: int):
		try:
			import biolib
		except ImportError:
			raise FeatureError("pybiolib is not installed; pip install pybiolib, or use the local tool")
		self.biolib = biolib
		self.app = app
		self.args_template = args_template
		self.batch_size = batch_size
		self.jobs: list[tuple[object, Path]] = []
		self.result_suffix = ""

	def start(self, sequences: dict[str, str], result_suffix: str) -> None:
		self.result_suffix = result_suffix
		application = self.biolib.load(self.app)
		items = list(sequences.items())
		for at in range(0, len(items), self.batch_size):
			tmp = Path(tempfile.mkdtemp(prefix="flags3_biolib_"))
			write_fasta(tmp / "query.fasta", dict(items[at:at + self.batch_size]))
			with self.LOCK:
				previous = os.getcwd()
				os.chdir(tmp)
				try:
					job = application.cli(args=self.args_template.format(fasta="query.fasta"), blocking=False)
				finally:
					os.chdir(previous)
			self.jobs.append((job, tmp))
			note("{}: {} sequences submitted as job {}, follow it at {}".format(
				self.app, min(self.batch_size, len(items) - at), _job_id(job), RESULTS_URL.format(_job_id(job))))

	def finish(self) -> str:
		texts = []
		for job, tmp in self.jobs:
			if hasattr(job, "wait"):
				job.wait()
			job.save_files(str(tmp / "out"))
			path = next(iter(sorted((tmp / "out").rglob("*" + self.result_suffix))), None)
			if path is None:
				raise FeatureError("{} job {} produced no {}".format(self.app, _job_id(job), self.result_suffix))
			texts.append(path.read_text())
		return "\n".join(texts)


def _job_id(job) -> str:
	for attr in ("id", "job_id", "uuid"):
		value = getattr(job, attr, None)
		if value:
			return str(value)
	return "?"


class TmScanner:
	APP = "DTU/DeepTMHMM2"
	ARGS = "{fasta} results --simplify-io"
	LEGACY_APP = "DTU/DeepTMHMM"
	LEGACY_ARGS = "--fasta {fasta}"
	RESULT = "predicted_topologies.3line"
	BATCH = 2000

	@classmethod
	def remote_args(cls, app: str) -> str:
		return cls.LEGACY_ARGS if app == cls.LEGACY_APP else cls.ARGS

	def __init__(self, runner, want_signal: bool):
		self.runner = runner
		self.want_signal = want_signal

	def start(self, sequences: dict[str, str]) -> None:
		if sequences:
			self.runner.start(sequences, self.RESULT)

	def finish(self) -> list[Region]:
		return parse_3line(self.runner.finish(), self.want_signal)


class SignalPScanner:
	APP = "DTU/SignalP-6"
	ARGS = "--fastafile {fasta} --output_dir output --organism other --format txt --mode fast"
	RESULT = "prediction_results.txt"
	BATCH = 1000

	def __init__(self, runner):
		self.runner = runner

	def start(self, sequences: dict[str, str]) -> None:
		if sequences:
			self.runner.start(sequences, self.RESULT)

	def finish(self) -> list[Region]:
		return parse_signalp(self.runner.finish())
