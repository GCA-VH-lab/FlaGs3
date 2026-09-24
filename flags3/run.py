import shutil
import time
from pathlib import Path
from typing import Iterator, Optional

RUN_FILE = "run.tsv"
CONFIG_FILE = "config.tsv"
STATUS_FILE = "status.tsv"
CONSOLE_LOG = "console.log"
INPUT_DIR = "input"
FIGURES_DIR = "figures"
REPORT_DIR = "report"


class RunError(Exception):
	pass


class KeyValueFile:
	def __init__(self, path: Path):
		self.file = path
		self.values: dict[str, str] = {}

	def load(self) -> "KeyValueFile":
		self.values = {}
		if self.file.is_file():
			with open(self.file, encoding="utf-8") as handle:
				for line in handle:
					line = line.rstrip("\n")
					if not line or line.startswith("#"):
						continue
					key, _, value = line.partition("\t")
					self.values[key] = value
		return self

	def save(self) -> None:
		self.file.parent.mkdir(parents=True, exist_ok=True)
		with open(self.file, "w", encoding="utf-8") as handle:
			for key in self.values:
				handle.write("{}\t{}\n".format(key, self.values[key]))

	def __contains__(self, key: str) -> bool:
		return key in self.values

	def __iter__(self) -> Iterator[str]:
		return iter(self.values)

	def get(self, key: str, default: Optional[str] = None) -> Optional[str]:
		return self.values.get(key, default)

	def set(self, key: str, value) -> None:
		self.values[key] = _stringify(value)

	def update(self, other: dict) -> None:
		for key, value in other.items():
			self.set(key, value)


class Config(KeyValueFile):
	def text(self, key: str, default: str = "") -> str:
		return self.values.get(key, default)

	def integer(self, key: str, default: Optional[int] = None) -> Optional[int]:
		raw = self.values.get(key)
		return default if raw in (None, "") else int(raw)

	def number(self, key: str, default: Optional[float] = None) -> Optional[float]:
		raw = self.values.get(key)
		return default if raw in (None, "") else float(raw)

	def flag(self, key: str, default: bool = False) -> bool:
		raw = self.values.get(key)
		return default if raw is None else raw.lower() in ("true", "1", "yes")

	def path(self, key: str) -> Optional[Path]:
		raw = self.values.get(key)
		return Path(raw) if raw else None

	def workers(self) -> int:
		import os
		return self.integer("cpu") or os.cpu_count() or 1

	def items(self, prefix: str) -> dict[str, str]:
		cut = len(prefix)
		return {k[cut:]: v for k, v in self.values.items() if k.startswith(prefix)}


class Status(KeyValueFile):
	OK = "ok"
	FAILED = "failed"
	RUNNING = "running"

	@property
	def state(self) -> Optional[str]:
		return self.values.get("status")

	@property
	def ok(self) -> bool:
		return self.state == self.OK


class RunDir:
	def __init__(self, path):
		self.path = Path(path).resolve()

	@property
	def name(self) -> str:
		return self.path.name

	@property
	def run_file(self) -> Path:
		return self.path / RUN_FILE

	@property
	def config_file(self) -> Path:
		return self.path / CONFIG_FILE

	@property
	def console_log(self) -> Path:
		return self.path / CONSOLE_LOG

	@property
	def input_dir(self) -> Path:
		return self.path / INPUT_DIR

	@property
	def figures_dir(self) -> Path:
		return self.path / FIGURES_DIR

	@property
	def report_dir(self) -> Path:
		return self.path / REPORT_DIR

	def exists(self) -> bool:
		return self.run_file.is_file()

	def create(self, version: str, command: str) -> "RunDir":
		if self.exists():
			raise RunError("{} is already a FlaGs3 run".format(self.path))
		self.path.mkdir(parents=True, exist_ok=True)
		self.input_dir.mkdir(exist_ok=True)
		info = KeyValueFile(self.run_file)
		info.set("version", version)
		info.set("started", _now())
		info.set("command", command)
		info.save()
		return self

	def require(self) -> "RunDir":
		if not self.exists():
			raise RunError("{} is not a FlaGs3 run: no {}".format(self.path, RUN_FILE))
		return self

	def info(self) -> KeyValueFile:
		return KeyValueFile(self.run_file).load()

	def config(self) -> Config:
		return Config(self.config_file).load()

	def stage_dir(self, stage: str) -> Path:
		return self.path / stage

	def stage_file(self, stage: str, name: str) -> Path:
		return self.stage_dir(stage) / name

	def status(self, stage: str) -> Status:
		return Status(self.stage_file(stage, STATUS_FILE)).load()

	def has(self, stage: str) -> bool:
		return self.status(stage).ok

	def reset(self, stage: str) -> Path:
		target = self.stage_dir(stage)
		if target.exists():
			shutil.rmtree(target)
		target.mkdir(parents=True)
		return target

	def stages(self) -> list[str]:
		found = []
		for child in sorted(self.path.iterdir()):
			if child.is_dir() and (child / STATUS_FILE).is_file():
				found.append(child.name)
		return found


def _stringify(value) -> str:
	if value is None:
		return ""
	if isinstance(value, bool):
		return "true" if value else "false"
	if isinstance(value, (list, tuple)):
		return ",".join(_stringify(v) for v in value)
	return str(value)


def _now() -> str:
	return time.strftime("%Y-%m-%d %H:%M:%S")
