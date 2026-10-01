import csv
import os
import shlex
import shutil
from importlib import resources
from pathlib import Path
from typing import Optional

TABLE_NAME = "tools_table.tsv"
COLUMNS = ("name", "command", "directory", "scan_range", "engine", "options")


class ToolError(Exception):
	pass


class Tool:
	def __init__(self, name: str):
		self.name = name
		self.command = ""
		self.directory = ""
		self.engine = ""
		self.options: dict[str, str] = {}
		self.scan_range: Optional[int] = None

	def update(self, row: dict[str, str], origin: str) -> None:
		if row.get("command"):
			self.command = row["command"]
		if row.get("directory"):
			self.directory = os.path.expanduser(row["directory"])
		if row.get("engine"):
			self.engine = row["engine"]
		if row.get("options"):
			self.options = self.parse_options(row["options"])
		span = row.get("scan_range", "")
		if span:
			if span.lower() == "genome":
				self.scan_range = 0
			elif span.isdigit():
				self.scan_range = int(span)
			else:
				raise ToolError("{}: scan_range for {} must be a number of bases or 'genome', got {!r}".format(
					origin, self.name, span))

	@staticmethod
	def parse_options(text: str) -> dict[str, str]:
		out = {}
		for item in text.split(";"):
			key, sep, value = item.partition("=")
			if sep:
				out[key.strip()] = value.strip()
		return out

	def argv(self, **values) -> list[str]:
		out = []
		for part in shlex.split(self.command):
			filled = part
			for key, value in values.items():
				filled = filled.replace("{%s}" % key, str(value))
			if not filled:
				continue
			if filled != part and " " in filled and part.startswith("{"):
				out.extend(shlex.split(filled))
			else:
				out.append(filled)
		return out

	def missing(self, values: dict) -> list[str]:
		return [k for k in values
			if "{" + k + "}" in self.command and not str(values.get(k, "")).strip()]

	def locate(self) -> tuple[bool, str]:
		if not self.command:
			return False, "no command"
		program = shlex.split(self.command)[0]
		if os.path.isabs(program):
			if os.path.isfile(program) and os.access(program, os.X_OK):
				return True, program
			return False, "{} is not an executable file".format(program)
		found = shutil.which(program)
		if found:
			return True, found
		return False, "{} not found on PATH".format(program)

	def environment(self) -> Optional[dict]:
		if not self.command:
			return None
		program = shlex.split(self.command)[0]
		if not os.path.isabs(program):
			return None
		bindir = os.path.dirname(program)
		if not os.path.isdir(bindir):
			return None
		env = dict(os.environ)
		parts = [p for p in env.get("PATH", "").split(os.pathsep) if p]
		if bindir in parts:
			return None
		env["PATH"] = os.pathsep.join([bindir] + parts)
		return env


class Tools:
	def __init__(self):
		self.tools: dict[str, Tool] = {}
		self.origins: list[str] = []
		self.shipped: dict[str, str] = {}

	@classmethod
	def load(cls, user_table: Optional[Path] = None) -> "Tools":
		tools = cls()
		with resources.as_file(resources.files("flags3.data") / TABLE_NAME) as shipped:
			tools.read(shipped, allow_new=True)
		tools.shipped = {name: tool.command for name, tool in tools.tools.items()}
		if user_table is not None:
			if not Path(user_table).is_file():
				raise ToolError("tool table not found: {}".format(user_table))
			tools.read(Path(user_table), allow_new=False)
		return tools

	def read(self, path: Path, allow_new: bool) -> None:
		with open(path, newline="", encoding="utf-8") as handle:
			lines = (ln for ln in handle if ln.strip() and not ln.startswith("##"))
			for raw in csv.DictReader(lines, delimiter="\t"):
				row = {(k or "").lstrip("#").strip(): (v or "").strip() for k, v in raw.items()}
				name = row.get("name", "").lower()
				if not name:
					continue
				if name not in self.tools:
					if not allow_new:
						raise ToolError("{}: unknown tool {!r}; known tools: {}".format(
							path, name, ", ".join(sorted(self.tools))))
					self.tools[name] = Tool(name)
				self.tools[name].update(row, str(path))
		self.origins.append(str(path))

	def __contains__(self, name: str) -> bool:
		return name in self.tools

	def __getitem__(self, name: str) -> Tool:
		try:
			return self.tools[name]
		except KeyError:
			raise ToolError("no tool named {!r}; known tools: {}".format(
				name, ", ".join(sorted(self.tools))))

	def template(self, name: str) -> str:
		command = self.shipped.get(name) or self[name].command
		return command.partition(" ")[2]

	def names(self, engine: Optional[str] = None) -> list[str]:
		if engine is None:
			return sorted(self.tools)
		return sorted(n for n, t in self.tools.items() if t.engine)

	def write(self, path: Path) -> None:
		with open(path, "w", encoding="utf-8") as out:
			out.write("#" + "\t".join(COLUMNS) + "\n")
			for name in sorted(self.tools):
				t = self.tools[name]
				span = "" if t.scan_range is None else ("genome" if t.scan_range == 0 else str(t.scan_range))
				options = ";".join("{}={}".format(k, v) for k, v in t.options.items())
				out.write("\t".join((name, t.command, t.directory, span, t.engine, options)) + "\n")


def brief(text, lines: int = 3, limit: int = 300) -> str:
	kept = [ln.strip() for ln in str(text or "").splitlines() if ln.strip()]
	if not kept:
		return "no output"
	joined = " | ".join(kept[-lines:])
	return joined if len(joined) <= limit else "..." + joined[-limit:]
