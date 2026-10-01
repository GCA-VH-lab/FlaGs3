import shutil
import subprocess
from io import StringIO
from pathlib import Path

from Bio import Phylo

from flags3 import fasta
from flags3.log import debug, record_command
from flags3.tools import Tool, Tools

GAPS = "-."
PRESETS = ("gappyout", "strict", "strictplus", "automated1", "nogaps", "noallgaps")


class TreeError(RuntimeError):
	pass


def run(tool: Tool, argv: list[str], **kwargs) -> subprocess.CompletedProcess:
	debug("running: " + " ".join(argv))
	kwargs.setdefault("env", tool.environment())
	kwargs.setdefault("cwd", tool.directory or None)
	try:
		done = subprocess.run(argv, **kwargs)
	except OSError as error:
		record_command(argv, "not run", stderr=str(error))
		raise TreeError("could not run {}: {}".format(argv[0], error))
	record_command(argv, done.returncode, stderr=getattr(done, "stderr", None))
	if done.returncode != 0:
		raise TreeError("{} exited {}: {}".format(argv[0], done.returncode,
			(getattr(done, "stderr", "") or "").strip()[-300:]))
	return done


class TreeBuilder:
	def __init__(self, tools: Tools, work: Path, threads: int = 0, engine: str = "veryfasttree",
			trimal_mode: str = "gt", trimal_value: float = 0.1, trimal_extra: str = ""):
		self.tools = tools
		self.work = work
		self.threads = threads
		self.engine = engine
		self.trimal_mode = (trimal_mode or "gt").lstrip("-")
		self.trimal_value = trimal_value
		self.trimal_extra = trimal_extra
		self.commands: list[str] = []
		self.raw_alignment: dict[str, str] = {}
		self.alignment: dict[str, str] = {}

	def build(self, sequences: dict[str, str]) -> str:
		if len(sequences) < 3:
			raise TreeError("a tree needs at least 3 query sequences, got {}".format(len(sequences)))
		self.work.mkdir(parents=True, exist_ok=True)
		query = self.work / "queries.fasta"
		fasta.write(query, sequences.items())
		aligned = self.align(query)
		trimmed = self.trim(aligned)
		return self.iqtree(trimmed, len(sequences)) if self.engine == "iqtree" else self.veryfasttree(trimmed)

	def align(self, query: Path) -> Path:
		tool = self.tools["mafft"]
		aligned = self.work / "alignment.aln"
		argv = tool.argv(threads=self.threads or 1, **{"in": query})
		with open(aligned, "w") as out:
			run(tool, argv, stdout=out, stderr=subprocess.PIPE, text=True)
		self.commands.append(" ".join(argv))
		self.raw_alignment = {name: seq for name, _, seq in fasta.read(aligned)}
		return aligned

	def trim(self, aligned: Path) -> Path:
		tool = self.tools["trimal"]
		trimmed = self.work / "trimmed.aln"
		spec = "-" + self.trimal_mode if self.trimal_mode in PRESETS else "-{} {}".format(self.trimal_mode, self.trimal_value)
		if self.trimal_extra:
			spec += " " + self.trimal_extra
		argv = tool.argv(mode=spec, **{"in": aligned, "out": trimmed})
		if tool.locate()[0]:
			run(tool, argv, capture_output=True, text=True)
			self.commands.append(" ".join(argv))
			self.alignment = {name: seq for name, _, seq in fasta.read(trimmed)}
		else:
			self.alignment = self.trim_internally()
			fasta.write(trimmed, self.alignment.items())
		return trimmed

	def trim_internally(self) -> dict[str, str]:
		if self.trimal_mode != "gt":
			print("Warning: trimal not found and -{} has no internal equivalent; using the untrimmed alignment.".format(self.trimal_mode))
			return dict(self.raw_alignment)
		print("Warning: trimal not found; trimming columns at gt {} internally instead.".format(self.trimal_value))
		rows = list(self.raw_alignment.values())
		width = len(rows[0])
		need = float(self.trimal_value) * len(rows)
		keep = [i for i in range(width) if sum(1 for r in rows if r[i] not in GAPS) >= need]
		if not keep:
			return dict(self.raw_alignment)
		return {name: "".join(seq[i] for i in keep) for name, seq in self.raw_alignment.items()}

	def veryfasttree(self, trimmed: Path) -> str:
		tool = self.tools["veryfasttree"]
		argv = tool.argv(**{"in": trimmed})
		done = run(tool, argv, capture_output=True, text=True)
		self.commands.append(" ".join(argv))
		return done.stdout.strip()

	def iqtree(self, trimmed: Path, taxa: int) -> str:
		tool = self.tools["iqtree"]
		prefix = self.work / "iq"
		argv = tool.argv(model="MFP", prefix=prefix, threads=self.threads or "AUTO", **{"in": trimmed})
		if shutil.which(argv[0]) is None:
			alternative = next((b for b in ("iqtree3", "iqtree2", "iqtree") if shutil.which(b)), None)
			if alternative is None:
				raise TreeError("no iqtree binary found")
			argv[0] = alternative
		if taxa >= 4:
			argv += ["-B", "1000"]
		run(tool, argv, capture_output=True, text=True)
		self.commands.append(" ".join(argv))
		return (self.work / "iq.treefile").read_text().strip()


def leaf_order(newick: str, ladderize: bool = True) -> list[str]:
	tree = Phylo.read(StringIO(newick), "newick")
	if ladderize:
		try:
			tree.root_at_midpoint()
		except Exception:
			pass
		tree.ladderize()
	return [tip.name for tip in tree.get_terminals()]
