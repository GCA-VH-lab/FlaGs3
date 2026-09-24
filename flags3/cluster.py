import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

from flags3.log import debug, record_command
from flags3.tools import Tool, ToolError, Tools

ENGINES = ("jackhmmer", "nhmmer", "mmseqs")


class ClusterError(RuntimeError):
	pass


def methods(tools: Tools) -> list[str]:
	return [n for n in tools.names() if tools[n].engine in ENGINES]


def method(tools: Tools, name: str) -> Tool:
	if name not in tools or tools[name].engine not in ENGINES:
		raise ClusterError("no clustering method {!r} in the tools table; available: {}".format(
			name, ", ".join(methods(tools)) or "none"))
	return tools[name]


def option(tool: Tool, key: str, default):
	raw = tool.options.get(key)
	if raw is None:
		return default
	try:
		return type(default)(raw)
	except (TypeError, ValueError):
		raise ClusterError("{}: option {}={!r} is not a {}".format(tool.name, key, raw, type(default).__name__))


def connected_components(adjacency: dict[str, set]) -> list[list[str]]:
	parent = {node: node for node in adjacency}

	def find(node):
		root = node
		while parent[root] != root:
			root = parent[root]
		while parent[node] != root:
			parent[node], node = root, parent[node]
		return root

	for node, hits in adjacency.items():
		root = find(node)
		for hit in hits:
			if hit not in parent:
				continue
			other = find(hit)
			if other != root:
				parent[other] = root
				root = find(root)
	groups: dict[str, list[str]] = {}
	for node in parent:
		groups.setdefault(find(node), []).append(node)
	families = [sorted(group) for group in groups.values()]
	families.sort(key=lambda family: (-len(family), family[0]))
	return families


def by_name(products: dict[str, str]) -> list[list[str]]:
	groups: dict[str, list[str]] = {}
	for accession, product in sorted(products.items()):
		key = " ".join((product or "").lower().split()) or accession
		groups.setdefault(key, []).append(accession)
	families = [sorted(members) for members in groups.values()]
	families.sort(key=lambda family: (-len(family), family[0]))
	return families


class Clusterer:
	def __init__(self, tool: Tool, workers: Optional[int], work: Path):
		self.tool = tool
		self.workers = workers or os.cpu_count() or 1
		self.work = work
		self.adjacency: dict[str, set] = {}

	def cluster(self, sequences: dict[str, str]) -> list[list[str]]:
		self.adjacency = {}
		if not sequences:
			return []
		self.adjacency = self.edges(sequences)
		for name in sequences:
			self.adjacency.setdefault(name, set())
		return connected_components(self.adjacency)

	def edges(self, sequences: dict[str, str]) -> dict[str, set]:
		raise NotImplementedError


class PyhmmerClusterer(Clusterer):
	alphabet_name = "amino"

	def chunk_size(self, total: int) -> int:
		cap = option(self.tool, "chunk", 100)
		per_worker = option(self.tool, "chunks_per_worker", 16)
		return max(1, min(cap, total // (max(self.workers, 1) * max(per_worker, 1))))

	def edges(self, sequences: dict[str, str]) -> dict[str, set]:
		from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
		alphabet = Alphabet.amino() if self.alphabet_name == "amino" else Alphabet.dna()
		digital = [(name, TextSequence(name=name.encode(), sequence=seq).digitize(alphabet))
			for name, seq in sorted(sequences.items())]
		block = DigitalSequenceBlock(alphabet, [d for _, d in digital])
		size = self.chunk_size(len(digital))
		chunks = [digital[at:at + size] for at in range(0, len(digital), size)]

		def run(chunk):
			results = self.search([query for _, query in chunk], block)
			return [(name, {_name(h) for h in self.hits_of(result) if h.included})
				for (name, _), result in zip(chunk, results)]

		adjacency: dict[str, set] = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			for part in pool.map(run, chunks):
				adjacency.update(part)
		return adjacency

	def search(self, queries, block):
		raise NotImplementedError

	@staticmethod
	def hits_of(result):
		return result


class JackhmmerClusterer(PyhmmerClusterer):
	@staticmethod
	def hits_of(result):
		return result.hits

	def search(self, queries, block):
		import pyhmmer
		return pyhmmer.hmmer.jackhmmer(queries, block, max_iterations=option(self.tool, "iterations", 3),
			incE=option(self.tool, "incE", 1e-3), cpus=1)


class NhmmerClusterer(PyhmmerClusterer):
	alphabet_name = "dna"

	def search(self, queries, block):
		import pyhmmer
		return pyhmmer.hmmer.nhmmer(queries, block, incE=option(self.tool, "incE", 1e-3), cpus=1)


class MmseqsClusterer(Clusterer):
	def edges(self, sequences: dict[str, str]) -> dict[str, set]:
		self.work.mkdir(parents=True, exist_ok=True)
		fasta = self.work / "input.fasta"
		hits = self.work / "hits.tsv"
		with open(fasta, "w") as handle:
			for name, seq in sorted(sequences.items()):
				handle.write(">{}\n{}\n".format(name, seq))
		self.run_search(fasta, hits, len(sequences))
		return self.read_pairs(hits, sequences)

	def command_line(self, fasta: Path, hits: Path, total: int) -> list[str]:
		if not self.tool.command:
			raise ClusterError("{}: no command in the tools table".format(self.tool.name))
		values = dict(self.tool.options)
		values.update({"in": fasta, "out": hits, "tmp": self.work / "tmp",
			"threads": self.workers, "maxseqs": max(total, 300)})
		command = self.tool.argv(**values)
		if any("{" in part for part in command):
			raise ClusterError("{}: command needs an option the table does not supply: {}".format(
				self.tool.name, " ".join(p for p in command if "{" in p)))
		return command

	def run_search(self, fasta: Path, hits: Path, total: int) -> None:
		command = self.command_line(fasta, hits, total)
		found, where = self.tool.locate()
		if not found:
			raise ClusterError("{} for clustering method {!r}: {}".format(command[0], self.tool.name, where))
		debug("cluster {}: {}".format(self.tool.name, " ".join(command)))
		try:
			done = subprocess.run(command, capture_output=True, text=True, env=self.tool.environment())
		except OSError as error:
			raise ClusterError("{}: {}".format(self.tool.name, error))
		record_command(command, done.returncode, done.stdout, done.stderr)
		if done.returncode != 0 or not hits.is_file():
			raise ClusterError("{} exited {}: {}".format(command[0], done.returncode,
				(done.stderr or done.stdout or "").strip()[-300:]))

	@staticmethod
	def read_pairs(hits: Path, sequences: dict[str, str]) -> dict[str, set]:
		adjacency: dict[str, set] = {name: set() for name in sequences}
		with open(hits, encoding="utf-8") as handle:
			for line in handle:
				cells = line.rstrip("\n").split("\t")
				if len(cells) >= 2 and cells[0] in adjacency and cells[1] in adjacency:
					adjacency[cells[0]].add(cells[1])
		return adjacency


ENGINE_CLASSES = {"jackhmmer": JackhmmerClusterer, "nhmmer": NhmmerClusterer, "mmseqs": MmseqsClusterer}


def build(tools: Tools, name: str, workers: Optional[int], work: Path) -> Clusterer:
	tool = method(tools, name)
	return ENGINE_CLASSES[tool.engine](tool, workers, work)


def _name(hit) -> str:
	raw = hit.name
	return raw.decode() if isinstance(raw, bytes) else raw
