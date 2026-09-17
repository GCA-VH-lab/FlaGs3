import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Sequence

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence

from flags_log import debug

ENGINES = ("jackhmmer", "nhmmer", "mmseqs")

import flags_tools


class MethodError(RuntimeError):
	pass


class Method:
	def __init__(self, name: str, engine: str, command: str, options: str):
		self.name = name
		self.engine = engine
		self.command = command.strip()
		self.options = self._parse(options)

	@staticmethod
	def _parse(text: str) -> Dict[str, str]:
		out: Dict[str, str] = {}
		for item in re.split(r"[;\s]+", (text or "").strip()):
			if not item:
				continue
			if "=" not in item:
				raise MethodError("option {!r} is not key=value".format(item))
			key, value = item.split("=", 1)
			out[key.strip()] = value.strip()
		return out

	def number(self, key: str, default):
		raw = self.options.get(key)
		if raw is None:
			return default
		try:
			return type(default)(raw)
		except (TypeError, ValueError):
			raise MethodError("{}: option {}={!r} is not a {}".format(
				self.name, key, raw, type(default).__name__))


def names() -> List[str]:
	return [n for n in flags_tools.clustering_names()
			if flags_tools.clustering(n).get("engine") in ENGINES]


def method(name: str) -> Method:
	entry = flags_tools.clustering(name)
	engine = entry.get("engine", "")
	if engine not in ENGINES:
		raise MethodError(
			"no clustering method {!r} in tools_table.tsv; available: {}".format(
				name, ", ".join(names()) or "none"))
	command, _ = flags_tools.get(name)
	return Method(name, engine, command or "", entry.get("options", ""))


def connected_components(adjacency: Dict[str, set]) -> List[List[str]]:
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

	groups: Dict[str, List[str]] = {}
	for node in parent:
		groups.setdefault(find(node), []).append(node)
	families = [sorted(group) for group in groups.values()]
	families.sort(key=lambda family: (-len(family), family[0]))
	return families


def hit_name(hit) -> str:
	raw = hit.name
	return raw.decode() if isinstance(raw, bytes) else raw


class Clusterer:
	alphabet_name = "amino"

	def __init__(self, setting: Method, workers: Optional[int] = None,
				 work_dir: Optional[str] = None):
		self.setting = setting
		self.workers = workers
		self.work_dir = work_dir
		self.alphabet = (Alphabet.amino() if self.alphabet_name == "amino"
						 else Alphabet.dna())
		self.adjacency: Dict[str, set] = {}
		self.note = ""

	@property
	def name(self) -> str:
		return self.setting.name

	def cluster(self, sequences: Dict[str, str]) -> List[List[str]]:
		self.adjacency = {}
		if not sequences:
			return []
		self.adjacency = self.edges(sequences)
		for name in sequences:
			self.adjacency.setdefault(name, set())
		return connected_components(self.adjacency)

	def edges(self, sequences: Dict[str, str]) -> Dict[str, set]:
		raise NotImplementedError

	def digitize(self, sequences: Dict[str, str]):
		return {name: TextSequence(name=name.encode(),
								   sequence=seq).digitize(self.alphabet)
				for name, seq in sorted(sequences.items())}


class PyhmmerClusterer(Clusterer):
	def search(self, queries, block):
		raise NotImplementedError

	@staticmethod
	def hits_of(result):
		return result

	def chunk_size(self, total: int) -> int:
		cap = self.setting.number("chunk", 100)
		per_worker = self.setting.number("chunks_per_worker", 16)
		slots = max(self.workers or 1, 1) * max(per_worker, 1)
		return max(1, min(cap, total // slots))

	def edges(self, sequences: Dict[str, str]) -> Dict[str, set]:
		digital = self.digitize(sequences)
		block = DigitalSequenceBlock(self.alphabet, list(digital.values()))
		items = sorted(digital.items())
		size = self.chunk_size(len(items))
		chunks = [items[at:at + size] for at in range(0, len(items), size)]

		def run(chunk):
			labels = [name for name, _ in chunk]
			results = self.search([query for _, query in chunk], block)
			return [(label, {hit_name(h) for h in self.hits_of(result)
							 if h.included})
					for label, result in zip(labels, results)]

		adjacency: Dict[str, set] = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			for part in pool.map(run, chunks):
				adjacency.update(part)
		return adjacency


class JackhmmerClusterer(PyhmmerClusterer):
	@staticmethod
	def hits_of(result):
		return result.hits

	def search(self, queries, block):
		return pyhmmer.hmmer.jackhmmer(
			queries, block,
			max_iterations=self.setting.number("iterations", 3),
			incE=self.setting.number("incE", 1e-3),
			cpus=1)


class NhmmerClusterer(PyhmmerClusterer):
	alphabet_name = "dna"

	def search(self, queries, block):
		return pyhmmer.hmmer.nhmmer(
			queries, block,
			incE=self.setting.number("incE", 1e-3),
			cpus=1)


class MmseqsClusterer(Clusterer):
	def edges(self, sequences: Dict[str, str]) -> Dict[str, set]:
		work = self.work_dir or tempfile.mkdtemp(prefix="flags3_cluster_")
		os.makedirs(work, exist_ok=True)
		fasta = os.path.join(work, "input.fasta")
		hits = os.path.join(work, "hits.tsv")
		try:
			with open(fasta, "w") as handle:
				for name, seq in sorted(sequences.items()):
					handle.write(">{}\n{}\n".format(name, seq))
			self.run_search(fasta, hits, work, len(sequences))
			return self.read_pairs(hits, sequences)
		finally:
			if not self.work_dir:
				shutil.rmtree(work, ignore_errors=True)

	def command_line(self, fasta: str, hits: str, work: str,
					 total: int) -> List[str]:
		if not self.setting.command:
			raise MethodError("{}: no command in tools_table.tsv".format(self.name))
		values = dict(self.setting.options)
		values.update({"in": fasta, "out": hits,
					   "tmp": os.path.join(work, "tmp"),
					   "threads": self.workers or 1,
					   "maxseqs": max(total, 300)})
		try:
			rendered = self.setting.command.format(**values)
		except KeyError as missing:
			raise MethodError("{}: command needs {} but no option supplies it"
							  .format(self.name, missing))
		return rendered.split()

	def run_search(self, fasta: str, hits: str, work: str, total: int):
		command = self.command_line(fasta, hits, work, total)
		binary = command[0]
		if not (os.path.isfile(binary) or shutil.which(binary)):
			raise MethodError(
				"{} not found for clustering method {!r}. Run "
				"mmseqs_installer.sh, or point tools_table.local.tsv at the "
				"binary.".format(binary, self.name))
		debug("cluster {}: {}".format(self.name, " ".join(command)))
		try:
			done = subprocess.run(command, capture_output=True, text=True)
		except OSError as error:
			raise MethodError("{}: {}".format(self.name, error))
		self.record(command, done)
		if done.returncode != 0 or not os.path.isfile(hits):
			raise MethodError("{} exited {}: {}".format(
				binary, done.returncode,
				(done.stderr or done.stdout or "").strip()[-300:]))

	@staticmethod
	def record(command: Sequence[str], done):
		try:
			import flags_log
			flags_log.record_command(list(command), done.returncode,
									 done.stdout, done.stderr)
		except (ImportError, AttributeError):
			pass

	@staticmethod
	def read_pairs(hits: str, sequences: Dict[str, str]) -> Dict[str, set]:
		adjacency: Dict[str, set] = {name: set() for name in sequences}
		with open(hits, encoding="utf-8") as handle:
			for line in handle:
				cells = line.rstrip("\n").split("\t")
				if len(cells) < 2:
					continue
				query, target = cells[0], cells[1]
				if query in adjacency and target in adjacency:
					adjacency[query].add(target)
		return adjacency


ENGINE_CLASSES = {
	"jackhmmer": JackhmmerClusterer,
	"nhmmer": NhmmerClusterer,
	"mmseqs": MmseqsClusterer,
}


def build(name: str, workers: Optional[int] = None,
		  work_dir: Optional[str] = None) -> Clusterer:
	setting = method(name)
	return ENGINE_CLASSES[setting.engine](setting, workers, work_dir)


def cluster_by_name(products: Dict[str, str]) -> List[List[str]]:
	groups: Dict[str, List[str]] = {}
	for accession, product in sorted(products.items()):
		key = " ".join((product or "").lower().split()) or accession
		groups.setdefault(key, []).append(accession)
	families = [sorted(members) for members in groups.values()]
	families.sort(key=lambda family: (-len(family), family[0]))
	return families
