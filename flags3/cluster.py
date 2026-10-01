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


LETTERS = "abcdefghijklmnopqrstuvwxyz"
UNIFORM = 0.95


def subfamilies(family: list[str], full: dict[str, set], contains: dict[str, set]) -> tuple[dict[str, str], list[list[str]]]:
	members = set(family)
	inner = {m: (full.get(m, set()) & members) for m in members}
	groups = [g for g in connected_components(inner) if len(g) > 1]
	if len(groups) < 2:
		return {}, []
	holders = {m: {h for h in members if m in contains.get(h, ())} for m in members}

	def held_by(m, group) -> bool:
		return len(holders[m] & set(group)) >= max(1, len(group) / 2)

	merged = True
	while merged and len(groups) > 1:
		merged = False
		for i, g in enumerate(groups):
			targets = [{j for j, h in enumerate(groups) if j != i and held_by(m, h)} for m in g]
			common = set.intersection(*targets) if targets else set()
			if len(common) == 1:
				groups[common.pop()].extend(g)
				del groups[i]
				merged = True
				break
	if len(groups) < 2:
		return {}, []
	groups.sort(key=lambda g: (-len(g), sorted(g)[0]))
	groups = [sorted(g) for g in groups]
	group_of = {m: i for i, g in enumerate(groups) for m in g}
	holds_any = {m: {group_of[h] for h in contains.get(m, ()) if h in group_of and h in members} - ({group_of[m]} if m in group_of else set())
		for m in members}
	uniform = {i: {j for j in range(len(groups)) if j != i
		and sum(1 for m in g if j in holds_any[m]) >= UNIFORM * len(g)} for i, g in enumerate(groups)}
	letter = lambda i: LETTERS[i % 26] * (i // 26 + 1)
	labels = {}
	for m in members:
		own = {group_of[m]} if m in group_of else set()
		holds = {j for j in holds_any[m] if not any(j in uniform[i] for i in own)}
		inside = set() if own else {j for j, g in enumerate(groups) if held_by(m, g)}
		letters = sorted(own | holds | inside)
		labels[m] = "".join(letter(i) for i in letters) if letters else "?"
	return labels, groups


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
		self.coverage: dict[tuple[str, str], tuple[float, float]] = {}

	def cluster(self, sequences: dict[str, str]) -> list[list[str]]:
		self.adjacency = {}
		self.coverage = {}
		if not sequences:
			return []
		self.adjacency = self.edges(sequences)
		for name in sequences:
			self.adjacency.setdefault(name, set())
		return connected_components(self.adjacency)

	def full_length(self, threshold: float) -> dict[str, set]:
		full: dict[str, set] = {name: set() for name in self.adjacency}
		if threshold <= 0:
			return full
		for (a, b), (cov_a, cov_b) in self.coverage.items():
			if a != b and cov_a >= threshold and cov_b >= threshold:
				full.setdefault(a, set()).add(b)
				full.setdefault(b, set()).add(a)
		return full

	def contains(self, threshold: float) -> dict[str, set]:
		out: dict[str, set] = {name: set() for name in self.adjacency}
		if threshold <= 0:
			return out
		for (a, b), (cov_a, cov_b) in self.coverage.items():
			if a == b:
				continue
			if cov_b >= threshold:
				out.setdefault(a, set()).add(b)
			if cov_a >= threshold:
				out.setdefault(b, set()).add(a)
		return out

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
			out = []
			for (name, query), result in zip(chunk, results):
				hits, covs = set(), {}
				for h in self.hits_of(result):
					if not h.included:
						continue
					target = _name(h)
					hits.add(target)
					covs[(name, target)] = _mutual_coverage(h, len(query))
				out.append((name, hits, covs))
			return out

		adjacency: dict[str, set] = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			for part in pool.map(run, chunks):
				for name, hits, covs in part:
					adjacency[name] = hits
					self.coverage.update(covs)
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
			incE=option(self.tool, "incE", 1e-10), incdomE=option(self.tool, "incdomE", 1e-10), cpus=1)


class NhmmerClusterer(PyhmmerClusterer):
	alphabet_name = "dna"

	def search(self, queries, block):
		import pyhmmer
		return pyhmmer.hmmer.nhmmer(queries, block, incE=option(self.tool, "incE", 1e-10), cpus=1)


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

	def read_pairs(self, hits: Path, sequences: dict[str, str]) -> dict[str, set]:
		adjacency: dict[str, set] = {name: set() for name in sequences}
		with open(hits, encoding="utf-8") as handle:
			for line in handle:
				cells = line.rstrip("\n").split("\t")
				if len(cells) >= 2 and cells[0] in adjacency and cells[1] in adjacency:
					adjacency[cells[0]].add(cells[1])
					if len(cells) >= 4:
						try:
							self.coverage[(cells[0], cells[1])] = (float(cells[2]), float(cells[3]))
						except ValueError:
							pass
		return adjacency


ENGINE_CLASSES = {"jackhmmer": JackhmmerClusterer, "nhmmer": NhmmerClusterer, "mmseqs": MmseqsClusterer}


def build(tools: Tools, name: str, workers: Optional[int], work: Path) -> Clusterer:
	tool = method(tools, name)
	return ENGINE_CLASSES[tool.engine](tool, workers, work)


def _mutual_coverage(hit, query_length: int) -> tuple[float, float]:
	domains = [d for d in getattr(hit, "domains", []) if getattr(d, "included", True)]
	if not domains:
		return 0.0, 0.0
	t_lo = min(d.alignment.target_from for d in domains)
	t_hi = max(d.alignment.target_to for d in domains)
	q_lo = min(d.alignment.hmm_from for d in domains)
	q_hi = max(d.alignment.hmm_to for d in domains)
	target_length = domains[0].alignment.target_length or 1
	query_length = domains[0].alignment.hmm_length or query_length or 1
	return (q_hi - q_lo + 1) / query_length, (t_hi - t_lo + 1) / target_length


def _name(hit) -> str:
	raw = hit.name
	return raw.decode() if isinstance(raw, bytes) else raw
