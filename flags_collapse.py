import os
import shutil
import subprocess
import tempfile
from typing import Dict, List, NamedTuple, Optional, Tuple

QUERY_SECONDS = 0.1365      # measured jackhmmer cost, fixed part
BLOCK_SECONDS = 4.365e-05   # measured jackhmmer cost, per sequence in the block


class CollapseResult(NamedTuple):
	representatives: Dict[str, str]
	members: Dict[str, List[str]]
	min_seq_id: float
	coverage: float
	command: str


def estimate_seconds(n: int) -> float:
	return n * (QUERY_SECONDS + BLOCK_SECONDS * n)


def format_estimate(n: int, workers: Optional[int] = None) -> str:
	seconds = estimate_seconds(n)
	if workers and workers > 1:
		seconds /= workers
	if seconds < 90:
		return "{:.0f} s".format(seconds)
	if seconds < 5400:
		return "{:.0f} min".format(seconds / 60)
	if seconds < 86400 * 2:
		return "{:.1f} h".format(seconds / 3600)
	return "{:.1f} days".format(seconds / 86400)


class Collapser:
	def __init__(self, min_seq_id: float = 0.9, coverage: float = 0.8,
				 threads: int = 0, keep_dir: Optional[str] = None):
		self.min_seq_id = min_seq_id
		self.coverage = coverage
		self.threads = threads
		self.keep_dir = keep_dir
		self.command = ""

	def available(self) -> Tuple[bool, str]:
		import flags_tools
		cmd, wd = flags_tools.command(
			"mmseqs", **{"in": "x", "out": "y", "tmp": "z",
						 "id": self.min_seq_id, "cov": self.coverage,
						 "threads": self.threads or 1})
		if not cmd:
			return False, "no mmseqs row in tools_table.tsv"
		if wd:
			candidate = os.path.join(wd, str(cmd[0]))
			if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
				return True, candidate
		return flags_tools.locate(cmd)

	def collapse(self, sequences: Dict[str, str]) -> CollapseResult:
		if not sequences:
			return CollapseResult({}, {}, self.min_seq_id, self.coverage, "")
		work = self.keep_dir or tempfile.mkdtemp(prefix="flags3_collapse_")
		os.makedirs(work, exist_ok=True)
		fasta = os.path.join(work, "input.fasta")
		with open(fasta, "w") as out:
			for name, seq in sequences.items():
				out.write(">{}\n{}\n".format(name, seq))

		import flags_tools
		prefix = os.path.join(work, "clu")
		cmd, wd = flags_tools.command(
			"mmseqs", **{"in": fasta, "out": prefix,
						 "tmp": os.path.join(work, "tmp"),
						 "id": self.min_seq_id, "cov": self.coverage,
						 "threads": self.threads or 1})
		self.command = " ".join(cmd)
		try:
			proc = subprocess.run(cmd, cwd=wd or None, capture_output=True, text=True,
								  env=flags_tools.env_for(cmd))
		except FileNotFoundError:
			raise FileNotFoundError(
				"{} not found; run mmseqs_installer.sh or drop --cluster_collapse."
				.format(cmd[0] if cmd else "mmseqs"))
		finally:
			self._record(cmd, locals().get("proc"))
		if proc.returncode != 0:
			raise RuntimeError("mmseqs exited {}: {}".format(
				proc.returncode, flags_tools.brief(proc.stderr or proc.stdout)))

		members = self._read_clusters(prefix + "_cluster.tsv", sequences)
		representatives = {rep: sequences[rep] for rep in members}
		if not self.keep_dir:
			shutil.rmtree(work, ignore_errors=True)
		return CollapseResult(representatives, members,
							  self.min_seq_id, self.coverage, self.command)

	@staticmethod
	def _record(cmd, proc):
		try:
			import flags_log
		except ImportError:
			return
		if proc is None:
			flags_log.record_command(cmd, "not run")
		else:
			flags_log.record_command(cmd, proc.returncode, proc.stdout, proc.stderr)

	@staticmethod
	def _read_clusters(path: str, sequences: Dict[str, str]) -> Dict[str, List[str]]:
		if not os.path.isfile(path):
			raise RuntimeError("mmseqs produced no cluster table at {}".format(path))
		members: Dict[str, List[str]] = {}
		seen = set()
		with open(path) as fh:
			for line in fh:
				parts = line.rstrip("\n").split("\t")
				if len(parts) < 2:
					continue
				rep, member = parts[0], parts[1]
				if rep not in sequences or member not in sequences:
					continue
				members.setdefault(rep, []).append(member)
				seen.add(member)
		for name in sequences:
			if name not in seen:
				members.setdefault(name, []).append(name)
		return members


def expand(families: List[List[str]], members: Dict[str, List[str]]) -> List[List[str]]:
	out = []
	for family in families:
		grown = set()
		for rep in family:
			grown.update(members.get(rep, [rep]))
		out.append(sorted(grown))
	out.sort(key=len, reverse=True)
	return out


def write_report(result: CollapseResult, path: str, fam_of=None):
	fam_of = fam_of or {}
	with open(path, "w") as out:
		out.write("# flanking proteins collapsed before clustering\n")
		out.write("# min_seq_id={} coverage={}\n".format(
			result.min_seq_id, result.coverage))
		out.write("# {}\n".format(result.command))
		out.write("# a member shares its representative's family; the representative's\n"
				  "# own evidence is in _jackhits.tsv\n")
		out.write("#representative\tfamily\tn_members\tmembers\n")
		for rep in sorted(result.members, key=lambda r: (-len(result.members[r]), r)):
			group = sorted(result.members[rep])
			out.write("{}\t{}\t{}\t{}\n".format(
				rep, fam_of.get(rep, "-"), len(group), ",".join(group)))
