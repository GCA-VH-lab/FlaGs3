import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from flags3 import home


class DomainError(RuntimeError):
	pass


@dataclass(frozen=True)
class DomainHit:
	protein: str
	name: str
	accession: str
	database: str
	group: str
	start: int
	end: int
	evalue: float


@dataclass(frozen=True)
class HmmSource:
	name: str
	path: Path
	query_cov: float = 0.0
	hmm_cov: float = 0.0
	group_sep: str = "__"

	@classmethod
	def parse(cls, spec: str, coverage: Optional[dict] = None) -> "HmmSource":
		name, _, path = spec.partition("=")
		if not path:
			name, path = "", spec
		path = Path(path).expanduser()
		if not name:
			name = path.name
			if name.lower().endswith(".hmm"):
				name = name[:-4]
		q, h = (coverage or {}).get(name, (coverage or {}).get("", (0.0, 0.0)))
		return cls(name, path, q, h)

	def files(self) -> list[Path]:
		if self.path.is_dir():
			found = sorted(p for p in self.path.iterdir() if p.suffix.lower() == ".hmm" and not p.name.startswith("."))
			if not found:
				raise DomainError("no .hmm files in {}".format(self.path))
			return found
		if not self.path.is_file():
			raise DomainError("HMM database not found: {} (flags3 install pfam sets up the default)".format(self.path))
		return [self.path]

	def group(self, name: str) -> str:
		return name.split(self.group_sep, 1)[0] if self.group_sep in name else name


def parse_coverage(specs) -> dict[str, tuple[float, float]]:
	out = {}
	for spec in specs or []:
		name, _, values = spec.rpartition("=")
		parts = [p for p in values.split(",") if p.strip()]
		try:
			nums = [float(p) for p in parts]
		except ValueError:
			raise DomainError("--hmm_coverage expects numbers, got {!r}".format(spec))
		if not nums or len(nums) > 2 or any(not 0 <= n <= 1 for n in nums):
			raise DomainError("--hmm_coverage takes one or two fractions between 0 and 1, got {!r}".format(spec))
		out[name] = (nums[0], nums[1] if len(nums) > 1 else 0.0)
	return out


class DomainScanner:
	def __init__(self, sources: list[HmmSource], evalue: float = 1e-3, cpus: int = 0):
		self.sources = sources
		self.evalue = evalue
		self.cpus = cpus
		self.counts: dict[str, int] = {}

	def scan(self, sequences: dict[str, str]) -> list[DomainHit]:
		from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
		if not sequences:
			return []
		alphabet = Alphabet.amino()
		block = DigitalSequenceBlock(alphabet, [
			TextSequence(name=name.encode(), sequence=seq).digitize(alphabet) for name, seq in sequences.items()])
		hits: list[DomainHit] = []
		for source in self.sources:
			before = len(hits)
			for path in source.files():
				for top_hits in self._search(path, block):
					for top in top_hits:
						hits.extend(self._collect(top, source, sequences))
			self.counts[source.name] = len(hits) - before
		return hits

	def _search(self, path: Path, block):
		import pyhmmer
		from pyhmmer.plan7 import HMMFile
		with HMMFile(str(path)) as probe:
			first = next(iter(probe), None)
			gathering = first is not None and first.cutoffs.gathering_available()
		handle = HMMFile(str(path))
		models = handle
		if handle.is_pressed:
			try:
				models = handle.optimized_profiles()
			except ValueError:
				handle.rewind()
		if gathering:
			yield pyhmmer.hmmer.hmmsearch(models, block, bit_cutoffs="gathering", cpus=self.cpus)
		else:
			yield pyhmmer.hmmer.hmmsearch(models, block, E=self.evalue, cpus=self.cpus)
		handle.close()

	def _collect(self, top, source: HmmSource, sequences) -> list[DomainHit]:
		name = _text(top.query.name)
		accession = _text(top.query.accession or b"")
		out = []
		for hit in top:
			protein = _text(hit.name)
			if protein not in sequences:
				continue
			for dom in hit.domains:
				if not dom.included:
					continue
				al = dom.alignment
				if not self._covered(al, source):
					continue
				out.append(DomainHit(protein, name, accession, source.name, source.group(name),
					al.target_from, al.target_to, dom.i_evalue))
		return out

	@staticmethod
	def _covered(al, source: HmmSource) -> bool:
		if source.query_cov > 0 and al.target_length:
			if (al.target_to - al.target_from) / al.target_length < source.query_cov:
				return False
		if source.hmm_cov > 0 and al.hmm_length:
			if (al.hmm_to - al.hmm_from) / al.hmm_length < source.hmm_cov:
				return False
		return True


def load_clans(path: Path) -> dict[str, str]:
	mapping: dict[str, str] = {}
	opener = gzip.open if str(path).endswith(".gz") else open
	with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			cols = line.rstrip("\n").split("\t")
			if len(cols) < 4 or not cols[1]:
				continue
			clan = cols[2] or cols[1]
			if cols[0]:
				mapping[cols[0]] = clan
			if cols[3]:
				mapping[cols[3]] = clan
	return mapping


class InterPro:
	COLUMNS = ("accession", "name", "type", "characterization_status", "informativeness_label", "interpretation")

	def __init__(self):
		self.by_pfam: dict[str, dict[str, str]] = {}
		self.collisions = 0

	@staticmethod
	def base(accession: str) -> str:
		return accession.split(".")[0].strip().upper()

	def load(self, path: Path) -> int:
		csv.field_size_limit(1 << 24)
		opener = gzip.open if str(path).endswith(".gz") else open
		with opener(path, "rt", encoding="utf-8", newline="") as handle:
			reader = csv.DictReader(handle, delimiter="\t")
			missing = [c for c in ("accession", "pfam_members") if c not in (reader.fieldnames or [])]
			if missing:
				raise DomainError("{} is not an InterPro metadata table (missing: {})".format(path, ", ".join(missing)))
			for row in reader:
				members = (row.get("pfam_members") or "").strip()
				if not members:
					continue
				record = {c: (row.get(c) or "").strip() for c in self.COLUMNS}
				for member in re.split(r"[;,\s]+", members):
					key = self.base(member)
					if not key:
						continue
					if key in self.by_pfam:
						self.collisions += 1
						continue
					self.by_pfam[key] = record
		return len(self.by_pfam)

	def get(self, accession: str) -> dict[str, str]:
		return self.by_pfam.get(self.base(accession), {})


def find_interpro(configured: Optional[Path]) -> Optional[Path]:
	candidates = [configured] if configured else [home.INTERPRO, home.INTERPRO.with_suffix(".tsv.gz")]
	for path in candidates:
		if path and path.is_file():
			return path
	if configured:
		raise DomainError("InterPro table not found: {}".format(configured))
	return None


def _text(value) -> str:
	return value.decode() if isinstance(value, (bytes, bytearray)) else value
