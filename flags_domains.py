import gzip
import os
import re
from typing import Dict, List, NamedTuple

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from pyhmmer.plan7 import HMMFile


class DomainHit(NamedTuple):
	protein: str
	name: str 
	start: int 
	end: int
	evalue: float
	accession: str = ""
	database: str = ""
	group: str = ""


class HmmSource(NamedTuple):
	name: str
	path: str
	query_cov: float = 0.0
	hmm_cov: float = 0.0
	group_sep: str = "__"

	@classmethod
	def parse(cls, spec: str, coverage=None, group_sep: str = "__"):
		name, _, path = spec.partition("=")
		if not path:
			name, path = "", spec
		path = os.path.expanduser(path)
		if not name:
			name = os.path.basename(os.path.normpath(path))
			for ext in (".hmm", ".HMM"):
				if name.endswith(ext):
					name = name[:-len(ext)]
		q, h = (coverage or {}).get(name, (coverage or {}).get("", (0.0, 0.0)))
		return cls(name=name, path=path, query_cov=q, hmm_cov=h,
				   group_sep=group_sep)

	def models(self):
		if os.path.isdir(self.path):
			files = sorted(os.path.join(self.path, f)
						   for f in os.listdir(self.path)
						   if f.lower().endswith(".hmm") and not f.startswith("."))
			if not files:
				raise ValueError("no .hmm files in {}".format(self.path))
		else:
			files = [self.path]
		models = []
		for path in files:
			with HMMFile(path) as handle:
				if len(files) == 1 and handle.is_pressed:
					try:
						models.extend(handle.optimized_profiles())
						continue
					except ValueError:
						handle.rewind()
				models.extend(handle)
		return models

	def label(self, name: str) -> str:
		if self.group_sep and self.group_sep in name:
			return name.split(self.group_sep, 1)[0]
		return name


class InterProAnnotator:

	COLUMNS = ("accession", "name", "type",
			   "characterization_status", "informativeness_label", "interpretation")

	def __init__(self):
		self.by_pfam: Dict[str, Dict[str, str]] = {}
		self.collisions = 0

	@staticmethod
	def _base(accession: str) -> str:
		return accession.split(".")[0].strip().upper()

	def load(self, path: str) -> int:
		import csv
		csv.field_size_limit(1 << 24)   # description fields are far past the default
		opener = gzip.open if path.endswith(".gz") else open
		with opener(path, "rt", encoding="utf-8", newline="") as fh:
			reader = csv.DictReader(fh, delimiter="\t")
			missing = [c for c in ("accession", "pfam_members") if c not in
					   (reader.fieldnames or [])]
			if missing:
				raise ValueError(
					"{} is not an InterPro metadata table (missing column{}: {})".format(
						path, "s" if len(missing) > 1 else "", ", ".join(missing)))
			for row in reader:
				members = (row.get("pfam_members") or "").strip()
				if not members:
					continue
				record = {c: (row.get(c) or "").strip() for c in self.COLUMNS}
				for member in re.split(r"[;,\s]+", members):
					key = self._base(member)
					if not key:
						continue
					if key in self.by_pfam:
						self.collisions += 1
						continue   # first entry wins, keeping the join deterministic
					self.by_pfam[key] = record
		return len(self.by_pfam)

	def get(self, accession: str, name: str = "") -> Dict[str, str]:
		return self.by_pfam.get(self._base(accession), {})


class DomainScanner: 
	def __init__(self, sources, evalue: float = 1e-10, cpus: int = 0):
		if isinstance(sources, str):
			sources = [HmmSource.parse(sources)]
		self.sources = list(sources)
		self.evalue = evalue
		self.cpus = cpus  
		self.alphabet = Alphabet.amino()
		self.counts: Dict[str, int] = {}

	def scan(self, sequences: Dict[str, str]) -> Dict[str, List[DomainHit]]:
		if not sequences:
			return {}
		block = DigitalSequenceBlock(self.alphabet, [
			TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
			for name, seq in sequences.items()])
		hits: Dict[str, List[DomainHit]] = {name: [] for name in sequences}
		for source in self.sources:
			found = 0
			for model in self._searches(source, block):
				for top in model:
					found += self._collect(top, source, hits)
			self.counts[source.name] = found
		return hits

	def _searches(self, source: HmmSource, block):
		models = source.models()
		gated = [m for m in models if getattr(m, "cutoffs", None) is not None
				 and m.cutoffs.gathering_available()]
		plain = [m for m in models if m not in gated] if gated else models
		runs = []
		if gated:
			runs.append(pyhmmer.hmmer.hmmsearch(
				gated, block, bit_cutoffs="gathering", cpus=self.cpus))
		if plain:
			runs.append(pyhmmer.hmmer.hmmsearch(
				plain, block, E=self.evalue, cpus=self.cpus))
		return runs

	def _collect(self, top, source: HmmSource, hits) -> int:
		name = self._decode(top.query.name)
		accession = self._decode(top.query.accession or b"")
		group = source.label(name)
		found = 0
		for hit in top:
			protein = self._decode(hit.name)
			if protein not in hits:
				continue
			for dom in hit.domains:
				if not dom.included:
					continue
				al = dom.alignment
				if not self._covered(al, source):
					continue
				hits[protein].append(DomainHit(
					protein=protein, name=name, accession=accession,
					start=al.target_from, end=al.target_to, evalue=dom.i_evalue,
					database=source.name, group=group))
				found += 1
		return found

	@staticmethod
	def _covered(al, source: HmmSource) -> bool:
		if source.query_cov > 0 and al.target_length:
			if (al.target_to - al.target_from) / al.target_length < source.query_cov:
				return False
		if source.hmm_cov > 0 and al.hmm_length:
			if (al.hmm_to - al.hmm_from) / al.hmm_length < source.hmm_cov:
				return False
		return True

	@staticmethod
	def _decode(value) -> str:
		return value.decode() if isinstance(value, (bytes, bytearray)) else value

	@staticmethod
	def write_report(hits: Dict[str, List[DomainHit]], path: str,
					 clans: Dict[str, str] = None,
					 families: Dict[str, str] = None,
					 interpro: "InterProAnnotator" = None):
		clans = clans or {}
		families = families or {}
		header = ["#protein", "family", "database", "domain", "group", "pfam",
				  "clan", "start", "end", "evalue"]
		if interpro:
			header += ["interpro", "interpro_name", "interpro_type",
					   "characterization", "informativeness", "interpretation"]
		with open(path, "w") as out:
			out.write("\t".join(header) + "\n")
			for protein in sorted(hits):
				for d in sorted(hits[protein], key=lambda x: x.start):
					row = [protein, families.get(protein, "-"), d.database or "-",
						   d.name, d.group or "-", d.accession or "-",
						   clans.get(d.name, "-"),
						   str(d.start), str(d.end), "{:.2e}".format(d.evalue)]
					if interpro:
						meta = interpro.get(d.accession, d.name)
						row += [meta.get("accession", "-") or "-",
								meta.get("name", "-") or "-",
								meta.get("type", "-") or "-",
								meta.get("characterization_status", "-") or "-",
								meta.get("informativeness_label", "-") or "-",
								meta.get("interpretation", "-") or "-"]
					out.write("\t".join(row) + "\n")

	@staticmethod
	def load_clans(path: str) -> Dict[str, str]:
		mapping: Dict[str, str] = {}
		opener = gzip.open if path.endswith(".gz") else open
		with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
			for line in fh:
				cols = line.rstrip("\n").split("\t")
				if len(cols) < 4:
					continue
				pfam_id, clan_id, clan_name, family_name = cols[0], cols[1], cols[2], cols[3]
				if not clan_id:
					continue 
				clan = clan_name or clan_id
				if pfam_id:
					mapping[pfam_id] = clan
				if family_name:
					mapping[family_name] = clan
		return mapping