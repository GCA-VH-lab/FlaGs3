import gzip
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
	def __init__(self, hmm_db: str, evalue: float = 1e-10, cpus: int = 0):
		self.hmm_db = hmm_db
		self.evalue = evalue
		self.cpus = cpus  
		self.alphabet = Alphabet.amino()

	def scan(self, sequences: Dict[str, str]) -> Dict[str, List[DomainHit]]:
		if not sequences:
			return {}
		block = DigitalSequenceBlock(self.alphabet, [
			TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
			for name, seq in sequences.items()])
		with HMMFile(self.hmm_db) as handle:
			if handle.is_pressed:
				profiles = list(handle.optimized_profiles())
			else:
				profiles = list(handle)

		hits: Dict[str, List[DomainHit]] = {name: [] for name in sequences}
		for top in pyhmmer.hmmer.hmmsearch(profiles, block,
										   E=self.evalue, cpus=self.cpus):
			domain_name = self._decode(top.query.name)
			domain_acc = self._decode(top.query.accession or b"") 
			for hit in top:
				protein = self._decode(hit.name)
				for dom in hit.domains:
					if not dom.included:
						continue
					al = dom.alignment
					hits[protein].append(DomainHit(
						protein=protein, name=domain_name, accession=domain_acc,
						start=al.target_from, end=al.target_to, evalue=dom.i_evalue))
		return hits

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
		header = ["#protein", "family", "domain", "pfam", "clan",
				  "start", "end", "evalue"]
		if interpro:
			header += ["interpro", "interpro_name", "interpro_type",
					   "characterization", "informativeness", "interpretation"]
		with open(path, "w") as out:
			out.write("\t".join(header) + "\n")
			for protein in sorted(hits):
				for d in sorted(hits[protein], key=lambda x: x.start):
					row = [protein, families.get(protein, "-"), d.name,
						   d.accession or "-", clans.get(d.name, "-"),
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