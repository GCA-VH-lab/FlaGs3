from flags_report import family_numbers
from flags_log import debug

DEFAULT_INTERPRO = "interpro_metadata_processed.tsv"

import gzip
import os
import sys
import time
from collections import Counter
from typing import Optional
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

	def single_file(self):
		return None if os.path.isdir(self.path) else self.path

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
		csv.field_size_limit(1 << 24)
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
						continue
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
		self._open: List = []

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
		for handle in self._open:
			handle.close()
		self._open = []
		return hits

	def _searches(self, source: HmmSource, block):
		single = source.single_file()
		if single:
			return [self._stream_search(single, source, block)]
		gated, plain = [], []
		for model in source.models():
			cutoffs = getattr(model, "cutoffs", None)
			if cutoffs is not None and cutoffs.gathering_available():
				gated.append(model)
			else:
				plain.append(model)
		runs = []
		if gated:
			runs.append(pyhmmer.hmmer.hmmsearch(
				gated, block, bit_cutoffs="gathering", cpus=self.cpus))
		if plain:
			runs.append(pyhmmer.hmmer.hmmsearch(
				plain, block, E=self.evalue, cpus=self.cpus))
		return runs

	def _stream_search(self, path: str, source: HmmSource, block):
		with HMMFile(path) as handle:
			first = next(iter(handle), None)
			cutoffs = getattr(first, "cutoffs", None)
			gathering = cutoffs is not None and cutoffs.gathering_available()
		handle = HMMFile(path)
		self._open.append(handle)
		models = handle
		if handle.is_pressed:
			try:
				models = handle.optimized_profiles()
			except ValueError:
				handle.rewind()
		if gathering:
			return pyhmmer.hmmer.hmmsearch(
				models, block, bit_cutoffs="gathering", cpus=self.cpus)
		return pyhmmer.hmmer.hmmsearch(models, block, E=self.evalue,
									   cpus=self.cpus)

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


def parse_coverage(specs):
	out = {}
	for spec in specs or []:
		name, _, values = spec.rpartition("=")
		parts = [p for p in values.split(",") if p.strip()]
		try:
			nums = [float(p) for p in parts]
		except ValueError:
			raise ValueError("--hmm_coverage expects numbers, got {!r}".format(spec))
		if not nums or len(nums) > 2 or any(not 0 <= n <= 1 for n in nums):
			raise ValueError(
				"--hmm_coverage takes one or two fractions between 0 and 1, "
				"got {!r}".format(spec))
		out[name] = (nums[0], nums[1] if len(nums) > 1 else 0.0)
	return out


def resolve_interpro(path: str) -> Optional[str]:
	candidates = [path]
	if not os.path.isabs(path):
		candidates.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
									   path))
		candidates.append(path + ".gz")
		candidates.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
									   path + ".gz"))
	for candidate in candidates:
		if os.path.isfile(candidate):
			return candidate
	if path != DEFAULT_INTERPRO:
		sys.exit("Error: --interpro file not found: {}".format(path))
	debug("no InterPro table at {}; domain table will omit those columns".format(path))
	return None


def scan_domains(args, extractor, families, all_neighborhoods, out_path,
				 timings, t0, want_domain_fig):
	domains, clan_map, domain_table_written = {}, None, False
	if want_domain_fig:
		domains = {}
		clan_map = None
		if args.domains and not args.hmmdb:
			print("Warning: --domains needs --hmmdb; drawing the figure without domains.")
		elif args.domains:
			try:
				import flags_domains as dom_mod
				if args.verbose:
					print(">> scanning {} proteins for domains...".format(
						len(extractor.sequences)), flush=True)
				sources = [dom_mod.HmmSource.parse(spec, args.hmm_coverage)
						   for spec in args.hmmdb]
				scanner = dom_mod.DomainScanner(sources, evalue=args.ethreshold,
												cpus=args.cpu or 0)
				domains = scanner.scan(extractor.sequences)
				clan_map = dom_mod.DomainScanner.load_clans(args.clans) if args.clans else None
				if args.verbose:
					for name, n in scanner.counts.items():
						print(">> {}: {} domain hits".format(name, n), flush=True)
			except Exception as e:
				print("Warning: could not read the HMM database, drawing the figure without domains ({}).".format(e))
		if domains:
			interpro = None
			interpro_path = resolve_interpro(args.interpro) if args.interpro else None
			if interpro_path:
				try:
					interpro = dom_mod.InterProAnnotator()
					n = interpro.load(interpro_path)
					if args.verbose:
						print(">> InterPro: {} Pfam entries mapped from {}{}".format(
							n, os.path.basename(interpro_path),
							", {} duplicate Pfam ids ignored".format(interpro.collisions)
							if interpro.collisions else ""), flush=True)
				except Exception as e:
					print("Warning: could not read the InterPro table, writing the "
						  "domain table without it ({}).".format(e))
					interpro = None
			try:
				dom_mod.DomainScanner.write_report(
					domains, out_path("_domains.tsv"), clans=clan_map,
					interpro=interpro,
					families=family_numbers(
						families,
						{g.accession for g in all_neighborhoods if g.is_rna},
						{g.accession for g in all_neighborhoods if g.offset == 0},
						Counter(g.accession for g in all_neighborhoods)))
				domain_table_written = True
			except Exception as e:
				print("Warning: could not write the domain table ({}).".format(e))
		timings["8_domains"] = time.perf_counter() - t0
	return domains, clan_map, domain_table_written
