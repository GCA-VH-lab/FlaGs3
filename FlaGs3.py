import argparse
import atexit
import gzip
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter
from concurrent.futures import (ThreadPoolExecutor, as_completed, wait,
								FIRST_COMPLETED)
from typing import Dict, List, Optional, Tuple, NamedTuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

import flags_log
from flags_log import debug, set_debug

NCBI_TOOL = "flags3"


DEFAULT_HMMDB = "./pfam_db/Pfam-A.hmm"


def parse_collapse(spec):
	parts = [p for p in str(spec).split(",") if p.strip()]
	try:
		nums = [float(p) for p in parts]
	except ValueError:
		raise ValueError("--cluster_collapse expects fractions, got {!r}".format(spec))
	if not nums or len(nums) > 2 or any(not 0 < n <= 1 for n in nums):
		raise ValueError(
			"--cluster_collapse takes an identity and optionally a coverage, "
			"each above 0 and at most 1, got {!r}".format(spec))
	return (nums[0], nums[1] if len(nums) > 1 else 0.8)


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
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock


def plural(n, word, plural_form=None):
	return "{} {}".format(n, word if n == 1 else (plural_form or word + "s"))


class AccessionListReader: 
	def __init__(self, paths):
		self.paths = [p for p in (paths or [])]

	BLAST_MARK = "blast"

	def read(self):
		proteins_assembly, proteins_only, blast_queries = [], [], []
		for path in self.paths:
			with open(path, 'r') as f:
				for n, line in enumerate(f, 1):
					line = line.split('#', 1)[0].strip()
					if not line:
						continue
					if '\t' in line:
						fields = [c for c in line.split('\t') if c]
						if len(fields) >= 2:
							if fields[1].strip().lower() == self.BLAST_MARK:
								blast_queries.append((fields[0], path, n))
							else:
								proteins_assembly.append([fields[0], fields[1]])
						else:
							print("Warning: line {} of {} is malformed, skipping it: "
								  "{!r}".format(n, path, line))
					else:
						proteins_only.append(line)
		seen = set()
		proteins_only = [p for p in proteins_only
						 if not (p in seen or seen.add(p))]
		return proteins_assembly, proteins_only, blast_queries


class ProteinAssemblyMapper: 
	REFSEQ_PROTEIN = re.compile(r"^[A-Z]{2}_")

	def __init__(self, email: str, max_assemblies: int = 1,
				 api_key: Optional[str] = None, ncbi_time: float = 0.4,
				 cross_db: bool = True):
		self.max_assemblies = max_assemblies
		self.ncbi_time = ncbi_time
		self.cross_db = cross_db
		self.accessions_in: Dict[str, Dict[str, set]] = {}
		self.dropped_cross_db: Dict[str, set] = {}
		self.unreachable = ""
		Entrez.email = email
		Entrez.tool = NCBI_TOOL
		if api_key:
			Entrez.api_key = api_key

	@classmethod
	def _same_db(cls, protein: str, assembly: str) -> bool:
		"""RefSeq proteins (WP_, NP_, ...) belong with GCF_; INSDC proteins with GCA_."""
		want = "GCF" if cls.REFSEQ_PROTEIN.match(protein) else "GCA"
		return assembly[:3] == want

	def map(self, proteins: List[str]) -> Dict[str, List[str]]:
		if not proteins:
			return {}
		xp = [p for p in proteins if p.startswith("XP_")]
		ipg = [p for p in proteins if not p.startswith("XP_")]
		found: Dict[str, set] = {}
		if ipg:
			try:
				found.update(self._map_ipg(ipg))
			except Exception as e:
				self.unreachable = "{}: {}".format(type(e).__name__, e)
				debug("IPG lookup failed for {} proteins".format(len(ipg)), exc=True)
		for acc in xp:
			try:
				found[acc] = self._map_xp(acc)
			except Exception as e:
				self.unreachable = "{}: {}".format(type(e).__name__, e)
				debug("BioProject lookup failed for {}".format(acc), exc=True)
		if not self.cross_db:
			for acc, asm in found.items():
				keep = {a for a in asm if self._same_db(acc, a)}
				if asm - keep:
					self.dropped_cross_db[acc] = asm - keep
				found[acc] = keep
		gcf_first = lambda a: (0 if a[:3] == "GCF" else 1, a)
		return {acc: sorted(asm, key=gcf_first)[:self.max_assemblies]
				for acc, asm in found.items()}

	IPG_CHUNK = 200        # NCBI's recommended ceiling for one efetch

	def _map_ipg(self, proteins: List[str]) -> Dict[str, set]:
		"""Resolve proteins to assemblies in chunks. One request for thousands of
		accessions returns a report large enough to time out or come back
		truncated, and a truncated IPG report is indistinguishable from those
		proteins having no assembly -- a silent wrong answer rather than an
		error. A failed chunk is reported and the rest still resolve."""
		found: Dict[str, set] = {acc: set() for acc in proteins}
		failed = 0
		for start in range(0, len(proteins), self.IPG_CHUNK):
			chunk = proteins[start:start + self.IPG_CHUNK]
			try:
				found.update(self._map_ipg_chunk(chunk))
			except Exception as e:
				failed += len(chunk)
				self.unreachable = "{}: {}".format(type(e).__name__, e)
				debug("IPG chunk {}-{} failed".format(
					start, start + len(chunk)), exc=True)
		if failed:
			print("Warning: {} of {} could not be resolved through IPG; they are "
				  "reported as unresolved in <prefix>_accessionIssues.txt.".format(
					  plural(failed, "protein"), len(proteins)))
		return found

	def _map_ipg_chunk(self, proteins: List[str]) -> Dict[str, set]:
		queries = set(proteins)
		time.sleep(self.ncbi_time)
		handle = Entrez.efetch(db="ipg", id=",".join(proteins),
							   rettype="ipg", retmode="text")
		data = handle.read()
		handle.close()
		if isinstance(data, bytes):
			data = data.decode("utf-8", errors="replace")
		groups: Dict[str, dict] = {}
		for line in data.splitlines():
			if line[0:2] == "Id" or not re.search(r"GC._\d*\.\d", line):
				continue
			fields = line.rstrip().split("\t")
			ipg_id, acc, assembly = fields[0], fields[6], fields[-1]
			grp = groups.setdefault(ipg_id, {"queries": set(), "by_asm": {}})
			grp["by_asm"].setdefault(assembly, set()).add(acc)
			if acc in queries:
				grp["queries"].add(acc)

		found = {acc: set() for acc in proteins}
		for grp in groups.values():
			for q in grp["queries"]:
				self.accessions_in.setdefault(q, {})
				for assembly, accs in grp["by_asm"].items():
					found[q].add(assembly)
					self.accessions_in[q].setdefault(assembly, set()).update(accs)
		return found

	def _map_xp(self, accession: str) -> set:
		try:
			time.sleep(self.ncbi_time)
			handle = Entrez.efetch(db="protein", id=accession,
								   rettype="gbwithparts", retmode="text")
			record = SeqIO.read(handle, "genbank")
			handle.close()
		except Exception:
			return set()
		bioprojects = [x.split(":", 1)[1] for x in record.dbxrefs
					   if x.split(":", 1)[0] == "BioProject"]
		assemblies: set = set()
		for bp in bioprojects:
			assemblies.update(self._assemblies_for_bioproject(bp))
		self.accessions_in.setdefault(accession, {})
		for asm in assemblies:
			self.accessions_in[accession].setdefault(asm, set()).add(accession)
		return assemblies

	def _assemblies_for_bioproject(self, bioproject: str) -> set:
		try:
			time.sleep(self.ncbi_time)
			search = Entrez.read(Entrez.esearch(db="bioproject", term=bioproject))
			ids = search.get("IdList", [])
			if not ids:
				return set()
			time.sleep(self.ncbi_time)
			links = Entrez.read(Entrez.elink(dbfrom="bioproject", db="assembly",
											 id=",".join(ids)))
			asm_ids = []
			for linkset in links:
				for db in linkset.get("LinkSetDb", []):
					asm_ids.extend(link["Id"] for link in db.get("Link", []))
			if not asm_ids:
				return set()
			time.sleep(self.ncbi_time)
			summary = Entrez.read(Entrez.esummary(db="assembly",
												  id=",".join(asm_ids)))
			docs = summary["DocumentSummarySet"]["DocumentSummary"]
			out = set()
			for d in docs:
				acc = d.get("AssemblyAccession", "")
				if re.match(r"GC._\d+\.\d", acc):
					out.add(acc)
			return out
		except Exception:
			return set()


class RateLimiter: 

	def __init__(self, rate: float):
		self.min_interval = 1.0 / rate if rate > 0 else 0.0
		self._lock = threading.Lock()
		self._next_slot = time.monotonic()

	def wait(self):
		if self.min_interval <= 0:
			return
		with self._lock:
			now = time.monotonic()
			start = max(now, self._next_slot)
			self._next_slot = start + self.min_interval
		delay = start - now
		if delay > 0:
			time.sleep(delay)


class GenomeFiles(NamedTuple):
	gff: Optional[str] = None
	faa: Optional[str] = None
	rna: Optional[str] = None       
	genome: Optional[str] = None    


class _GenomeDownloader: 
	SUFFIX: Dict[str, str] = {}

	def __init__(self, out_dir: Optional[str] = None, workers: int = 8,
				 rate: float = 5.0, want_rna: bool = False, want_genome: bool = False):
		self.out_dir = out_dir or tempfile.gettempdir()
		self.workers = workers
		self.want_rna = want_rna
		self.want_genome = want_genome
		self.limiter = RateLimiter(rate)
		self.failures: Dict[str, str] = {}  
		os.makedirs(self.out_dir, exist_ok=True)
		self.session = requests.Session()
		self.session.headers["User-Agent"] = "{}/{}".format(NCBI_TOOL, VERSION)
		retry = Retry(total=5, backoff_factor=0.5, respect_retry_after_header=True,
					  status_forcelist=[429, 500, 502, 503, 504],
					  allowed_methods=frozenset(["GET"]))
		size = workers * 4
		self.session.mount("https://", HTTPAdapter(
			max_retries=retry, pool_connections=size, pool_maxsize=size))

	def download_many(self, assemblies: List[str], progress_cb=None) -> Dict[str, GenomeFiles]:
		if not assemblies:
			return {}
		results = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			futures = {pool.submit(self._fetch_one, a): a for a in assemblies}
			for done, fut in enumerate(as_completed(futures), 1):
				results[futures[fut]] = fut.result()
				if progress_cb:
					progress_cb(done, len(assemblies))
		return results

	def _slots(self) -> List[str]:
		slots = ["gff", "faa"]
		if self.want_rna and "rna" in self.SUFFIX:
			slots.append("rna")
		if self.want_genome:
			slots.append("genome")
		return slots

	def _fetch_files(self, jobs: Dict[str, Tuple[Optional[str], str]]) -> GenomeFiles:

		with ThreadPoolExecutor(max_workers=max(len(jobs), 1)) as pool:
			done = {slot: pool.submit(self._stream, url, local)
					for slot, (url, local) in jobs.items()}
			return GenomeFiles(**{slot: jobs[slot][1]
								  for slot, fut in done.items() if fut.result()})

	def _stream(self, url: Optional[str], local: str) -> bool:
		if not url:
			return False
		for attempt in range(3):
			try:
				t_wait = time.monotonic()
				self.limiter.wait()
				t_send = time.monotonic()
				with self.session.get(url, stream=True, timeout=120) as r:
					t_head = time.monotonic()
					if r.status_code != 200:
						note = "HTTP {}".format(r.status_code)
						retry_after = r.headers.get("Retry-After")
						if retry_after:
							note += " Retry-After={}".format(retry_after)
						self.failures[url] = note
						debug("{} {}".format(note, url))
						return False
					size = 0
					with open(local, "wb") as fout:
						for chunk in r.iter_content(chunk_size=1 << 16):
							fout.write(chunk)
							size += len(chunk)
				t_done = time.monotonic()
				debug("got {:>9,d} B in {:5.2f}s (limiter {:4.2f}s, ttfb {:5.2f}s, "
					  "body {:5.2f}s, {:6.1f} kB/s) {}".format(
						  size, t_done - t_wait, t_send - t_wait, t_head - t_send,
						  t_done - t_head,
						  size / 1024.0 / max(t_done - t_head, 1e-6),
						  os.path.basename(local)))
				if size > 0:
					self.failures.pop(url, None)
					return True
				self.failures[url] = "empty response body"
				debug("empty body {}".format(url))
			except Exception as e:
				self.failures[url] = "{}: {}".format(type(e).__name__, e)
				debug("attempt {}/3 failed for {}: {!r}".format(attempt + 1, url, e),
					  exc=True)
				time.sleep(0.5)
		return False


class AssemblyDownloader(_GenomeDownloader): 
	BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
	SUFFIX = {"gff": "_genomic.gff.gz", "faa": "_protein.faa.gz",
			  "rna": "_rna_from_genomic.fna.gz", "genome": "_genomic.fna.gz"}

	def _partition_url(self, assembly: str) -> str:
		prefix, digits = assembly.split("_")
		digits = digits.split(".")[0]
		return "{}/{}/{}/{}/{}".format(
			self.BASE, prefix, digits[0:3], digits[3:6], digits[6:9])

	def _versioned_dir(self, assembly: str) -> Optional[str]:
		for _ in range(3):
			try:
				self.limiter.wait()
				r = self.session.get(self._partition_url(assembly) + "/", timeout=30)
				r.raise_for_status()
				self.failures.pop(assembly, None)
				for name in re.findall(r'href="([^"/]+)/"', r.text):
					if name.startswith(assembly):
						return name
				return None  
			except Exception as e:
				self.failures[assembly] = "{}: {}".format(type(e).__name__, e)
				time.sleep(0.5)
		return None

	def _fetch_one(self, assembly: str) -> GenomeFiles:
		vdir = self._versioned_dir(assembly)
		if vdir is None:
			return GenomeFiles()
		base = "{}/{}/{}".format(self._partition_url(assembly), vdir, vdir)
		return self._fetch_files({
			slot: (base + self.SUFFIX[slot],
				   os.path.join(self.out_dir, assembly + self.SUFFIX[slot]))
			for slot in self._slots()})


class MgnifyGenomeDownloader(_GenomeDownloader): 
	API = "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/{}"
	API_DOWNLOADS = "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/{}/downloads"
	ACCESSION_RE = re.compile(r"^MGYG\d+$")
	SUFFIX = {"gff": ".gff", "faa": ".faa", "genome": ".fna"}

	@classmethod
	def is_mgnify_accession(cls, assembly: str) -> bool:
		return bool(cls.ACCESSION_RE.match(assembly))

	@staticmethod
	def _download_url(entry) -> Optional[str]:
		"""Pull a URL out of one download record without pinning field names."""
		if isinstance(entry, str):
			return entry
		if not isinstance(entry, dict):
			return None
		for key in ("url", "link", "self", "href", "download_url"):
			value = entry.get(key)
			if isinstance(value, str) and value.startswith("http"):
				return value
		links = entry.get("links")
		if isinstance(links, dict):
			for key in ("self", "download", "related"):
				value = links.get(key)
				if isinstance(value, str) and value.startswith("http"):
					return value
		return None

	@classmethod
	def _url_map(cls, payload) -> Dict[str, str]:
		if isinstance(payload, dict):
			records = payload.get("downloads")
			if records is None:
				records = payload.get("items", payload.get("data", []))
		else:
			records = payload or []
		urls = {}
		for entry in records if isinstance(records, list) else []:
			url = cls._download_url(entry)
			if not url:
				continue
			name = os.path.basename(url.split("?", 1)[0].rstrip("/"))
			if name:
				urls[name] = url
		return urls

	def _lookup(self, assembly: str) -> Dict[str, str]:
		t0 = time.monotonic()
		r = self.session.get(self.API.format(assembly), timeout=30)
		debug("MGnify API v2 {} -> HTTP {} in {:.2f}s".format(
			assembly, r.status_code, time.monotonic() - t0))
		r.raise_for_status()
		payload = r.json()
		if flags_log.DEBUG and isinstance(payload, dict):
			debug("  genome keys: {}".format(sorted(payload)[:20]))
		urls = self._url_map(payload)
		if not urls:
			self.limiter.wait()
			r = self.session.get(self.API_DOWNLOADS.format(assembly), timeout=30)
			debug("MGnify downloads {} -> HTTP {}".format(assembly, r.status_code))
			r.raise_for_status()
			urls = self._url_map(r.json())
		return urls

	def _fetch_one(self, assembly: str) -> GenomeFiles:
		try:
			self.limiter.wait()
			urls = self._lookup(assembly)
			self.failures.pop(assembly, None)
		except Exception as e:
			self.failures[assembly] = "{}: {}".format(type(e).__name__, e)
			debug("MGnify lookup failed for {}: {!r}".format(assembly, e), exc=True)
			return GenomeFiles()
		if flags_log.DEBUG:
			debug("MGnify {} offers: {}".format(assembly, sorted(urls)))

		jobs = {}
		for slot in self._slots():
			base = assembly + self.SUFFIX[slot]
			name = next((n for n in (base + ".gz", base) if n in urls), None)
			if name is None:
				debug("MGnify {} {} NOT LISTED".format(assembly, slot))
				jobs[slot] = (None, os.path.join(self.out_dir, base))
				continue
			debug("MGnify {} {} -> {}".format(assembly, slot, urls[name]))
			jobs[slot] = (urls[name], os.path.join(self.out_dir, name))
		return self._fetch_files(jobs)


class LocalGenomeResolver: 
	GFF_EXT = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
	FAA_EXT = (".faa", ".faa.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz")
	RNA_EXT = (".rna.fna", ".rna.fna.gz", ".rna.fa", ".rna.fa.gz",
			   "_rna_from_genomic.fna", "_rna_from_genomic.fna.gz")
	GENOME_EXT = (".fna", ".fna.gz")

	def __init__(self, directory: str):
		self.genomes = self._pair_files(directory)
		self._protein_index = None

	def _pair_files(self, directory):
		found = {"gff": {}, "faa": {}, "rna": {}, "genome": {}}
		try:
			names = sorted(os.listdir(directory))
		except OSError:
			return {}
		claimed = set()

		for slot, exts in (("rna", self.RNA_EXT), ("gff", self.GFF_EXT),
						   ("faa", self.FAA_EXT), ("genome", self.GENOME_EXT)):
			for name in names:
				path = os.path.join(directory, name)
				if name in claimed or not os.path.isfile(path):
					continue
				base = self._basename(name, exts)
				if base is not None:
					found[slot][base] = path
					claimed.add(name)

		return {b: GenomeFiles(found["gff"][b], found["faa"][b],
							   found["rna"].get(b), found["genome"].get(b))
				for b in found["gff"] if b in found["faa"]}

	INFIX = ("_genomic", "_protein", "_rna_from_genomic", "_cds_from_genomic")

	@staticmethod
	def _basename(name, exts):
		for ext in sorted(exts, key=len, reverse=True):
			if name.endswith(ext):
				stem = name[:-len(ext)]
				for infix in LocalGenomeResolver.INFIX:
					if stem.endswith(infix):
						stem = stem[:-len(infix)]
						break
				return stem
		return None

	def _build_protein_index(self):
		index = {}
		for base, files in self.genomes.items():
			with NeighborhoodExtractor._open(files.faa) as fh:
				for rec in SeqIO.parse(fh, "fasta"):
					index.setdefault(rec.id, base)
		self._protein_index = index

	def resolve_pair(self, assembly: str) -> Optional[Tuple[str, GenomeFiles]]:
		if assembly in self.genomes:
			return assembly, self.genomes[assembly]
		for base in self.genomes:
			if base.startswith(assembly) or assembly.startswith(base):
				return base, self.genomes[base]
		return None

	def resolve_protein(self, protein: str) -> Optional[Tuple[str, GenomeFiles]]:
		if self._protein_index is None:
			self._build_protein_index()
		base = self._protein_index.get(protein)
		return (base, self.genomes[base]) if base else None


class FlankingGene(NamedTuple):
	accession: str   # protein/RNA accession, or '<biotype>*' for unidentified non-coding
	strand: str      # '+'/'-', normalized relative to the query gene
	start: int
	end: int
	product: str
	offset: int      # 0 = query, negative = upstream, positive = downstream
	query: str       # the query protein this neighbor belongs to
	is_rna: bool = False   # True for tRNA/rRNA/ncRNA genes (keep the RNA outline)
	contig: str = ""       # contig/sequence id the gene lies on (for locus matching)


class RangeInfo(NamedTuple):
	row: str
	query: str
	assembly: str
	contig: str
	contig_length: int
	q_start: int
	q_end: int
	q_strand: str
	up_available: int      # bp of contig upstream of the query, in query orientation
	down_available: int
	up_reached: int        # bp actually spanned by the extracted genes
	down_reached: int
	genes_up: int
	genes_down: int
	scan_start: int
	scan_end: int


class NeighborhoodExtractor: 
	def __init__(self, flank: int = 4, label_assembly: bool = False,
				 range_bp: Optional[int] = None, scan_range: Optional[int] = None):
		self.flank = flank
		self.range_bp = range_bp
		self.scan_range = scan_range
		self.label_assembly = label_assembly
		self.sequences: Dict[str, str] = {}        # all flanking proteins: accession -> sequence
		self.query_sequences: Dict[str, str] = {}  # query proteins only: accession -> sequence
		self.row_sequences: Dict[str, str] = {}    # row id -> query sequence (one leaf per assembly-row)
		self.rna_sequences: Dict[str, str] = {}    # flanking RNAs: accession -> nucleotide sequence
		self.rna_products: Dict[str, str] = {}     # flanking RNAs: accession -> product name
		self.species: Dict[str, str] = {}          # row id -> organism name
		self.row_label: Dict[str, str] = {}        # row id -> display label
		self.ranges: Dict[str, RangeInfo] = {}     # row id -> window actually obtained
		self.drawn_accessions: Dict[str, set] = {} # row id -> the genes actually drawn
		self.window_genes: Dict[str, list] = {}    # row id -> every gene in the scan window
		self.window_sequences: Dict[str, str] = {} # those genes' proteins, by accession
		self._gff_cache: Dict[str, List[dict]] = {}
		self._faa_cache: Dict[str, Dict[str, Tuple[str, str]]] = {}
		self._rna_cache: Dict[str, Dict[str, str]] = {}
		self._genome_cache: Tuple[Optional[str], Dict[str, str]] = (None, {})
		self._contig_len: Dict[str, Dict[str, int]] = {}
		self._contig_span: Dict[str, Dict[str, Tuple[int, int]]] = {}
		self._reach: Dict[str, List[int]] = {}

	def extract(self, assembly: str, gff_path: str, faa_path: str,
				query: str, acceptable: Optional[set] = None,
				rna_path: Optional[str] = None,
				genome_path: Optional[str] = None) -> List[FlankingGene]:
		acceptable = acceptable or {query}
		genes = self._genes(assembly, gff_path)
		idx = next((i for i, g in enumerate(genes) if g["accession"] in acceptable), None)
		if idx is None:
			return []

		faa = self._faa(assembly, faa_path)
		rna = self._rna(assembly, rna_path) if rna_path else {}
		contig = genes[idx]["contig"]
		qstrand = genes[idx]["strand"]
		lo, hi = self._window(assembly, genes, idx, contig)
		q_acc = genes[idx]["accession"]
		q_organism = faa.get(q_acc, (None, ""))[1]
		row_id = "{}|{}".format(query, assembly)
		name = "{}|{}".format(query, assembly) if self.label_assembly else query
		self.row_label[row_id] = "{}  {}".format(name, q_organism) if q_organism else name

		neighborhood = []
		for j in range(lo, hi):
			g = genes[j]
			if g["contig"] != contig:
				continue
			acc = g["accession"]
			if g["is_rna"]:
				self.rna_products[acc] = g["product"]
				rseq = rna.get(g["locus_tag"]) or rna.get(acc)
				if rseq is None and genome_path:

					rseq = self._slice_genome(
						self._genome(assembly, genome_path),
						g["contig"], g["start"], g["end"], g["strand"])
				if rseq is not None:
					self.rna_sequences[acc] = rseq
			else:
				seq, organism = faa.get(acc, (None, ""))
				if seq is not None:
					self.sequences[acc] = seq
					if j == idx:
						self.query_sequences[acc] = seq
						self.row_sequences[row_id] = seq   
						self.species[row_id] = organism
			offset = j - idx
			if qstrand == "-":
				offset = -offset
			neighborhood.append(FlankingGene(
				accession=acc,
				strand=self._norm_strand(qstrand, g["strand"]),
				start=g["start"], end=g["end"],
				product=g["product"],
				offset=offset,
				query=row_id,
				is_rna=g["is_rna"],
				contig=g["contig"],
			))
		if neighborhood:
			self.drawn_accessions[row_id] = {g.accession for g in neighborhood}
			self.ranges[row_id] = self._range_info(
				row_id, query, assembly, genes, idx, contig, qstrand, neighborhood)
			if self.scan_range:
				self._collect_window(row_id, assembly, genes, idx, contig,
									 qstrand, faa)
		return neighborhood

	def _collect_window(self, row_id, assembly, genes, idx, contig, qstrand, faa):
		lo, hi = self._window(assembly, genes, idx, contig, self.scan_range)
		found = []
		for j in range(lo, hi):
			g = genes[j]
			if g["contig"] != contig or g["is_rna"]:
				continue
			acc = g["accession"]
			seq = faa.get(acc, (None, None))[0]
			if not seq:
				continue
			self.window_sequences.setdefault(acc, seq)
			found.append(FlankingGene(
				accession=acc,
				strand=self._norm_strand(qstrand, g["strand"]),
				start=g["start"], end=g["end"], product=g["product"],
				offset=(idx - j) if qstrand == "-" else (j - idx),
				query=row_id, is_rna=False, contig=contig))
		self.window_genes[row_id] = found

	def forget(self, assembly: str):
		"""Drop an assembly's parsed GFF and proteins once nothing else needs
		them. Contig lengths stay: they are small and the scanning tools and the
		range report both read them after extraction is over."""
		self._gff_cache.pop(assembly, None)
		self._faa_cache.pop(assembly, None)
		self._rna_cache.pop(assembly, None)
		self._reach.pop(assembly, None)
		if self._genome_cache[0] == assembly:
			self._genome_cache = (None, {})

	def contig_lengths(self, assembly: str) -> Dict[str, int]:
		return self._contig_len.get(assembly, {})

	def _window(self, assembly: str, genes: List[dict], idx: int,
				contig: str, span: Optional[int] = None) -> Tuple[int, int]:
		c_lo, c_hi = self._contig_span.get(assembly, {}).get(
			contig, (0, len(genes)))
		reach_bp = span if span is not None else self.range_bp
		if not reach_bp:
			return max(c_lo, idx - self.flank), min(c_hi, idx + self.flank + 1)
		lo_bp = genes[idx]["start"] - reach_bp
		hi_bp = genes[idx]["end"] + reach_bp
		reach = self._reach.get(assembly) or []
		lo = idx
		while lo > c_lo and (reach[lo - 1] if reach else genes[lo - 1]["end"]) >= lo_bp:
			lo -= 1
		hi = idx + 1
		while hi < c_hi and genes[hi]["start"] <= hi_bp:
			hi += 1
		return lo, hi

	def _range_info(self, row_id, query, assembly, genes, idx, contig,
					qstrand, neighborhood) -> RangeInfo:
		q = genes[idx]
		length = self._contig_len.get(assembly, {}).get(contig)
		if not length:
			c_lo, c_hi = self._contig_span.get(assembly, {}).get(contig, (0, len(genes)))
			length = max((g["end"] for g in genes[c_lo:c_hi]), default=q["end"])
		left_avail, right_avail = q["start"] - 1, max(0, length - q["end"])
		lo_bp = min(g.start for g in neighborhood)
		hi_bp = max(g.end for g in neighborhood)
		left_reach, right_reach = q["start"] - lo_bp, hi_bp - q["end"]
		before = sum(1 for g in neighborhood if g.offset < 0)
		after = sum(1 for g in neighborhood if g.offset > 0)
		if qstrand == "-":
			left_avail, right_avail = right_avail, left_avail
			left_reach, right_reach = right_reach, left_reach
		span = self.scan_range
		scan_start = max(1, q["start"] - span) if span else 0
		scan_end = min(length, q["end"] + span) if span else 0
		return RangeInfo(
			row=row_id, query=query, assembly=assembly, contig=contig,
			contig_length=length, q_start=q["start"], q_end=q["end"],
			q_strand=qstrand, up_available=left_avail, down_available=right_avail,
			up_reached=left_reach, down_reached=right_reach,
			genes_up=before, genes_down=after,
			scan_start=scan_start, scan_end=scan_end)

	@staticmethod
	def _open(path: str):
		if path.endswith(".gz"):
			return gzip.open(path, "rt", encoding="utf-8", errors="replace")
		return open(path, "rt", encoding="utf-8", errors="replace")

	def _genes(self, assembly: str, gff_path: str) -> List[dict]:
		if assembly in self._gff_cache:
			return self._gff_cache[assembly]
		genes = []
		lengths: Dict[str, int] = {}
		with self._open(gff_path) as fh:
			for raw in fh:
				if raw.startswith("#"):
					if raw.startswith("##sequence-region"):
						parts = raw.split()
						if len(parts) >= 4:
							try:
								lengths[parts[1]] = int(parts[3])
							except ValueError:
								pass
					continue
				col = raw.rstrip("\n").split("\t")
				if len(col) < 9:
					continue
				feature, attrs = col[2], self._attrs(col[8])
				if feature == "region":
					if col[0] not in lengths:
						try:
							lengths[col[0]] = int(col[4])
						except ValueError:
							pass
				elif feature.endswith("gene"):
					genes.append(self._record(col, None, "",
						biotype=self._attr(attrs, "gene_biotype") or "",
						locus_tag=self._attr(attrs, "locus_tag") or "", is_rna=False))
				elif feature == "CDS":
					locus_tag = self._attr(attrs, "locus_tag") or ""
					accession = self._cds_accession(attrs, locus_tag)
					product = self._attr(attrs, "product") or ""
					pseudo = (self._attr(attrs, "pseudo") or "").lower() == "true"
					if genes and genes[-1]["accession"] is None:
						if pseudo or genes[-1]["biotype"] == "pseudogene":
							genes[-1]["biotype"] = "pseudogene"
						else:
							genes[-1]["accession"] = accession
						genes[-1]["product"] = product
					else:
						genes.append(self._record(col, accession, product,
							biotype="pseudogene" if pseudo else "protein_coding",
							locus_tag=locus_tag, is_rna=False))
				elif feature.endswith("RNA"):

					locus_tag = self._attr(attrs, "locus_tag") or ""
					accession = (self._attr(attrs, "Name") or self._attr(attrs, "transcript_id")
								 or self._strip_id_prefix(self._attr(attrs, "ID")) or locus_tag)
					product = self._attr(attrs, "product") or ""
					if genes and genes[-1]["accession"] is None:
						genes[-1]["accession"] = accession
						genes[-1]["product"] = product
						genes[-1]["is_rna"] = True
						if not genes[-1]["locus_tag"]:
							genes[-1]["locus_tag"] = locus_tag
					else:
						genes.append(self._record(col, accession, product,
							biotype=feature, locus_tag=locus_tag, is_rna=True))
		for g in genes:
			if g["biotype"] == "pseudogene":
				g["accession"] = "pseudogene*"
			elif not g["accession"]:
				g["accession"] = (g["biotype"] or "noProtein") + "*"

		genes.sort(key=lambda g: (g["contig"], g["start"]))
		spans: Dict[str, Tuple[int, int]] = {}
		reach: List[int] = []
		best = 0
		for i, g in enumerate(genes):
			c = g["contig"]
			if c in spans:
				spans[c] = (spans[c][0], i + 1)
				best = max(best, g["end"])
			else:
				spans[c] = (i, i + 1)
				best = g["end"]
			reach.append(best)
		self._contig_len[assembly] = lengths
		self._contig_span[assembly] = spans
		self._reach[assembly] = reach
		self._gff_cache[assembly] = genes
		return genes

	@staticmethod
	def _record(col: List[str], accession: Optional[str], product: str,
				biotype: str, locus_tag: str, is_rna: bool) -> dict:
		return {
			"contig": col[0], "start": int(col[3]), "end": int(col[4]),
			"strand": col[6], "accession": accession, "product": product,
			"biotype": biotype, "locus_tag": locus_tag, "is_rna": is_rna,
		}

	def _faa(self, assembly: str, faa_path: str) -> Dict[str, Tuple[str, str]]:
		if assembly in self._faa_cache:
			return self._faa_cache[assembly]
		table = {}
		with self._open(faa_path) as fh:
			for rec in SeqIO.parse(fh, "fasta"):
				m = re.search(r"\[([^\]]+)\]\s*$", rec.description)
				organism = m.group(1) if m else ""
				table[rec.id] = (str(rec.seq), organism)
		self._faa_cache[assembly] = table
		return table

	def _rna(self, assembly: str, rna_path: str) -> Dict[str, str]:

		if assembly in self._rna_cache:
			return self._rna_cache[assembly]
		table = {}
		with self._open(rna_path) as fh:
			for rec in SeqIO.parse(fh, "fasta"):
				seq = str(rec.seq)
				m = re.search(r"\[locus_tag=([^\]]+)\]", rec.description)
				if m:
					table[m.group(1)] = seq
				table[rec.id] = seq
		self._rna_cache[assembly] = table
		return table

	def _genome(self, assembly: str, genome_path: str) -> Dict[str, str]:

		cached_for, table = self._genome_cache
		if cached_for == assembly:
			return table
		table = {}
		with self._open(genome_path) as fh:
			for rec in SeqIO.parse(fh, "fasta"):
				table[rec.id] = str(rec.seq)
		self._genome_cache = (assembly, table)
		return table

	@staticmethod
	def _slice_genome(genome_seqs: Dict[str, str], contig: str,
					   start: int, end: int, strand: str) -> Optional[str]:
		seq = genome_seqs.get(contig)
		if seq is None:
			return None
		sub = seq[start - 1:end] 
		if strand == "-":
			sub = str(Seq(sub).reverse_complement())
		return sub or None

	@staticmethod
	def _attr(attributes, key: str) -> Optional[str]:
		if isinstance(attributes, dict):
			return attributes.get(key)
		return NeighborhoodExtractor._attrs(attributes).get(key)

	@staticmethod
	def _attrs(attributes: str) -> Dict[str, str]:
		out: Dict[str, str] = {}
		for field in attributes.split(";"):
			key, sep, value = field.partition("=")
			if sep and key not in out:
				out[key] = value       # first wins, as the old regex did
		return out

	@classmethod
	def _cds_accession(cls, attrs, locus_tag: str) -> Optional[str]:
		return (cls._attr(attrs, "protein_id")
				or cls._strip_id_prefix(cls._attr(attrs, "ID"))
				or locus_tag
				or cls._attr(attrs, "Name"))

	@staticmethod
	def _strip_id_prefix(id_attr: Optional[str]) -> Optional[str]:
		if not id_attr:
			return None
		for prefix in ("cds-", "rna-"):
			if id_attr.startswith(prefix):
				return id_attr[len(prefix):]
		return id_attr

	@staticmethod
	def _norm_strand(query_strand: str, gene_strand: str) -> str:
		if query_strand == "+":
			return gene_strand
		return "-" if gene_strand == "+" else "+"


class NeighborhoodClusterer:
	CHUNK = 100               # most queries in one jackhmmer call
	CHUNKS_PER_WORKER = 16    # so the tail of the run is not one coarse round

	def __init__(self, iterations: int = 3, incE: float = 1e-3,
				 workers: Optional[int] = None):
		self.iterations = iterations
		self.incE = incE
		self.workers = workers
		self.alphabet = Alphabet.amino()
		self.adjacency: Dict[str, set] = {}

	def cluster(self, sequences: Dict[str, str]) -> List[List[str]]:
		if not sequences:
			return []

		digital = {name: TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
				   for name, seq in sequences.items()}
		block = DigitalSequenceBlock(self.alphabet, list(digital.values()))

		def search_chunk(chunk):
			names = [n for n, _ in chunk]
			results = pyhmmer.hmmer.jackhmmer(
				[q for _, q in chunk], block,
				max_iterations=self.iterations,
				incE=self.incE,
				cpus=1,
			)
			return [(name, {self._name(h) for h in r.hits if h.included})
					for name, r in zip(names, results)]

		items = sorted(digital.items())
		# Chunking exists only to keep the number of jackhmmer calls down,
		# because each call leaves about 2 MB behind and thousands of them
		# exhaust memory. It costs time at the tail: the last round runs with
		# only (chunks mod workers) workers busy while the rest idle for a whole
		# chunk, so coarse chunks waste more. Enough chunks per worker keeps
		# that under about 5%, and below a few hundred queries it falls to one
		# query per chunk, which is exactly what the run did before chunking
		# existed and is the fastest thing at that size.
		slots = max(self.workers or 1, 1) * self.CHUNKS_PER_WORKER
		size = max(1, min(self.CHUNK, len(items) // slots))
		chunks = [items[i:i + size] for i in range(0, len(items), size)]
		adjacency = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			for part in pool.map(search_chunk, chunks):
				adjacency.update(part)

		self.adjacency = adjacency
		return self._connected_components(adjacency)

	@staticmethod
	def _name(hit) -> str:
		n = hit.name
		return n.decode() if isinstance(n, bytes) else n

	@staticmethod
	def _connected_components(adjacency: Dict[str, set]) -> List[List[str]]:
		symmetric: Dict[str, set] = {node: set() for node in adjacency}
		for node, hits in adjacency.items():
			for hit in hits:
				if hit in symmetric:
					symmetric[node].add(hit)
					symmetric[hit].add(node)
		seen, families = set(), []
		for node in symmetric:
			if node in seen:
				continue
			stack, component = [node], set()
			while stack:
				x = stack.pop()
				if x in component:
					continue
				component.add(x)
				seen.add(x)
				stack.extend(symmetric.get(x, set()) - component)
			families.append(sorted(component))
		families.sort(key=lambda f: (-len(f), f[0]))
		return families


class RnaClusterer: 
	def __init__(self, incE: float = 1e-3, workers: Optional[int] = None):
		self.incE = incE
		self.workers = workers
		self.alphabet = Alphabet.dna()
		self.adjacency: Dict[str, set] = {}

	def cluster(self, sequences: Dict[str, str]) -> List[List[str]]:
		if not sequences:
			return []
		digital = {name: TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
				   for name, seq in sequences.items()}
		block = DigitalSequenceBlock(self.alphabet, list(digital.values()))

		def search_one(item):
			name, query = item
			hits = list(pyhmmer.hmmer.nhmmer([query], block, incE=self.incE, cpus=1))[0]
			return name, {NeighborhoodClusterer._name(h) for h in hits if h.included}

		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			adjacency = dict(pool.map(search_one, digital.items()))
		self.adjacency = adjacency
		return NeighborhoodClusterer._connected_components(adjacency)

	@staticmethod
	def cluster_by_name(products: Dict[str, str]) -> List[List[str]]:
		groups: Dict[str, List[str]] = {}
		for acc, product in products.items():
			key = RnaClusterer._normalise(product) or acc
			groups.setdefault(key, []).append(acc)
		families = [sorted(v) for v in groups.values()]
		families.sort(key=len, reverse=True)
		return families

	@staticmethod
	def _normalise(product: str) -> str:
		return " ".join(product.lower().split())


class ReportWriter: 
	def __init__(self, neighborhoods, families, species,
				 queries, protein_to_assemblies, matched,
				 order=None, adjacency=None, sequences=None, row_sequences=None,
				 ranges=None, requested=None):
		self.ranges = ranges or {}
		self.requested = requested
		self.neighborhoods = neighborhoods
		self.families = families
		self.species = species
		self.queries = queries
		self.protein_to_assemblies = protein_to_assemblies
		self.matched = matched
		self.adjacency = adjacency or {}
		self.sequences = sequences or {}
		self.row_sequences = row_sequences or {}
		rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
		queries = {g.accession for g in neighborhoods if g.offset == 0}
		self.occurrences = Counter(g.accession for g in neighborhoods)
		self.fam_of = family_numbers(families, rna_accessions, queries,
									 self.occurrences)
		self.by_query = {}
		for g in neighborhoods:
			self.by_query.setdefault(g.query, []).append(g)
		if order:
			rank = {row: i for i, row in enumerate(order)}
			self.by_query = {row: self.by_query[row] for row in
							 sorted(self.by_query, key=lambda r: rank.get(r, len(rank)))}
		self.products = {}
		for g in neighborhoods:
			self.products.setdefault(g.accession, g.product)

	def write_all(self, out_path):
		self.operon_tsv(out_path("_operon.tsv"))
		self.clusters_tsv(out_path("_clusters.tsv"))
		self.outdesc_txt(out_path("_outdesc.txt"))
		self.species_info(out_path("_speciesInfo.txt"))
		self.query_status(out_path("_QueryStatus.txt"))
		self.flankgene_report(out_path("_flankgene_Report.log"))
		self.fasta_outputs(out_path)
		if self.ranges:
			self.range_report(out_path("_rangeReport.tsv"))
		if self.adjacency:
			self.jackhits_tsv(out_path("_jackhits.tsv"))
		return self.accession_issues(out_path("_accessionIssues.txt"))

	@staticmethod
	def _split_row(row_id):
		query, _, assembly = row_id.partition("|")
		return query, assembly

	def operon_tsv(self, path):
		with open(path, "w") as out:
			out.write("#query\tassembly\tspecies\tfamily\tstrand\toffset\t"
					  "start\tend\tlength\tcontig\tis_rna\taccession\tproduct\n")
			for row_id in self.by_query:
				query, assembly = self._split_row(row_id)
				sp = self.species.get(row_id, "")
				for g in sorted(self.by_query[row_id], key=lambda x: x.offset):
					out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
						query, assembly or "-", sp, self.fam_of.get(g.accession, "-"),
						g.strand, g.offset, g.start, g.end, g.end - g.start + 1,
						g.contig or "-", "True" if g.is_rna else "False",
						g.accession, g.product))

	def clusters_tsv(self, path):
		with open(path, "w") as out:
			out.write("#family\tsize\tmembers\n")
			for fam in self.families:
				label = (self.fam_of.get(fam[0], "-")
						 if family_shared(fam, self.occurrences) else "-")
				out.write("{}\t{}\t{}\n".format(label, len(fam), ",".join(fam)))

	def outdesc_txt(self, path):
		blocks = [fam for fam in self.families
				  if family_shared(fam, self.occurrences)]
		blocks.sort(key=lambda fam: (-sum(self.occurrences.get(a, 0) for a in fam),
									 fam[0]))
		with open(path, "w") as out:
			for fam in blocks:
				label = self.fam_of.get(fam[0], "-")
				for acc in fam:
					out.write("{}({})\t{}\t{}\n".format(
						label, self.occurrences.get(acc, 0), acc,
						self.products.get(acc, "")))
				out.write("\n\n")

	def fasta_outputs(self, out_path):
		if not (self.sequences or self.row_sequences):
			return
		query_accs = {g.accession for g in self.neighborhoods if g.offset == 0}
		flanking = [a for a in sorted(self.sequences) if a not in query_accs]

		def label(acc):
			product = self.products.get(acc, "")
			return "{}|{}".format(acc, product) if product else acc

		rows_out = [(row, self.row_sequences[row]) for row in sorted(self.row_sequences)]
		flank_out = [(label(a), self.sequences[a]) for a in sorted(flanking)]
		self._write_fasta(out_path("_tree.fasta"), rows_out)
		self._write_fasta(out_path("_flankgene.fasta"), flank_out)
		self._write_fasta(out_path("_all.fasta"), rows_out + flank_out)

	def range_report(self, path):
		want = self.requested
		with open(path, "w") as out:
			out.write("#query\tassembly\tcontig\tcontig_length\tquery_start\t"
					  "query_end\tquery_strand\trequested_bp\tup_available\t"
					  "down_available\tup_reached\tdown_reached\ttruncated\t"
					  "genes_up\tgenes_down\tgenes_total\tscan_start\tscan_end\t"
					  "scan_span\n")
			for row_id in self.by_query:
				info = self.ranges.get(row_id)
				if info is None:
					continue
				sides = []
				if want:
					if info.up_available < want:
						sides.append("up")
					if info.down_available < want:
						sides.append("down")
				out.write("\t".join(str(v) for v in (
					info.query, info.assembly or "-", info.contig or "-",
					info.contig_length, info.q_start, info.q_end, info.q_strand,
					want if want else "-",
					info.up_available, info.down_available,
					info.up_reached, info.down_reached,
					",".join(sides) if sides else "-",
					info.genes_up, info.genes_down,
					info.genes_up + info.genes_down + 1,
					info.scan_start or "-", info.scan_end or "-",
					(info.scan_end - info.scan_start + 1) if info.scan_start else "-",
					)) + "\n")

	@staticmethod
	def _write_fasta(path, records):
		with open(path, "w") as out:
			for name, seq in records:
				out.write(">{}\n{}\n".format(name, seq))

	def jackhits_tsv(self, path):
		with open(path, "w") as out:
			out.write("#accession\tfamily\tn_hits\thits\n")
			for acc in sorted(self.adjacency):
				hits = sorted(self.adjacency[acc])
				out.write("{}\t{}\t{}\t{}\n".format(
					acc, self.fam_of.get(acc, "-"), len(hits), ";".join(hits)))

	def species_info(self, path):
		with open(path, "w") as out:
			out.write("#query\tassembly\tspecies\n")
			for row_id in sorted(self.species):
				query, assembly = self._split_row(row_id)
				out.write("{}\t{}\t{}\n".format(query, assembly or "-",
											    self.species[row_id]))

	def query_status(self, path):
		with open(path, "w") as out:
			out.write("#query\tassemblies\tflanking_genes_found\n")
			for q in self.queries:
				asms = self.protein_to_assemblies.get(q, [])
				status = "Yes" if q in self.matched else "No"
				out.write("{}\t{}\t{}\n".format(q, ";".join(asms) if asms else "-", status))

	def flankgene_report(self, path):
		with open(path, "w") as out:
			for query in self.by_query:
				genes = sorted(self.by_query[query], key=lambda x: x.offset)
				chain = " ".join("{}({})".format(g.accession, self.fam_of.get(g.accession, "-"))
								 for g in genes)
				out.write("{}\t{}\n".format(query, chain))

	def accession_issues(self, path):
		lines = []
		for q in self.queries:
			if not self.protein_to_assemblies.get(q):
				lines.append("{}\tno assembly resolved (not found locally or via NCBI)".format(q))
			elif q not in self.matched:
				lines.append("{}\tassembly resolved but no flanking neighborhood extracted".format(q))
		with open(path, "w") as out:
			out.write("#query\tissue\n")
			for line in lines:
				out.write(line + "\n")
		return len(lines)

def family_shared(fam, occurrences=None):
	if occurrences:
		return sum(occurrences.get(acc, 0) for acc in fam) > 1
	return len(fam) > 1


def family_numbers(families, rna_accessions=None, query_accessions=None,
				   occurrences=None):
	rna_accessions = rna_accessions or set()
	query_accessions = query_accessions or set()
	number = {}
	prot_n, rna_n, query_n = 0, 0, 0
	shared = [fam for fam in families if family_shared(fam, occurrences)]
	if occurrences:
		shared.sort(key=lambda fam: (-sum(occurrences.get(a, 0) for a in fam),
									 fam[0]))
	for fam in shared:
		if fam[0] in rna_accessions:
			rna_n += 1
			label = "R{}".format(rna_n)
		elif query_accessions.intersection(fam):
			query_n += 1
			label = "Q{}".format(query_n)
		else:
			prot_n += 1
			label = str(prot_n)
		for acc in fam:
			number[acc] = label
	return number


VERSION = "2.2.2"

DEFAULT_INTERPRO = "interpro_metadata_processed.tsv"


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

class InstanceLock:

	def __init__(self, path: str):
		self.path = path
		self.acquired = False

	def acquire(self) -> Optional[str]:
		"""Return None on success, or a description of the holder."""
		for attempt in (1, 2):
			try:
				fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
			except FileExistsError:
				holder = self._read()
				if attempt == 1 and self._is_stale(holder):
					debug("removing stale lock {} ({})".format(self.path, holder))
					try:
						os.unlink(self.path)
						continue
					except OSError:
						pass
				return holder
			except OSError as e:
				debug("could not create lock {}: {!r}".format(self.path, e))
				return None   # unwritable location is not a reason to refuse to run
			with os.fdopen(fd, "w") as out:
				out.write("pid {}\nhost {}\nstarted {}\n".format(
					os.getpid(), platform.node(),
					time.strftime("%Y-%m-%d %H:%M:%S")))
			self.acquired = True
			debug("acquired lock {}".format(self.path))
			return None
		return self._read()

	def release(self):
		if not self.acquired:
			return
		try:
			os.unlink(self.path)
			debug("released lock {}".format(self.path))
		except OSError as e:
			debug("could not remove lock {}: {!r}".format(self.path, e))
		self.acquired = False

	def _read(self) -> str:
		try:
			with open(self.path) as fh:
				return " / ".join(fh.read().split("\n")).strip(" /")
		except OSError:
			return "unknown"

	def _is_stale(self, holder: str) -> bool:
		match = re.search(r"pid (\d+)", holder)
		host = re.search(r"host (\S+)", holder)
		if not match or (host and host.group(1) != platform.node()):
			return False   # another machine: cannot tell, so assume it is live
		try:
			os.kill(int(match.group(1)), 0)
		except ProcessLookupError:
			return True
		except PermissionError:
			return False   # exists, owned by someone else
		except OSError:
			return False
		return False

SECRET_ARGS = ("api_key",)
SECRET_FLAGS = ("--api_key",)


def _redact_command_line(argv):
	"""Mask secret values in both '--flag value' and '--flag=value' forms."""
	out, skip = [], False
	for i, token in enumerate(argv):
		if skip:
			out.append("<given>")
			skip = False
			continue
		flag = token.split("=", 1)[0]
		if flag in SECRET_FLAGS:
			out.append(flag + "=<given>" if "=" in token else token)
			skip = "=" not in token
		else:
			out.append(token)
	return " ".join(out)


def note_skipped(args, tool, reason):
	path = os.path.join(args.output, os.path.basename(
		os.path.normpath(args.output)) + "_runinfo.txt")
	try:
		with open(path, "a") as out:
			out.write("\nskipped {}\n  {}\n".format(tool, reason))
	except OSError:
		pass


def write_run_info(args, parser):
	os.makedirs(args.output, exist_ok=True)
	prefix = os.path.basename(os.path.normpath(args.output))
	out_path = lambda suffix: os.path.join(args.output, prefix + suffix)

	defaults = {a.dest: a.default for a in parser._actions}
	given, default = [], []
	for dest in sorted(vars(args)):
		value = getattr(args, dest)
		shown = "<given>" if dest in SECRET_ARGS and value else value
		(given if value != defaults.get(dest) else default).append(
			"  {:18} {}".format(dest, shown))

	with open(out_path("_runinfo.txt"), "w") as out:
		out.write("FlaGs3 {}\n".format(VERSION))
		out.write("run started   {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S %Z")))
		out.write("host          {}\n".format(platform.node()))
		out.write("python        {}\n".format(sys.version.split()[0]))
		out.write("working dir   {}\n".format(os.getcwd()))
		out.write("\ncommand line\n  {}\n".format(_redact_command_line(sys.argv)))
		out.write("\noptions set explicitly\n")
		out.write("\n".join(given) + "\n" if given else "  (none)\n")
		out.write("\noptions left at default\n")
		out.write("\n".join(default) + "\n" if default else "  (none)\n")

	sources = [(p, "_input.txt" if i == 0 else "_input{}.txt".format(i + 1))
			   for i, p in enumerate(list(args.input_list or []))]
	if args.blast_input:
		sources.append((args.blast_input, "_blast_input.txt"))
	args.input_copies = [suffix for _, suffix in sources]
	for source, suffix in sources:
		try:
			shutil.copyfile(source, out_path(suffix))
		except OSError as e:
			print("Warning: could not copy {} into the output directory "
				  "({}).".format(source, e))


def build_parser():
	usage = ''' Description: Identify flanking genes and cluster them based on similarity. Requirement= Python3, BioPython, pyhmmer, requests. '''
	parser = argparse.ArgumentParser(description=usage)
	parser.add_argument("-i", "--input_list", action="append", metavar="FILE", help=" Protein Accession eg. WP_047256880.1, optionally tab-separated with an assembly Identifier eg. GCF_000001765.3 (NCBI RefSeq/GenBank) or MGYG000454827 (MGnify Genomes catalogue -- the protein accession must then be the exact locus tag used in that genome's annotation). One per line. May be given more than once. Optional if --blast_input or a legacy list is given. ")
	parser.add_argument("-bi", "--blast_input", help=" File holding ONE starting point: a RefSeq protein accession (WP_/NP_/YP_/XP_/AP_), or a protein sequence as FASTA or bare residues over any number of lines. BlastP finds its homologues and they become the queries, appended to anything in -i. ")
	parser.add_argument("-bm", "--blast_mode", choices=("remote", "local"), default="remote", help=" Where to run BlastP. 'remote' uses NCBI QBLAST -- no install, but a search takes minutes. 'local' runs blastp from NCBI BLAST+ against a local database, which is far faster but needs the binary and the database. Default = remote ")
	parser.add_argument("-bd", "--blast_db", default="refseq_select", help=" Database to search. Aliases: refseq_select (representative RefSeq proteins), refseq_protein (full RefSeq), genbank (nr), swissprot. Each resolves to the right name for the chosen mode. Any other value passes through unchanged, which is how a local database path is given. Default = refseq_select ")
	parser.add_argument("-be", "--blast_evalue", type=float, default=1e-5, help=" E-value cutoff for BlastP. Default = 1e-5 ")
	parser.add_argument("-bw", "--blast_wait", type=float, default=60, help=" Give up waiting on NCBI's queue after this many minutes. The job usually still finishes on NCBI's side and the printed link stays valid. Default = 60 ")
	parser.add_argument("-bh", "--blast_hits", type=int, default=50, help=" How many BlastP hits to carry forward as queries. Each one becomes a genome download and a row in the figure, so this is the main control on how big the run gets. No upper limit, but webFlaGs caps at 200 and remote QBLAST may return fewer than asked. Default = 50 ")
	parser.add_argument("-u", "--user_email", required=True, help=" User Email Address (required by NCBI Entrez). ")
	parser.add_argument("-api", "--api_key", help=" NCBI API Key. ")
	parser.add_argument("-g", "--gene", type=int, default=4, help=" Number of flanking genes up/downstream. Default = 4 ")
	parser.add_argument("-r", "--range", type=int, metavar="BP", help=" Take the neighbourhood as every gene within this many bases of the query gene, instead of a fixed number of genes. Measured outwards from the query gene's own start and end, and a gene straddling the edge is included. Overrides -g. Contigs shorter than the window are reported in <dir>_rangeReport.tsv rather than silently truncated. ")
	parser.add_argument("-sm", "--scan_margin", type=int, default=10000, metavar="BP", help=" Extra sequence handed to the scanning tools beyond the analysis range, so a system straddling the edge is still called whole rather than cut in half. Hits reaching into the margin are kept and marked 'partial' in the coverage column. Set 0 to scan exactly the analysis range. Default = 10000 ")
	parser.add_argument("-sr", "--scan_range", type=int, metavar="BP", help=" Genomic span around the query gene given to the tools that scan sequence rather than genes (Sismis, geNomad). Without it they get the whole genome, as before. Independent of -g/-r, so -g 5 -sr 50000 clusters 11 genes and scans 100 kb. A tool can override it in the scan_range column of tools_table.tsv. ")
	parser.add_argument("-m", "--max_assemblies", type=int, default=1, help=" Max assemblies per protein. Default = 1 ")
	parser.add_argument("-rm", "--remap", action="store_true", help=" Look up again, through IPG, any protein whose paired assembly in the input produced nothing. A protein given with an assembly otherwise never goes near IPG, so a withdrawn assembly, a failed download, or an accession that is simply not in that assembly leaves it unresolved with no second attempt. Off by default because it costs extra NCBI requests and downloads; worth it when the assembly column was compiled some time ago. ")
	parser.add_argument("-nc", "--no_cross_db", action="store_true", help=" Keep protein and genome in the same database: RefSeq proteins (WP_, NP_, ...) resolve only to GCF_ assemblies and INSDC proteins only to GCA_. Proteins whose only assemblies are in the other database are then reported as unresolved. Does not affect assemblies given explicitly in the input file. ")
	parser.add_argument("-e", "--ethreshold", type=float, default=1e-3, help=" Jackhmmer inclusion E-value threshold for clustering flanking genes. Default = 1e-3 ")
	parser.add_argument("-n", "--number", type=int, default=3, help=" Number of jackhmmer iterations for clustering. Default = 3 ")
	parser.add_argument("-cc", "--cluster_collapse", nargs="?", const="0.9,0.8", metavar="ID[,COV]", help=" Collapse near-identical flanking proteins with MMseqs2 before clustering, run jackhmmer on the representatives only, and give every member its representative's family. Clustering cost grows with the square of the number of proteins, so this is what makes a run of thousands of queries finish. Takes a minimum sequence identity and optionally a coverage, both fractions; bare -cc means 0.9,0.8, which only merges sequences jackhmmer would certainly have joined. Every representative and its members are written to <dir>_collapse.tsv. Needs mmseqs on PATH; run mmseqs_installer.sh. ")
	parser.add_argument("-c", "--cpu", type=int, help=" Max parallel CPU workers (default: auto-detect). ")
	parser.add_argument("-tmp", "--temporary", default="./genomes", help=" Temporary directory for downloaded assemblies; deleted at the end. Default = ./genomes ")
	parser.add_argument("-o", "--output", default="output", help=" Directory for result files; its name is also the file prefix. A YYYYMMDD_HHMMSS stamp of the run start is appended, so repeated runs do not overwrite each other. Default = output ")
	parser.add_argument("-nt", "--no_timestamp", action="store_true", help=" Use -O/--output verbatim instead of appending a date-time stamp. Repeated runs then overwrite each other; useful for scripted pipelines that need a fixed path. ")
	parser.add_argument("-t", "--tree", action="store_true", help=" Also build a phylogenetic tree and write it as <dir>_tree.nwk. Figures that use it are controlled by the figure table. ")
	parser.add_argument("-tm", "--trimal_mode", default="gt", help=" trimal column-filter mode: gt, cons, st, or a preset such as gappyout, strict, strictplus, automated1, nogaps, noallgaps. Default = gt ")
	parser.add_argument("-tv", "--trimal_value", type=float, default=0.1, help=" Value for the trimal mode that takes one (gt, cons, st). Default = 0.1 ")
	parser.add_argument("-tx", "--trimal_extra", default="", help=" Extra trimal arguments, passed through verbatim, e.g. \"-w 3\". ")
	parser.add_argument("-iq", "--iqtree", action="store_true", help=" Build the tree with IQ-TREE (ModelFinder plus 1000 ultrafast bootstrap replicates) instead of VeryFastTree. Implies --tree. Much slower but gives model selection and branch support. Needs iqtree on PATH. ")
	parser.add_argument("-to", "--tree_order", action="store_true", help=" Order the neighbours output by tree leaf order (implies --tree). ")
	parser.add_argument("-d", "--domains", action="store_true", help=" Scan flanking proteins for domains and write <dir>_domains.tsv (requires --hmmdb). ")
	parser.add_argument("-db", "--hmmdb", action="append", metavar="[NAME=]PATH", help=" HMM database for domain scanning: a .hmm file, or a directory of .hmm files such as DefenseFinder's profiles/. Repeat for several. Prefix with NAME= to label it in the outputs, otherwise the file or directory name is used. Models carrying a gathering threshold are scored by it; the rest use --ethreshold. Default = ./pfam_db/Pfam-A.hmm ")
	parser.add_argument("-hc", "--hmm_coverage", action="append", metavar="[NAME=]Q[,H]", help=" Minimum fraction of the protein (Q) and of the model (H) an alignment must span, dropping partial hits. Give NAME= to apply it to one database, or omit NAME to apply it to all. Sensible for full-length protein models such as DefenseFinder (e.g. 0.7,0.5); leave off for domain databases like Pfam, where partial coverage is normal. ")
	parser.add_argument("-pdf", "--pdf", action="store_true", help=" Also write a PDF beside every figure. Needs one of cairosvg, svglib, rsvg-convert or inkscape; without one the figures are still written as SVG. ")
	parser.add_argument("-nf", "--no_figures", action="store_true", help=" Run the analysis and write the tables, but draw nothing. Figures can be produced later with flags_redraw.py. ")
	parser.add_argument("-fh", "--figure_height", type=int, default=16383, metavar="PX", help=" Split a figure into parts when it would be taller than this, writing <name>_part1.svg and so on. The default is the canvas limit of Illustrator and librsvg, past which a figure opens in a browser but nowhere else. Raise it if you only ever view figures in a browser. Figures carrying a tree panel are never split, since the tree spans every row. ")
	parser.add_argument("-no", "--no_overlaps", action="store_true", help=" Leave a family number out of a gene too small to hold it. By default the number is drawn anyway, since knowing which family a small gene belongs to usually matters more than the overlap. ")
	parser.add_argument("-f", "--figures", metavar="TSV", help=" Figure table controlling which figures are drawn and every parameter of how. Default: visualisation_table.tsv, written into the output directory for you to edit and re-apply with flags_redraw.py. ")
	parser.add_argument("-cl", "--clans", help=" Pfam-A.clans.tsv(.gz): colour domains by clan instead of family. ")
	parser.add_argument("-ip", "--interpro", default=DEFAULT_INTERPRO, help=" InterPro metadata table (.tsv or .tsv.gz) with 'accession' and 'pfam_members' columns. Adds the InterPro entry, its name and type, and short characterisation/informativeness summaries to _domains.tsv, joined on the Pfam accession. Looked for in the working directory and next to FlaGs3.py; if it is not there the domain table is written without those columns. Default = " + DEFAULT_INTERPRO + " ")
	parser.add_argument("-lth", "--local_tmhmm", action="store_true", help=" Predict transmembrane regions with DeepTMHMM on this machine rather than in the BioLib cloud. Implies --tmhmm. The command and its directory come from the deeptmhmm row of tools_table.tsv; run deeptmhmm_installer.sh to set that up. ")
	parser.add_argument("-lsp", "--local_signalp", action="store_true", help=" Predict signal peptides with SignalP on this machine rather than in the BioLib cloud. Implies --signalp. The command and its directory come from the signalp row of tools_table.tsv; run signalp_installer.sh to set that up. ")
	parser.add_argument("--tools", metavar="TSV", help=" Table of external tool commands (mafft, trimal, VeryFastTree, IQ-TREE, blastp, sismis, DeepTMHMM, SignalP). Default: the tools_table.tsv next to FlaGs3.py. ")
	parser.add_argument("-th", "--tmhmm", action="store_true", help=" Predict transmembrane regions (DeepTMHMM, via BioLib cloud) and draw them as double red dotted lines in the domain figure. Off by default; needs pybiolib and network. ")
	parser.add_argument("-sp", "--signalp", action="store_true", help=" Predict signal peptides (SignalP-6, via BioLib cloud) and draw them as black triangles in the domain figure. Off by default; needs pybiolib and network. ")
	parser.add_argument("-df", "--defensefinder", action="store_true", help=" Call anti-phage defence systems with DefenseFinder and draw each one as a band spanning its genes, labelled with the system name. Runs on the neighbourhood genes, so it follows -g/-r and ignores -sr. Needs defensefinder_installer.sh. ")
	parser.add_argument("-pl", "--padloc", action="store_true", help=" Call anti-phage defence systems with PadLoc, drawn the same way as DefenseFinder's. Both can run together; a system both tools agree on is drawn once and the called_by column names both. Needs padloc_installer.sh. ")
	parser.add_argument("-gn", "--genomad", action="store_true", help=" Scan for proviruses and plasmids with geNomad and draw them as bands like Sismis' secretion systems. Each band is labelled with what geNomad actually called it -- the virus taxon, or plasmid/conjugative plasmid -- not a generic 'mobile element'. Uses -sr like Sismis does, so it sees a window around each query rather than whole genomes. Needs genomad and its database; run genomad_installer.sh then genomad_loader.sh. ")
	parser.add_argument("-gdb", "--genomad_db", metavar="DIR", help=" geNomad database directory, overriding the db column of the genomad row in tools_table.tsv. ")
	parser.add_argument("-ss", "--sismis", action="store_true", help=" Scan each assembly's genomic FASTA for secretion systems using Sismis (github.com/lmc297/Sismis) and write <dir>_secretion.tsv, noting which query neighborhoods (if any) each hit overlaps. Downloads the genomic FASTA per assembly. Off by default; needs sismis installed (pip install sismis). ")
	parser.add_argument("-k", "--keep", action="store_true", help=" Keep the downloaded assemblies instead of deleting the temporary directory at the end. ")
	parser.add_argument("-ul", "--use_local", metavar="DIR", help=" Directory of local .gff/.faa genome files to search before falling back to NCBI. Files may be gzipped; a genome is a .gff and .faa sharing a basename. ")
	parser.add_argument("-cr", "--cluster_rna", action="store_true", help=" Cluster flanking RNA genes into families too (off by default). Uses RNA sequences (nhmmer) when available, otherwise groups by product name. ")
	parser.add_argument("-vb", "--verbose", action="store_true", help=" Print progress and a stage-by-stage funnel. ")
	parser.add_argument("-dbg", "--debug", action="store_true", help=" Print diagnostics to stderr: HTTP status codes and Retry-After headers, external command lines and exit codes, and full tracebacks for errors that are otherwise only summarised. Implies --verbose. ")
	parser.add_argument("-nl", "--no_lock", action="store_true", help=" Skip the lock that stops two FlaGs3 runs sharing one -tmp directory. Two runs writing the same genome files also double the request rate against NCBI and EBI, which gets you throttled, so only use this with separate -tmp directories. ")
	parser.add_argument("-v", "--version", action="version", version="FlaGs3 " + VERSION)
	return parser


class Scans(NamedTuple):
	module: object
	hits: list
	rows: dict
	statuses: dict
	features: dict
	scanner: object = None
	genomad: object = None
	genomad_hits: list = ()
	genomad_statuses: dict = {}
	genomad_scanner: object = None
	defence: object = None
	defence_hits: list = ()
	defence_statuses: dict = {}
	defence_replicons: int = 0


def flags_tools_span(args, tool):
	import flags_tools
	span = flags_tools.scan_range(tool, args.scan_range)
	return "{}bp".format(span) if span else "genome"


def scan_windows(args, extractor, mod, tool):
	import flags_tools
	import flags_scan
	span = flags_tools.scan_range(tool, args.scan_range)
	if not span:
		return {}
	spans = {}
	for row_id, info in extractor.ranges.items():
		assembly = row_id.rsplit("|", 1)[-1]
		query = row_id.rsplit("|", 1)[0]
		# use the span the range report recorded, which is already clipped to
		# the contig; recomputing it here claimed an analysed range running
		# past the end of a short contig
		if info.scan_start and info.scan_end:
			lo, hi = info.scan_start, info.scan_end
		else:
			lo = max(1, info.q_start - span)
			hi = min(info.q_end + span, info.contig_length or info.q_end + span)
		spans.setdefault(assembly, []).append((info.contig, lo, hi))
		# so each cut window can name the queries it was cut for
		flags_scan.WINDOW_QUERIES.setdefault(
			(assembly, info.contig, lo, hi), []).append(query)
		flags_scan.WINDOW_QUERY_POS[query] = (info.q_start, info.q_end,
											  info.q_strand)
	return {assembly: flags_scan.merge_windows(
			entries, args.scan_margin, extractor.contig_lengths(assembly))
			for assembly, entries in spans.items()}


def scan_outcome(statuses):
	ok, failed, skipped = 0, [], []
	for assembly, status in statuses.items():
		if status.startswith("error"):
			failed.append((assembly, status))
		elif status.startswith("skipped"):
			skipped.append((assembly, status))
		else:
			ok += 1
	return ok, failed, skipped


def report_scan(label, statuses, hits, noun, verbose):
	ok, failed, skipped = scan_outcome(statuses)
	if failed and not ok:
		print("Warning: {} failed on every genome. First error: {}".format(
			label, failed[0][1]))
		print("         See <prefix>_{}_diagnostics.txt, and the "
			  "'--- {} (exit N) ---' block in <prefix>_console.log for the "
			  "tool's own message.".format(label.lower(), label.lower()))
	elif failed:
		print("Warning: {} failed on {} of {}; the rest were scanned. "
			  "See <prefix>_{}_diagnostics.txt.".format(
				  label, len(failed), plural(len(statuses), "genome"),
				  label.lower()))
	if verbose and ok:
		print(">> {}: scanned {}, found {}{}".format(
			label, plural(ok, "genome"), plural(len(hits), noun),
			"" if not skipped else
			" ({} skipped)".format(len(skipped))), flush=True)


def row_spans(all_neighborhoods):
	rows = {}
	for g in all_neighborhoods:
		assembly = g.query.rsplit("|", 1)[-1]
		if g.query in rows:
			_, contig, lo, hi = rows[g.query]
			rows[g.query] = (assembly, contig, min(lo, g.start), max(hi, g.end))
		else:
			rows[g.query] = (assembly, g.contig, g.start, g.end)
	return rows


def run_background_scans(args, extractor, downloaded, all_neighborhoods, timings):
	genomad_box = {"mod": None, "hits": [], "statuses": {}, "scanner": None}
	defence_box = {"mod": None, "hits": [], "statuses": {}, "replicons": 0}
	sismis_scanner = []
	sismis_mod = None
	sismis_hits: List = []
	sismis_rows: Dict[str, Tuple[str, str, int, int]] = {}
	sismis_statuses: Dict[str, str] = {}
	features = {}

	def _run_sismis():
		nonlocal sismis_mod
		t = time.perf_counter()
		rows, hits, statuses = {}, [], {}
		if not (args.sismis and all_neighborhoods):
			return hits, rows, statuses, 0.0
		try:
			import flags_secretion as mod
		except ImportError:
			print("Warning: --sismis needs the sismis package (pip install sismis); skipping secretion-system detection.")
			return hits, rows, statuses, time.perf_counter() - t
		for g in all_neighborhoods:
			assembly = g.query.rsplit("|", 1)[-1]
			if g.query in rows:
				_, contig, lo, hi = rows[g.query]
				rows[g.query] = (assembly, contig, min(lo, g.start), max(hi, g.end))
			else:
				rows[g.query] = (assembly, g.contig, g.start, g.end)
		span = flags_tools_span(args, "sismis")
		scanner = mod.SismisScanner(out_dir=os.path.join(args.output, "sismis"))
		scanner.windows_dir = os.path.join(args.output, "windows")
		scanner.window_key = span
		windows = scan_windows(args, extractor, mod, "sismis")
		jobs = []
		for assembly in sorted({asm for asm, _, _, _ in rows.values()}):
			genome_path = downloaded.get(assembly, GenomeFiles()).genome
			if not genome_path:
				statuses[assembly] = "skipped: no genomic FASTA downloaded"
				continue
			jobs.append((assembly, genome_path, windows.get(assembly, [])))
		if jobs:
			hits = scanner.scan(jobs, statuses)
		sismis_mod = mod
		sismis_scanner.append(scanner)
		return hits, rows, statuses, time.perf_counter() - t

	def _run_defence():
		t = time.perf_counter()
		wanted = [n for n, flag in (("defensefinder", args.defensefinder),
									("padloc", args.padloc)) if flag]
		if not (wanted and all_neighborhoods):
			return [], {}, 0, 0.0
		try:
			import flags_defence as mod
		except ImportError as e:
			print("Warning: defence system calling could not be loaded ({}); "
				  "skipping.".format(e))
			return [], {}, 0, time.perf_counter() - t
		wide = [g for row in extractor.window_genes.values() for g in row]
		replicons = mod.build_replicons(wide or all_neighborhoods)
		scanner = mod.DefenceScanner(out_dir=os.path.join(args.output, "defence"),
									 threads=args.cpu or 0)
		hits = []
		for tool in wanted:
			ok, where = scanner.available(tool)
			if not ok:
				print("Warning: --{} needs {}; skipping.".format(tool, where))
				note_skipped(args, tool, where)
				scanner.statuses[tool] = "skipped: {}".format(where)
				continue
			try:
				hits.extend(scanner.run(
					tool, replicons,
					extractor.window_sequences or extractor.sequences))
			except Exception as e:
				scanner.statuses[tool] = "error: {}".format(e)
		defence_box["mod"] = mod
		return (mod.merge_calls(hits), scanner.statuses, len(replicons),
				time.perf_counter() - t)

	def _run_genomad():
		t = time.perf_counter()
		statuses, hits = {}, []
		if not (args.genomad and all_neighborhoods):
			return hits, statuses, 0.0
		try:
			import flags_genomad as mod
		except ImportError as e:
			print("Warning: --genomad could not be loaded ({}); skipping.".format(e))
			return hits, statuses, time.perf_counter() - t
		import flags_tools
		database = args.genomad_db or flags_tools.get("genomad")[1]
		scanner = mod.GenomadScanner(
			out_dir=os.path.join(args.output, "genomad"),
			database=database or "", threads=args.cpu or 0,
			windows_dir=os.path.join(args.output, "windows"),
			window_key=flags_tools_span(args, "genomad"))
		ok, where = scanner.available()
		if not ok:
			print("Warning: --genomad needs geNomad and {}; skipping.".format(where))
			note_skipped(args, "genomad", where)
			return hits, statuses, time.perf_counter() - t
		windows = scan_windows(args, extractor, mod, "genomad")
		jobs = []
		for assembly in sorted({asm for asm, _, _, _ in row_spans(all_neighborhoods).values()}):
			genome_path = downloaded.get(assembly, GenomeFiles()).genome
			if not genome_path:
				statuses[assembly] = "skipped: no genomic FASTA downloaded"
				continue
			jobs.append((assembly, genome_path, windows.get(assembly, [])))
		if jobs:
			hits = scanner.scan(jobs, statuses)
			if args.verbose:
				print(">> geNomad: {} over {}".format(
					plural(getattr(scanner, "batches", 1), "invocation"),
					plural(len(jobs), "genome")), flush=True)
		genomad_box["mod"] = mod
		genomad_box["scanner"] = scanner
		return hits, statuses, time.perf_counter() - t

	def _run_tmhmm():
		t = time.perf_counter()
		if not args.tmhmm:
			return {}, 0.0
		try:
			import flags_features as feat_mod
			tm = feat_mod.TMScanner(local=args.local_tmhmm).scan(
				extractor.sequences, want_signal=not args.signalp)
		except ImportError:
			print("Warning: --tmhmm needs pybiolib (pip install pybiolib); skipping transmembrane prediction.")
			return {}, time.perf_counter() - t
		except Exception as e:
			print("Warning: DeepTMHMM did not finish, skipping transmembrane regions ({}).".format(e))
			note_skipped(args, "tmhmm", e)
			return {}, time.perf_counter() - t
		return tm, time.perf_counter() - t

	def _run_signalp():
		t = time.perf_counter()
		if not args.signalp:
			return {}, 0.0
		try:
			import flags_features as feat_mod
			sp = feat_mod.SignalPScanner(local=args.local_signalp).scan(
				extractor.sequences)
		except ImportError:
			print("Warning: --signalp needs pybiolib (pip install pybiolib); skipping signal-peptide prediction.")
			return {}, time.perf_counter() - t
		except Exception as e:
			print("Warning: SignalP did not finish, skipping signal peptides ({}).".format(e))
			note_skipped(args, "signalp", e)
			return {}, time.perf_counter() - t
		return sp, time.perf_counter() - t

	tm, sp = {}, {}
	active = [n for n, flag in (("sismis", args.sismis), ("genomad", args.genomad),
								 ("defence", args.defensefinder or args.padloc),
								 ("tmhmm", args.tmhmm),
								 ("signalp", args.signalp)) if flag]
	if active:
		if args.verbose:
			local = [n for n, on in (("tmhmm", args.local_tmhmm),
									 ("signalp", args.local_signalp))
					 if on and n in active]
			where = "locally: {}".format(", ".join(local)) if local else \
				"cloud/subprocess, not local CPU"
			print(">> running {} in the background ({})...".format(
				", ".join(active), where), flush=True)
		if ((args.tmhmm and not args.local_tmhmm)
				or (args.signalp and not args.local_signalp)):
			import flags_features as feat_mod
			feat_mod.warm_up()
		task = {"sismis": _run_sismis, "genomad": _run_genomad,
				"defence": _run_defence,
				"tmhmm": _run_tmhmm, "signalp": _run_signalp}
		with ThreadPoolExecutor(max_workers=len(task)) as pool:
			futures = {pool.submit(task[name]): name for name in active}
			for fut in as_completed(futures):
				name = futures[fut]
				if name == "sismis":
					sismis_hits, sismis_rows, sismis_statuses, elapsed = fut.result()
					timings["sismis_scan"] = elapsed
					if sismis_mod:
						report_scan("sismis", sismis_statuses, sismis_hits,
									"secretion system", args.verbose)
				elif name == "genomad":
					genomad_box["hits"], genomad_box["statuses"], elapsed = fut.result()
					timings["genomad_scan"] = elapsed
					report_scan("geNomad", genomad_box["statuses"],
								genomad_box["hits"], "mobile element", args.verbose)
				elif name == "defence":
					(defence_box["hits"], defence_box["statuses"],
					 defence_box["replicons"], elapsed) = fut.result()
					timings["defence_scan"] = elapsed
					bad = [t for t, v in defence_box["statuses"].items()
						   if v.startswith("error")]
					if bad and not defence_box["hits"]:
						print("Warning: {} produced nothing but an error. See "
							  "<prefix>_defence_diagnostics.txt and the console "
							  "log.".format(" and ".join(sorted(bad))))
					elif args.verbose:
						print(">> defence systems: {} across {}".format(
							plural(len(defence_box["hits"]), "system"),
							plural(defence_box["replicons"], "neighbourhood")), flush=True)
				elif name == "tmhmm":
					tm, elapsed = fut.result()
					timings["tmhmm_scan"] = elapsed
					if args.verbose:
						print(">> DeepTMHMM: features on {} proteins".format(len(tm)), flush=True)
				elif name == "signalp":
					sp, elapsed = fut.result()
					timings["signalp_scan"] = elapsed
					if args.verbose:
						print(">> SignalP: signal peptides on {} proteins".format(len(sp)), flush=True)
		for acc, regs in tm.items():
			features.setdefault(acc, []).extend(regs)
		for acc, regs in sp.items():
			features.setdefault(acc, []).extend(regs)

	return Scans(sismis_mod, sismis_hits, sismis_rows, sismis_statuses, features,
				 sismis_scanner[0] if sismis_scanner else None,
				 genomad_box["mod"], genomad_box["hits"],
				 genomad_box["statuses"], genomad_box["scanner"],
				 defence_box["mod"], defence_box["hits"],
				 defence_box["statuses"], defence_box["replicons"])


def resolve_blast(args, proteins_assembly, proteins_only, inline, timings, t0):
	blast_hits, blast_mod, queries = [], None, []
	if not args.blast_input and not inline:
		return blast_hits, queries, blast_mod
	import flags_blast as blast_mod

	if args.blast_input:
		try:
			queries.append(blast_mod.read_query(args.blast_input))
		except (OSError, ValueError) as e:
			sys.exit("Error: {}".format(e))
	for text, path, n in inline:
		try:
			queries.append(blast_mod.parse_query(
				[text], "line {} of {}".format(n, path)))
		except ValueError as e:
			sys.exit("Error: {}".format(e))

	if args.blast_mode == "remote":
		print(">> waiting on NCBI's queue; this often takes several minutes "
			  "and can be much longer when NCBI is busy.", flush=True)
	Entrez.email = args.user_email
	Entrez.tool = NCBI_TOOL
	if args.api_key:
		Entrez.api_key = args.api_key
	searcher = blast_mod.BlastSearcher(
		mode=args.blast_mode, database=args.blast_db,
		evalue=args.blast_evalue, max_hits=args.blast_hits,
		threads=args.cpu or 0, email=args.user_email,
		max_wait=args.blast_wait * 60,
		report=lambda msg: print(">> {}".format(msg), flush=True))

	for query in queries:
		label = query.accession or "the supplied sequence"
		if args.verbose:
			print(">> BlastP ({}) for {} against {}...".format(
				args.blast_mode, label, args.blast_db), flush=True)
		try:
			hits = searcher.search(query)
		except Exception as e:
			debug("blast failed", exc=True)
			sys.exit("Error: BlastP search failed for {}: {}".format(label, e))
		if not hits:
			print("Warning: BlastP returned no hits for {}. Try a larger "
				  "--blast_evalue or a fuller --blast_db.".format(label))
			continue
		already = {p.split(".")[0] for p, _ in proteins_assembly}
		already |= {p.split(".")[0] for p in proteins_only}
		added = [h.accession for h in hits
				 if h.accession.split(".")[0] not in already]
		proteins_only.extend(added)
		blast_hits.extend(hits)
		if args.verbose:
			short = (" (asked for {}; the database had no more above the E-value "
					 "cutoff)".format(args.blast_hits)
					 if len(hits) < args.blast_hits else "")
			print(">> BlastP: {} hits, {} new queries{}".format(
				len(hits), len(added), short), flush=True)

	if queries and not blast_hits:
		sys.exit("Error: BlastP returned no hits for any query.")
	timings["1b_blast"] = time.perf_counter() - t0
	return blast_hits, queries, blast_mod


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


def write_window_tables(extractor, out_path):
	"""Every gene inside the scan window, not only the drawn neighbourhood. The
	scanning tools work on the whole window, so a defence system can name genes
	that appear nowhere else in the output; without this the hit says which
	accessions it spans and nothing can say what they are."""
	with open(out_path("_window.tsv"), "w") as out:
		out.write("#query\tassembly\tcontig\taccession\tstart\tend\tstrand\t"
				  "offset\tin_neighbourhood\tproduct\n")
		for row in sorted(extractor.window_genes):
			info = extractor.ranges.get(row)
			assembly = row.rsplit("|", 1)[-1]
			inside = extractor.drawn_accessions.get(row, set())
			for g in sorted(extractor.window_genes[row], key=lambda g: g.start):
				out.write("\t".join(str(v) for v in (
					row.rsplit("|", 1)[0], assembly, g.contig, g.accession,
					g.start, g.end, g.strand, g.offset,
					"yes" if g.accession in inside else "no",
					g.product or "")) + "\n")
	with open(out_path("_window.fasta"), "w") as out:
		for row in sorted(extractor.window_genes):
			info = extractor.ranges.get(row)
			query = row.rsplit("|", 1)[0]
			assembly = row.rsplit("|", 1)[-1]
			inside = extractor.drawn_accessions.get(row, set())
			for g in sorted(extractor.window_genes[row], key=lambda g: g.start):
				seq = extractor.window_sequences.get(g.accession)
				if not seq:
					continue
				lo, hi = relative_span(g, info)
				out.write(">{} {}\n{}\n".format(
					g.accession,
					window_header(query, assembly, g, len(seq), lo, hi,
								  g.accession in inside),
					seq))


def relative_span(gene, info):
	"""Where the gene sits relative to the query gene, in query orientation:
	negative upstream, positive downstream, so it reads the same way as offset
	and as the up/down columns of the range report."""
	if info is None:
		return gene.start, gene.end
	if info.q_strand == "-":
		return info.q_end - gene.end, info.q_end - gene.start
	return gene.start - info.q_start, gene.end - info.q_start


def window_header(query, assembly, gene, length, lo, hi, drawn):
	return ("query={} assembly={} contig={} len={} start={} end={} strand={} "
			"offset={} in_neighbourhood={} product={}").format(
				query, assembly, gene.contig, length, lo, hi, gene.strand,
				gene.offset, "yes" if drawn else "no",
				gene.product or "hypothetical protein")


def print_summary(args, prefix, extractor, families, rna_families, figures_written,
				  tree_written, want_tree, domain_table_written, features, sismis_mod,
				  blast_hits, collapse_written=False, genomad_written=False,
				  defence_written=False):
	print("\n{} -> {}".format(
		plural(len(extractor.sequences), "flanking protein"),
		plural(len(families) - len(rna_families), "family", "families")))
	print("\noutputs in {}/".format(os.path.relpath(args.output)
									if args.output.startswith(os.getcwd() + os.sep)
									else args.output))
	for name in figures_written:
		print("  {}".format(name))
	if tree_written:
		print("  tree/ ({}_tree.nwk, alignments, commands)".format(prefix))
	elif want_tree:
		print("  (tree skipped: fewer than 3 query sequences)")
	if domain_table_written:
		print("  {}_domains.tsv".format(prefix))
	if features:
		print("  {}_features.tsv".format(prefix))
	if args.sismis and sismis_mod:
		print("  {}_secretion.tsv / {}_sismis_diagnostics.txt".format(prefix, prefix))
	if genomad_written:
		print("  {}_genomad.tsv / {}_genomad_diagnostics.txt".format(prefix, prefix))
	if defence_written:
		print("  {}_defence.tsv / {}_defence_diagnostics.txt".format(prefix, prefix))
	for suffix in ("_operon.tsv", "_clusters.tsv", "_outdesc.txt", "_speciesInfo.txt",
				   "_QueryStatus.txt", "_flankgene_Report.log", "_jackhits.tsv",
				   "_accessionIssues.txt", "_tree.fasta", "_flankgene.fasta", "_all.fasta",
				   "_runinfo.txt"):
		print("  {}{}".format(prefix, suffix))
	if extractor.ranges:
		print("  {}_rangeReport.tsv".format(prefix))
	if extractor.window_genes:
		print("  {}_window.tsv / {}_window.fasta".format(prefix, prefix))
	if collapse_written:
		print("  {}_collapse.tsv".format(prefix))
	if flags_log.transcript_path():
		print("  {}_console.log".format(prefix))
	for suffix in getattr(args, "input_copies", []):
		print("  {}{}".format(prefix, suffix))
	if blast_hits:
		print("  {}_blast_hits.tsv / {}_blast_accessions.txt".format(prefix, prefix))
		if args.blast_input:
			print("  {}_blast_input.txt".format(prefix))



def main():
	flags_log.start_buffering()
	atexit.register(flags_log.close)
	parser = build_parser()
	args = parser.parse_args()

	args.input_list = list(args.input_list or [])
	all_lists = args.input_list
	if not all_lists and not args.blast_input:
		sys.exit("Error: give -i/--input_list, -bi/--blast_input, or both.")
	for path in all_lists:
		if not os.path.isfile(path):
			sys.exit("Error: input list not found: {}".format(path))
	if args.blast_input and not os.path.isfile(args.blast_input):
		sys.exit("Error: --blast_input file not found: {}".format(args.blast_input))
	if args.blast_input and args.blast_hits < 2:
		sys.exit("Error: --blast_hits must be at least 2.")
	if args.use_local and not os.path.isdir(args.use_local):
		sys.exit("Error: --use_local directory not found: {}".format(args.use_local))
	args.hmmdb = list(args.hmmdb or [])
	if args.domains and not args.hmmdb:
		args.hmmdb = [DEFAULT_HMMDB]
	for spec in args.hmmdb:
		path = spec.partition("=")[2] or spec
		if not os.path.exists(os.path.expanduser(path)):
			sys.exit("Error: --hmmdb path not found: {}".format(path))
	try:
		args.hmm_coverage = parse_coverage(args.hmm_coverage)
	except ValueError as e:
		sys.exit("Error: {}".format(e))
	if args.clans and not os.path.isfile(args.clans):
		sys.exit("Error: --clans file not found: {}".format(args.clans))
	if args.range is not None and args.range < 1:
		sys.exit("Error: --range must be a positive number of bases.")
	if args.scan_range is not None and args.scan_range < 1:
		sys.exit("Error: --scan_range must be a positive number of bases.")
	if args.scan_margin < 0:
		sys.exit("Error: --scan_margin cannot be negative.")
	if args.gene < 0:
		sys.exit("Error: --gene cannot be negative.")
	if args.cluster_collapse:
		try:
			args.cluster_collapse = parse_collapse(args.cluster_collapse)
		except ValueError as e:
			sys.exit("Error: {}".format(e))

	if not args.no_timestamp:
		args.output = "{}_{}".format(os.path.normpath(args.output),
									 time.strftime("%Y%m%d_%H%M%S"))

	args.output = os.path.abspath(args.output)
	args.temporary = os.path.abspath(args.temporary)

	try:
		flags_log.attach(os.path.join(
			args.output,
			os.path.basename(os.path.normpath(args.output)) + "_console.log"))
	except OSError as e:
		print("Warning: could not open the console log ({}); the run continues "
			  "without one.".format(e))

	if args.local_tmhmm:
		args.tmhmm = True
	if args.local_signalp:
		args.signalp = True
	set_debug(args.debug)
	import flags_tools
	try:
		flags_tools.load(args.tools)
	except (OSError, ValueError) as e:
		sys.exit("Error: {}".format(e))
	if args.debug:
		args.verbose = True
		debug("FlaGs3 {} on {} / python {}".format(
			VERSION, platform.platform(), sys.version.split()[0]))

	lock = InstanceLock(args.temporary + ".lock")
	if not args.no_lock:
		holder = lock.acquire()
		if holder:
			sys.exit(
				"Error: another FlaGs3 run is using {}\n"
				"  holder: {}\n"
				"Two runs sharing one temporary directory overwrite each other's "
				"genome files and double the request rate against NCBI and EBI, which "
				"gets the host throttled.\n"
				"Give this run its own directory with -tmp, or pass --no_lock to "
				"override. If no such run exists, delete {}".format(
					args.temporary, holder, lock.path))
		atexit.register(lock.release)

	write_run_info(args, parser)

	timings = {}
	t_start = time.perf_counter()
	t0 = t_start

	proteins_assembly, proteins_only, inline_blast = (
		AccessionListReader(args.input_list).read() if all_lists
		else ([], [], []))
	proteins_only = list(proteins_only)
	timings["1_read_input"] = time.perf_counter() - t0; t0 = time.perf_counter()

	blast_hits, blast_queries, blast_mod = resolve_blast(
		args, proteins_assembly, proteins_only, inline_blast, timings, t0)
	t0 = time.perf_counter()

	all_queries = [p for p, _ in proteins_assembly] + list(proteins_only)
	if args.verbose:
		print(">> read {} queries ({} paired, {} protein-only)".format(
			len(all_queries), len(proteins_assembly), len(proteins_only)), flush=True)

	local = LocalGenomeResolver(args.use_local) if args.use_local else None
	local_files: Dict[str, GenomeFiles] = {}   # source id -> files
	protein_to_assemblies = {}                 # protein -> [source ids]
	local_acceptable = {}                      # protein -> {source id: {accession}}
	pending_pairs = []                         # (protein, assembly) still needing NCBI
	pending_only = []                          # bare proteins still needing NCBI

	if local:
		for protein, assembly in proteins_assembly:
			hit = local.resolve_pair(assembly)
			if hit:
				base, files = hit
				local_files[base] = files
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_pairs.append([protein, assembly])
		for protein in proteins_only:
			hit = local.resolve_protein(protein)
			if hit:
				base, files = hit
				local_files[base] = files
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_only.append(protein)
		if args.verbose:
			print(">> local: {} genomes indexed; resolved {} of {} queries locally".format(
				len(local.genomes), len(protein_to_assemblies),
				len(proteins_assembly) + len(proteins_only)), flush=True)
	else:
		pending_pairs = proteins_assembly
		pending_only = proteins_only

	mapper = ProteinAssemblyMapper(email=args.user_email, api_key=args.api_key,
								   max_assemblies=args.max_assemblies,
								   cross_db=not args.no_cross_db)
	ncbi_map = mapper.map(pending_only)
	for protein, asms in ncbi_map.items():
		protein_to_assemblies.setdefault(protein, []).extend(asms)
	for protein, assembly in pending_pairs:
		protein_to_assemblies.setdefault(protein, []).append(assembly)
	timings["2_ipg_mapping"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose and pending_only:
		print(">> NCBI IPG: mapped {} of {} remaining proteins to assemblies".format(
			len(ncbi_map), len(pending_only)), flush=True)
	if mapper.unreachable:
		print("Warning: could not reach NCBI to resolve {} remaining protein(s) "
			  "({}). They are reported as unresolved; the rest of the run is "
			  "unaffected.".format(len(pending_only), mapper.unreachable))
	if args.verbose and mapper.dropped_cross_db:
		lost = sum(1 for a in mapper.dropped_cross_db if not ncbi_map.get(a))
		print(">> --no_cross_db: excluded cross-database assemblies for {} proteins"
			  "{}".format(len(mapper.dropped_cross_db),
						  ", {} left with none".format(lost) if lost else ""), flush=True)

	assemblies = sorted({asm for asms in protein_to_assemblies.values()
						 for asm in asms if asm not in local_files})
	mgnify_assemblies = [a for a in assemblies if MgnifyGenomeDownloader.is_mgnify_accession(a)]
	ncbi_assemblies = [a for a in assemblies if a not in set(mgnify_assemblies)]
	dl_workers = args.cpu if args.cpu else min(max(len(assemblies), 1), 10)
	dl_rate = 10.0 if args.api_key else 5.0

	def progress(label):
		if not args.verbose:
			return None
		return lambda done, total: print(
			">> {} download: {}/{}".format(label, done, total), flush=True)

	if args.verbose and assemblies:
		print(">> downloading {}...".format(plural(len(assemblies), "genome")), flush=True)
	downloaded: Dict[str, GenomeFiles] = {}
	failures: Dict[str, str] = {}
	if ncbi_assemblies:
		dl = AssemblyDownloader(out_dir=args.temporary, workers=dl_workers, rate=dl_rate,
								want_rna=args.cluster_rna,
								want_genome=args.sismis or args.genomad)
		downloaded.update(dl.download_many(ncbi_assemblies, progress("NCBI")))
		failures.update(dl.failures)
		timings["3a_download_ncbi"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if mgnify_assemblies:
		mg = MgnifyGenomeDownloader(out_dir=args.temporary, workers=dl_workers, rate=dl_rate,
									want_genome=args.sismis or args.genomad or args.cluster_rna)
		downloaded.update(mg.download_many(mgnify_assemblies, progress("MGnify")))
		failures.update(mg.failures)
		timings["3b_download_mgnify"] = time.perf_counter() - t0; t0 = time.perf_counter()
	downloaded.update(local_files)
	if args.verbose:
		ready = sum(1 for f in downloaded.values() if f.gff and f.faa)
		print(">> genomes ready: {} of {} usable ({} NCBI, {} MGnify, {} local)".format(
			ready, len(downloaded), len(ncbi_assemblies), len(mgnify_assemblies),
			len(local_files)), flush=True)
		if failures:

			print(">> {}, first few:".format(plural(len(failures), "download error")), flush=True)
			for key, msg in list(failures.items())[:5]:
				print("     {}: {}".format(key.rsplit("/", 1)[-1] or key, msg), flush=True)

	extractor = NeighborhoodExtractor(flank=args.gene, range_bp=args.range,
									  scan_range=args.scan_range,
									  label_assembly=args.max_assemblies > 1)
	all_neighborhoods = []
	matched = set()
	# Plan the work in input order, then do it grouped by assembly. Each GFF and
	# protein table is parsed once either way, but grouping lets the parsed copy
	# be freed as soon as that assembly is finished, which is the difference
	# between a bounded footprint and tens of gigabytes on a large input.
	plan = []
	for protein, asms in protein_to_assemblies.items():
		for asm in asms:
			files = downloaded.get(asm, GenomeFiles())
			if not (files.gff and files.faa):
				continue
			if asm in local_files:
				acceptable = local_acceptable.get(protein, {}).get(asm)
			else:
				acceptable = mapper.accessions_in.get(protein, {}).get(asm)
			plan.append((protein, asm, acceptable, files))

	grouped: Dict[str, List[int]] = {}
	for i, entry in enumerate(plan):
		grouped.setdefault(entry[1], []).append(i)
	extracted: Dict[int, list] = {}
	for asm, indexes in grouped.items():
		for i in indexes:
			protein, _, acceptable, files = plan[i]
			extracted[i] = extractor.extract(
				asm, files.gff, files.faa, protein, acceptable,
				rna_path=files.rna if args.cluster_rna else None,
				genome_path=files.genome if args.cluster_rna else None)
		extractor.forget(asm)
	for i, (protein, asm, _, _) in enumerate(plan):
		rows = extracted.get(i)
		if rows:
			matched.add(protein)
			all_neighborhoods.extend(rows)

	# A protein paired with an assembly in the input never went through IPG, so
	# a withdrawn assembly, a failed download, or an accession that simply is
	# not in that assembly left it silently unresolved. Retry those through IPG
	# once, which is the only route that can find where the protein actually
	# lives.
	paired = {p for p, _ in proteins_assembly}
	stranded = sorted(p for p in paired if p not in matched)
	if stranded and not args.remap and args.verbose:
		print(">> {} whose paired assembly gave nothing; --remap looks them up "
			  "again through IPG".format(
				  plural(len(stranded), "protein")), flush=True)
	if stranded and args.remap:
		if args.verbose:
			print(">> retrying {} through IPG; their paired assembly gave "
				  "nothing".format(plural(len(stranded), "protein")), flush=True)
		retry = mapper.map(stranded)
		fresh = sorted({a for asms in retry.values() for a in asms}
					   if retry else [])
		fresh = [a for a in fresh if a not in downloaded and a not in local_files]
		if fresh:
			rdl = AssemblyDownloader(out_dir=args.temporary, workers=dl_workers,
									 rate=dl_rate, want_rna=args.cluster_rna,
									 want_genome=args.sismis or args.genomad)
			downloaded.update(rdl.download_many(
				[a for a in fresh if not MgnifyGenomeDownloader.is_mgnify_accession(a)],
				progress("IPG retry")))
			failures.update(rdl.failures)
		regained = 0
		for protein in stranded:
			for asm in retry.get(protein, []):
				files = downloaded.get(asm, GenomeFiles())
				if not (files.gff and files.faa):
					continue
				rows = extractor.extract(
					asm, files.gff, files.faa, protein,
					mapper.accessions_in.get(protein, {}).get(asm),
					rna_path=files.rna if args.cluster_rna else None,
					genome_path=files.genome if args.cluster_rna else None)
				extractor.forget(asm)
				if rows:
					matched.add(protein)
					all_neighborhoods.extend(rows)
					regained += 1
					protein_to_assemblies.setdefault(protein, []).append(asm)
					break
		timings["4b_ipg_retry"] = time.perf_counter() - t0; t0 = time.perf_counter()
		if args.verbose:
			print(">> IPG retry recovered {} of {}".format(
				regained, len(stranded)), flush=True)
	timings["4_extract_neighbors"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		print(">> extracted {} flanking-gene records; {} of {} queries matched".format(
			len(all_neighborhoods), len(matched), len(all_queries)), flush=True)
	if args.range and extractor.ranges:
		short = [r for r in extractor.ranges.values()
				 if r.up_available < args.range or r.down_available < args.range]
		if short:
			print("Note: {} of {} rows sit closer than {} bp to a contig end, so their "
				  "window is truncated; see <prefix>_rangeReport.tsv.".format(
					  len(short), len(extractor.ranges), args.range))

	scans = run_background_scans(args, extractor, downloaded, all_neighborhoods,
								 timings)
	sismis_mod, sismis_hits = scans.module, scans.hits
	sismis_rows, sismis_statuses = scans.rows, scans.statuses
	features = scans.features
	if not args.keep:
		shutil.rmtree(args.temporary, ignore_errors=True)

	if not all_neighborhoods:
		os.makedirs(args.output, exist_ok=True)
		prefix = os.path.basename(os.path.normpath(args.output))
		reporter = ReportWriter(all_neighborhoods, [], extractor.species,
								all_queries, protein_to_assemblies, matched)
		issues_path = os.path.join(args.output, prefix + "_accessionIssues.txt")
		reporter.accession_issues(issues_path)
		reporter.query_status(os.path.join(args.output, prefix + "_QueryStatus.txt"))
		sys.exit("No flanking neighbourhoods could be extracted for any query. "
				 "See {} for per-query details.".format(issues_path))

	t0 = time.perf_counter()
	collapse_mod, collapse_result = None, None
	cluster_input = extractor.sequences
	if args.cluster_collapse:
		import flags_collapse as collapse_mod
		ident, cov = args.cluster_collapse
		collapser = collapse_mod.Collapser(min_seq_id=ident, coverage=cov,
										   threads=args.cpu or 0)
		ok, where = collapser.available()
		if not ok:
			sys.exit(
				"Error: --cluster_collapse needs MMseqs2 and {}.\n"
				"Run mmseqs_installer.sh to install it, or drop --cluster_collapse "
				"-- but note that clustering {} proteins without collapsing is "
				"estimated at {} on {}.".format(
					where, len(extractor.sequences),
					collapse_mod.format_estimate(len(extractor.sequences), args.cpu),
					plural(args.cpu or 1, "worker")))
		if args.verbose:
			print(">> collapsing {} flanking proteins at {:.0%} identity / {:.0%} "
				  "coverage ({})...".format(len(extractor.sequences), ident, cov,
										    where), flush=True)
		try:
			collapse_result = collapser.collapse(extractor.sequences)
		except (OSError, RuntimeError) as e:
			sys.exit("Error: MMseqs2 collapsing failed ({}).".format(e))
		cluster_input = collapse_result.representatives
		timings["5a_collapse"] = time.perf_counter() - t0; t0 = time.perf_counter()
		before, after = len(extractor.sequences), len(cluster_input)
		print(">> collapsed {} proteins into {} representatives ({:.2f}x); "
			  "clustering the representatives".format(before, after,
													  before / max(after, 1)))
		print(">> estimated clustering: {} instead of {} on {}".format(
			collapse_mod.format_estimate(after, args.cpu),
			collapse_mod.format_estimate(before, args.cpu),
			plural(args.cpu or 1, "worker")))
	if args.verbose:
		print(">> clustering {} flanking proteins...".format(len(cluster_input)), flush=True)
	clusterer = NeighborhoodClusterer(iterations=args.number, incE=args.ethreshold,
									  workers=args.cpu)
	families = clusterer.cluster(cluster_input)
	if collapse_result:
		rep_families = len(families)
		families = collapse_mod.expand(families, collapse_result.members)
		if args.verbose:
			print(">> expanded {} representative families to cover {} proteins".format(
				rep_families, sum(len(f) for f in families)), flush=True)

	rna_families = []
	if args.cluster_rna:
		rna = RnaClusterer(incE=args.ethreshold, workers=args.cpu)
		have_seq = set(extractor.rna_sequences)
		all_rna = set(extractor.rna_products)
		missing = all_rna - have_seq
		rna_families = rna.cluster(extractor.rna_sequences)
		if missing:
			fallback = {acc: extractor.rna_products[acc] for acc in missing}
			rna_families += rna.cluster_by_name(fallback)
			print("Warning: {} of {} flanking RNAs had no nucleotide sequence available, "
				  "so they were grouped by product name rather than by sequence."
				  .format(len(missing), len(all_rna)))
		families = families + rna_families
	timings["5_clustering"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		msg = ">> clustered {} into {}".format(
			plural(len(extractor.sequences), "flanking protein"),
			plural(len(families) - len(rna_families), "family", "families"))
		if args.cluster_rna:
			msg += "; {} into {}".format(
				plural(len(extractor.rna_products), "RNA"),
				plural(len(rna_families), "family", "families"))
		print(msg, flush=True)

	want_tree = args.tree or args.tree_order or args.iqtree
	newick = ""
	leaf_order = None
	tree_mod = None
	builder = None
	if want_tree:
		if args.verbose:
			print(">> building tree with {} ({} leaves)...".format(
				"IQ-TREE" if args.iqtree else "VeryFastTree",
				len(extractor.row_sequences)), flush=True)
		import flags_tree as tree_mod
		builder = tree_mod.TreeBuilder(
			threads=args.cpu or 0,
			engine="iqtree" if args.iqtree else "veryfasttree",
			trimal_mode=args.trimal_mode, trimal_value=args.trimal_value,
			trimal_extra=args.trimal_extra)
		newick, _ = builder.build(extractor.row_sequences)
		if newick:
			leaf_order = tree_mod.ladderized_leaf_order(newick)
	if want_tree:
		timings["6_tree"] = time.perf_counter() - t0
	t0 = time.perf_counter()
	if args.verbose and want_tree:
		print(">> tree: {}".format("built ({} leaves)".format(len(leaf_order))
			  if newick else "skipped (fewer than 3 query sequences)"), flush=True)

	os.makedirs(args.output, exist_ok=True)
	prefix = os.path.basename(os.path.normpath(args.output))
	def out_path(suffix):
		return os.path.join(args.output, prefix + suffix)
	order = leaf_order if args.tree_order else None
	tree_written = False
	if newick:
		tree_dir = os.path.join(args.output, "tree")
		os.makedirs(tree_dir, exist_ok=True)
		tree_path = lambda suffix: os.path.join(tree_dir, prefix + suffix)
		with open(tree_path("_tree.nwk"), "w") as out:
			out.write(newick + "\n")
		if builder:
			for suffix, aln in (("_alignment.aln", builder.raw_alignment),
								("_trimmed.aln", builder.alignment)):
				if not aln:
					continue
				with open(tree_path(suffix), "w") as out:
					for name, seq in aln.items():
						out.write(">{}\n{}\n".format(name, seq))
			if builder.commands:
				with open(tree_path("_commands.txt"), "w") as out:
					out.write("\n".join(builder.commands) + "\n")
				with open(out_path("_runinfo.txt"), "a") as out:
					out.write("\ntree commands\n")
					for c in builder.commands:
						out.write("  {}\n".format(c))
		tree_written = True
	timings["6b_tree_files"] = time.perf_counter() - t0; t0 = time.perf_counter()
	want_domain_fig = args.domains or features
	domains, clan_map, domain_table_written = scan_domains(
		args, extractor, families, all_neighborhoods, out_path, timings, t0,
		want_domain_fig)
	t0 = time.perf_counter()
	if features:
		try:
			import flags_features as feat_mod
			feat_mod.write_report(features, out_path("_features.tsv"))
		except Exception as e:
			print("Warning: could not write the feature table ({}).".format(e))
	t0 = time.perf_counter()

	if scans.defence:
		t0 = time.perf_counter()
		dmod = scans.defence
		dmatches = dmod.match_rows(scans.defence_hits, row_spans(all_neighborhoods))
		dmod.write_report(scans.defence_hits, dmatches, out_path("_defence.tsv"))
		dmod.write_diagnostics(scans.defence_statuses, scans.defence_replicons,
							   out_path("_defence_diagnostics.txt"))
		timings["9c_defence_report"] = time.perf_counter() - t0

	if args.genomad and scans.genomad:
		t0 = time.perf_counter()
		gmod = scans.genomad
		gmatches = gmod.match_rows(scans.genomad_hits, row_spans(all_neighborhoods))
		gmod.write_report(scans.genomad_hits, gmatches, out_path("_genomad.tsv"))
		gscanner = scans.genomad_scanner
		gmod.write_diagnostics(
			scans.genomad_statuses, out_path("_genomad_diagnostics.txt"),
			scanned=getattr(gscanner, "scanned_bases", None),
			windows=getattr(gscanner, "window_count", None),
			genome_size={a: sum(extractor.contig_lengths(a).values())
						 for a in scans.genomad_statuses if extractor.contig_lengths(a)})
		timings["9b_genomad_report"] = time.perf_counter() - t0

	if args.sismis and sismis_mod:
		t0 = time.perf_counter()
		matches = sismis_mod.match_rows(sismis_hits, sismis_rows)
		sismis_mod.write_report(sismis_hits, matches, out_path("_secretion.tsv"))
		scanner = scans.scanner
		sismis_mod.write_diagnostics(
			sismis_statuses, out_path("_sismis_diagnostics.txt"),
			scanned=getattr(scanner, "scanned_bases", None),
			windows=getattr(scanner, "window_count", None),
			genome_size={a: sum(extractor.contig_lengths(a).values())
						 for a in sismis_statuses if extractor.contig_lengths(a)})
		timings["9_sismis_report"] = time.perf_counter() - t0

	adjacency = dict(clusterer.adjacency)
	if args.cluster_rna:
		adjacency.update(rna.adjacency)
	reporter = ReportWriter(all_neighborhoods, families, extractor.species,
							all_queries, protein_to_assemblies, matched,
							order=order, adjacency=adjacency,
							sequences=extractor.sequences,
							row_sequences=extractor.row_sequences,
							ranges=extractor.ranges, requested=args.range)
	n_issues = reporter.write_all(out_path)
	if extractor.window_genes:
		write_window_tables(extractor, out_path)
	if collapse_result:
		collapse_mod.write_report(collapse_result, out_path("_collapse.tsv"),
								  reporter.fam_of)
	if blast_hits:
		blast_mod.write_report(blast_hits, blast_queries[0],
							   out_path("_blast_hits.tsv"))
		blast_mod.write_accessions(blast_hits, blast_queries,
								   out_path("_blast_accessions.txt"))
	if args.verbose:
		print(">> wrote data tables and reports ({} queries with issues)".format(n_issues),
			  flush=True)

	t0 = time.perf_counter()
	figures_written = []
	if not args.no_figures:
		redraw = os.path.join(os.path.dirname(os.path.abspath(__file__)),
							  "flags_redraw.py")
		cmd = [sys.executable, redraw, "--data", args.output, "--prefix", prefix,
			   "--max_height", str(args.figure_height)]
		if args.no_overlaps:
			cmd.append("--no_overlaps")
		if args.figures:
			cmd += ["--format", args.figures]
		if args.pdf:
			cmd.append("--pdf")
		if args.verbose:
			cmd.append("--verbose")
		debug("running: {}".format(" ".join(cmd)))
		try:
			done = subprocess.run(cmd, capture_output=True, text=True)
			debug("flags_redraw exit {}".format(done.returncode))
			flags_log.record("--- flags_redraw.py (exit {}) ---".format(done.returncode))
			flags_log.record(done.stdout)
			flags_log.record(done.stderr, "stderr")
			if done.returncode != 0:
				print("Warning: figures were not drawn ({}).".format(
					(done.stderr or done.stdout or "").strip()[:300]))
			else:
				for line in done.stdout.splitlines():
					name = line.strip()
					if name.endswith(".svg"):
						figures_written.append(name)
		except OSError as e:
			print("Warning: could not run {} ({}).".format(redraw, e))
	timings["7_visualize"] = time.perf_counter() - t0

	print_summary(args, prefix, extractor, families, rna_families, figures_written,
				  tree_written, want_tree, domain_table_written, features,
				  sismis_mod, blast_hits, collapse_written=bool(collapse_result),
				  genomad_written=bool(scans.genomad),
				  defence_written=bool(scans.defence))
	if args.verbose:
		print("\n--- timing (seconds) ---")
		for stage in sorted(timings):
			print("  {:24s} {:8.2f}".format(stage, timings[stage]))
		print("  {:24s} {:8.2f}".format("TOTAL", time.perf_counter() - t_start))


if __name__ == '__main__':
	try:
		main()
	except FileNotFoundError as e:
		sys.exit("Error: file not found - {}".format(e))
	except KeyboardInterrupt:
		sys.exit("\nInterrupted.")
	except Exception as e:
		sys.exit("Error: {}".format(e))