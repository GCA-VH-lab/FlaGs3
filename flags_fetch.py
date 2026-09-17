import os
import re
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, NamedTuple, Optional, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import Entrez, SeqIO

import flags_log
from flags_log import debug

NCBI_TOOL = "flags3"
from flags_extract import NeighborhoodExtractor
from flags_report import VERSION, plural


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

	IPG_CHUNK = 200

	def _map_ipg(self, proteins: List[str]) -> Dict[str, set]:
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
