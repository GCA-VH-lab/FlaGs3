import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO

from flags3.log import debug
from flags3.net import Downloader

ENTREZ_TOOL = "flags3"
REFSEQ_PROTEIN = re.compile(r"^[A-Z]{2}_")
ASSEMBLY = re.compile(r"GC[AF]_\d+\.\d+")
SUFFIX = {"gff": "_genomic.gff.gz", "faa": "_protein.faa.gz",
	"rna": "_rna_from_genomic.fna.gz", "genome": "_genomic.fna.gz"}
BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


class Resolution:
	def __init__(self):
		self.assemblies: dict[str, list[str]] = {}
		self.accessions: dict[str, dict[str, set[str]]] = {}
		self.dropped: dict[str, set[str]] = {}
		self.errors: dict[str, str] = {}

	def aliases(self, protein: str, assembly: str) -> set[str]:
		return set(self.accessions.get(protein, {}).get(assembly, ()))


class IpgMapper:
	CHUNK = 200

	def __init__(self, email: str, api_key: Optional[str] = None, max_assemblies: int = 1,
			cross_db: bool = True, pause: float = 0.4):
		self.max_assemblies = max_assemblies
		self.cross_db = cross_db
		self.pause = pause
		Entrez.email = email
		Entrez.tool = ENTREZ_TOOL
		if api_key:
			Entrez.api_key = api_key

	def map(self, proteins: list[str]) -> Resolution:
		result = Resolution()
		ipg = [p for p in proteins if not p.startswith("XP_")]
		for start in range(0, len(ipg), self.CHUNK):
			chunk = ipg[start:start + self.CHUNK]
			try:
				self._ipg_chunk(chunk, result)
			except Exception as error:
				for p in chunk:
					result.errors[p] = "IPG: {}: {}".format(type(error).__name__, error)
				debug("IPG chunk {}-{} failed".format(start, start + len(chunk)), exc=True)
		for protein in proteins:
			if protein.startswith("XP_"):
				try:
					self._xp(protein, result)
				except Exception as error:
					result.errors[protein] = "BioProject: {}: {}".format(type(error).__name__, error)
		self._rank(result, proteins)
		return result

	def _rank(self, result: Resolution, proteins: list[str]) -> None:
		gcf_first = lambda a: (0 if a.startswith("GCF") else 1, a)
		for protein in proteins:
			found = set(result.accessions.get(protein, {}))
			if not self.cross_db:
				want = "GCF" if REFSEQ_PROTEIN.match(protein) else "GCA"
				keep = {a for a in found if a.startswith(want)}
				if found - keep:
					result.dropped[protein] = found - keep
				found = keep
			result.assemblies[protein] = sorted(found, key=gcf_first)[:self.max_assemblies]

	def _fetch_ipg(self, proteins: list[str]) -> str:
		time.sleep(self.pause)
		handle = Entrez.efetch(db="ipg", id=",".join(proteins), rettype="ipg", retmode="text")
		data = handle.read()
		handle.close()
		return data.decode("utf-8", errors="replace") if isinstance(data, bytes) else data

	def _ipg_chunk(self, proteins: list[str], result: Resolution) -> None:
		queries = set(proteins)
		groups: dict[str, dict] = {}
		for line in self._fetch_ipg(proteins).splitlines():
			if line.startswith("Id") or not ASSEMBLY.search(line):
				continue
			fields = line.rstrip().split("\t")
			ipg_id, accession, assembly = fields[0], fields[6], fields[-1]
			group = groups.setdefault(ipg_id, {"queries": set(), "by_assembly": {}})
			group["by_assembly"].setdefault(assembly, set()).add(accession)
			if accession in queries:
				group["queries"].add(accession)
		for group in groups.values():
			for q in group["queries"]:
				entry = result.accessions.setdefault(q, {})
				for assembly, accessions in group["by_assembly"].items():
					entry.setdefault(assembly, set()).update(accessions)

	def _xp(self, accession: str, result: Resolution) -> None:
		time.sleep(self.pause)
		handle = Entrez.efetch(db="protein", id=accession, rettype="gbwithparts", retmode="text")
		record = SeqIO.read(handle, "genbank")
		handle.close()
		projects = [x.split(":", 1)[1] for x in record.dbxrefs if x.startswith("BioProject:")]
		entry = result.accessions.setdefault(accession, {})
		for project in projects:
			for assembly in self._project_assemblies(project):
				entry.setdefault(assembly, set()).add(accession)

	def _project_assemblies(self, project: str) -> set[str]:
		time.sleep(self.pause)
		ids = Entrez.read(Entrez.esearch(db="bioproject", term=project)).get("IdList", [])
		if not ids:
			return set()
		time.sleep(self.pause)
		links = Entrez.read(Entrez.elink(dbfrom="bioproject", db="assembly", id=",".join(ids)))
		asm_ids = [link["Id"] for linkset in links for db in linkset.get("LinkSetDb", [])
			for link in db.get("Link", [])]
		if not asm_ids:
			return set()
		time.sleep(self.pause)
		summary = Entrez.read(Entrez.esummary(db="assembly", id=",".join(asm_ids)))
		docs = summary["DocumentSummarySet"]["DocumentSummary"]
		return {d.get("AssemblyAccession", "") for d in docs if ASSEMBLY.match(d.get("AssemblyAccession", ""))}


class NcbiGenomes(Downloader):
	FTP_RATE = 10.0

	def __init__(self, directory: Path, rate: float, workers: int):
		super().__init__(self.FTP_RATE, workers)
		self.directory = Path(directory)
		self.workers = workers

	@staticmethod
	def partition(assembly: str) -> str:
		prefix, digits = assembly.split("_")
		digits = digits.split(".")[0]
		return "{}/{}/{}/{}/{}".format(BASE, prefix, digits[0:3], digits[3:6], digits[6:9])

	def listing(self, url: str) -> str:
		r = self.get(url + "/")
		r.raise_for_status()
		return r.text

	def versioned(self, assembly: str) -> Optional[str]:
		for _ in range(3):
			try:
				for name in re.findall(r'href="([^"/]+)/"', self.listing(self.partition(assembly))):
					if name.startswith(assembly):
						return name
				self.failures[assembly] = "not on the NCBI FTP site"
				return None
			except Exception as error:
				self.failures[assembly] = "{}: {}".format(type(error).__name__, error)
				time.sleep(0.5)
		return None

	def fetch(self, assembly: str, slots: list[str]) -> dict[str, Path]:
		name = self.versioned(assembly)
		if name is None:
			return {}
		base = "{}/{}/{}".format(self.partition(assembly), name, name)
		got = {}
		with ThreadPoolExecutor(max_workers=max(len(slots), 1)) as pool:
			jobs = {slot: pool.submit(self.stream, base + SUFFIX[slot], self.directory / (assembly + SUFFIX[slot]))
				for slot in slots}
			for slot, job in jobs.items():
				if job.result():
					got[slot] = self.directory / (assembly + SUFFIX[slot])
		return got

	def fetch_many(self, wanted: dict[str, list[str]], progress=None) -> dict[str, dict[str, Path]]:
		results = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			jobs = {pool.submit(self.fetch, assembly, slots): assembly for assembly, slots in wanted.items()}
			for done, job in enumerate(as_completed(jobs), 1):
				results[jobs[job]] = job.result()
				if progress:
					progress(done, len(jobs))
		return results
