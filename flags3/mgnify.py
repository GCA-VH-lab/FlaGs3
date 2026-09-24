import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from flags3.log import debug
from flags3.net import Downloader

API = "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/{}"
API_DOWNLOADS = "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/{}/downloads"
ACCESSION = re.compile(r"^MGYG\d+$")
SUFFIX = {"gff": ".gff", "faa": ".faa", "genome": ".fna"}


def is_mgnify(assembly: str) -> bool:
	return bool(ACCESSION.match(assembly))


def _url_of(entry) -> Optional[str]:
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


def url_map(payload) -> dict[str, str]:
	if isinstance(payload, dict):
		records = payload.get("downloads")
		if records is None:
			records = payload.get("items", payload.get("data", []))
	else:
		records = payload or []
	urls = {}
	for entry in records if isinstance(records, list) else []:
		url = _url_of(entry)
		if url:
			name = os.path.basename(url.split("?", 1)[0].rstrip("/"))
			if name:
				urls[name] = url
	return urls


class MgnifyGenomes(Downloader):
	API_RATE = 5.0

	def __init__(self, directory: Path, rate: float, workers: int):
		super().__init__(self.API_RATE, workers)
		self.directory = Path(directory)
		self.workers = workers

	def lookup(self, assembly: str) -> dict[str, str]:
		clock = time.perf_counter()
		r = self.get(API.format(assembly))
		r.raise_for_status()
		urls = url_map(r.json())
		debug("mgnify {}: genome record in {:.1f} s, {} download urls".format(assembly, time.perf_counter() - clock, len(urls)))
		if not urls:
			clock = time.perf_counter()
			r = self.get(API_DOWNLOADS.format(assembly))
			r.raise_for_status()
			urls = url_map(r.json())
			debug("mgnify {}: downloads record in {:.1f} s, {} urls".format(assembly, time.perf_counter() - clock, len(urls)))
		return urls

	def fetch(self, assembly: str, slots: list[str]) -> dict[str, Path]:
		try:
			urls = self.lookup(assembly)
		except Exception as error:
			self.failures[assembly] = "{}: {}".format(type(error).__name__, error)
			debug("MGnify lookup failed for {}: {!r}".format(assembly, error))
			return {}
		jobs = {}
		for slot in slots:
			if slot not in SUFFIX:
				continue
			base = assembly + SUFFIX[slot]
			name = next((n for n in (base + ".gz", base) if n in urls), None)
			if name is None:
				self.failures["{}:{}".format(assembly, slot)] = "not offered by MGnify"
				continue
			jobs[slot] = (urls[name], self.directory / name)
		got = {}

		def timed(slot, url, local):
			clock = time.perf_counter()
			ok = self.stream(url, local)
			size = local.stat().st_size if ok and local.is_file() else 0
			debug("mgnify {} {}: {} from {} in {:.1f} s ({:,d} B)".format(
				assembly, slot, "ok" if ok else "FAILED", url.split("/")[2], time.perf_counter() - clock, size))
			return ok

		with ThreadPoolExecutor(max_workers=max(len(jobs), 1)) as pool:
			running = {slot: pool.submit(timed, slot, url, local) for slot, (url, local) in jobs.items()}
			for slot, job in running.items():
				if job.result():
					got[slot] = jobs[slot][1]
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
