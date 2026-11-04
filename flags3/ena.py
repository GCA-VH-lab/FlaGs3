import gzip
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional
from urllib.parse import quote

from Bio import SeqIO

from flags3 import fasta
from flags3.log import debug
from flags3.net import Downloader

API = "https://www.ebi.ac.uk/ena/browser/api/embl/{}"
WGS = "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/{}/{}.dat.gz"
SEQUENCE = re.compile(r"^[A-Z]{1,2}\d{5,8}(\.\d+)?$")
CONTIG = re.compile(r"^(?!MGYG)[A-Z]{4,6}\d{8,11}(\.\d+)?$")
WGS_SET = re.compile(r"^([A-Z]{4,6}\d{2})(0{6,9})?$")
SUFFIX = {"gff": ".gff.gz", "faa": ".faa.gz", "rna": ".rna.fna.gz", "genome": ".fna.gz"}
RNA_TYPES = ("tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA")
GFF_SAFE = "-_.:^*$@!+/()[]{}' "
WIDTH = 60


def is_ena(assembly: str) -> bool:
	return bool(WGS_SET.match(assembly) or SEQUENCE.match(assembly) or CONTIG.match(assembly))


def wgs_set(assembly: str) -> Optional[str]:
	match = WGS_SET.match(assembly)
	return match.group(1) if match else None


def cache_name(assembly: str) -> str:
	return wgs_set(assembly) or assembly


def url_of(assembly: str) -> str:
	prefix = wgs_set(assembly)
	if prefix:
		return WGS.format(prefix[:3].lower(), prefix)
	return API.format(assembly)


class EmblRecords:
	def __init__(self, flat: Path):
		self.flat = Path(flat)
		self.records = 0
		self.counts = {"faa": 0, "rna": 0}

	def write(self, directory: Path, name: str) -> dict[str, Path]:
		paths = {slot: directory / (name + suffix) for slot, suffix in SUFFIX.items()}
		parts = {slot: path.with_name(path.name + ".part") for slot, path in paths.items()}
		handles = {slot: gzip.open(part, "wt", encoding="utf-8") for slot, part in parts.items()}
		orphans = []
		try:
			handles["gff"].write("##gff-version 3\n")
			with fasta.open_text(self.flat) as source:
				for record in SeqIO.parse(source, "embl"):
					orphans.extend(self._record(record, handles))
			for line in orphans:
				handles["gff"].write(line)
		finally:
			for handle in handles.values():
				handle.close()
		kept = {slot for slot in paths if self.records > 0 and (slot != "faa" or self.counts["faa"] > 0)}
		got = {}
		for slot, part in parts.items():
			if slot in kept:
				os.replace(part, paths[slot])
				got[slot] = paths[slot]
			else:
				part.unlink(missing_ok=True)
		return got

	def _record(self, record, handles) -> list[str]:
		self.records += 1
		gff, contig, length = handles["gff"], record.id, len(record.seq)
		organism = record.annotations.get("organism", "")
		attrs = {"ID": contig}
		if record.annotations.get("topology", "") == "circular":
			attrs["Is_circular"] = "true"
		gff.write("##sequence-region {} 1 {}\n".format(contig, length))
		gff.write(self._line(contig, "region", 1, length, "+", attrs))
		children: list[tuple[str, int, int]] = []
		pending = []
		for feature in record.features:
			if feature.type == "CDS":
				self._cds(feature, contig, organism, handles, children)
			elif feature.type in RNA_TYPES:
				self._rna(feature, contig, record, handles, children)
			elif feature.type == "gene" and _pseudo(feature):
				pending.append(feature)
		orphans = []
		for feature in pending:
			start, end = _span(feature)
			locus = _first(feature, "locus_tag")
			claimed = any(c[0] == locus for c in children) if locus else any(c[1:] == (start, end) for c in children)
			if not claimed:
				orphans.append(self._line(contig, "gene", start, end, _strand(feature),
					{"ID": "gene-" + (locus or "{}:{}".format(start, end)), "locus_tag": locus, "gene_biotype": "pseudogene"}))
		self._fasta(handles["genome"], contig, record.description, str(record.seq))
		return orphans

	def _cds(self, feature, contig, organism, handles, children) -> None:
		start, end = _span(feature)
		locus = _first(feature, "locus_tag")
		protein = _first(feature, "protein_id")
		product = _first(feature, "product")
		pseudo = _pseudo(feature)
		accession = protein or locus
		children.append((locus, start, end))
		attrs = {"ID": "cds-" + (accession or "{}:{}".format(start, end)), "protein_id": protein,
			"locus_tag": locus, "product": product}
		if pseudo:
			attrs["pseudo"] = "true"
		handles["gff"].write(self._line(contig, "CDS", start, end, _strand(feature), attrs))
		translation = _first(feature, "translation")
		if accession and translation and not pseudo:
			description = product + (" [{}]".format(organism) if organism else "")
			self._fasta(handles["faa"], accession, description.strip(), translation)
			self.counts["faa"] += 1

	def _rna(self, feature, contig, record, handles, children) -> None:
		start, end = _span(feature)
		locus = _first(feature, "locus_tag")
		product = _first(feature, "product")
		name = locus or "{}_{}_{}".format(feature.type, start, end)
		children.append((locus, start, end))
		handles["gff"].write(self._line(contig, feature.type, start, end, _strand(feature),
			{"ID": "rna-" + name, "locus_tag": locus, "product": product}))
		seq = str(feature.extract(record.seq))
		if seq:
			self._fasta(handles["rna"], name, "{} [locus_tag={}]".format(product, locus).strip(), seq)
			self.counts["rna"] += 1

	@staticmethod
	def _fasta(handle, name: str, description: str, seq: str) -> None:
		handle.write((">{} {}".format(name, description)).rstrip() + "\n")
		for i in range(0, len(seq), WIDTH):
			handle.write(seq[i:i + WIDTH] + "\n")

	@staticmethod
	def _line(contig, kind, start, end, strand, attrs) -> str:
		text = ";".join("{}={}".format(k, quote(v, safe=GFF_SAFE)) for k, v in attrs.items() if v)
		return "\t".join((contig, "ENA", kind, str(start), str(end), ".", strand, ".", text)) + "\n"


def _first(feature, key: str) -> str:
	values = feature.qualifiers.get(key)
	return str(values[0]).strip() if values else ""


def _pseudo(feature) -> bool:
	return "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers


def _span(feature) -> tuple[int, int]:
	return int(feature.location.start) + 1, int(feature.location.end)


def _strand(feature) -> str:
	return "-" if feature.location.strand == -1 else "+"


class EnaGenomes(Downloader):
	API_RATE = 5.0

	def __init__(self, directory: Path, rate: float, workers: int):
		super().__init__(self.API_RATE, workers)
		self.directory = Path(directory)
		self.workers = workers

	def fetch(self, assembly: str, slots: list[str]) -> dict[str, Path]:
		name = cache_name(assembly)
		url = url_of(assembly)
		flat = self.directory / (name + (".dat.gz" if url.endswith(".gz") else ".embl"))
		clock = time.perf_counter()
		if not self.stream(url, flat):
			self.failures[assembly] = self.failures.pop(url, "download failed")
			return {}
		debug("ena {}: flat file in {:.1f} s ({:,d} B)".format(assembly, time.perf_counter() - clock, flat.stat().st_size))
		try:
			got = EmblRecords(flat).write(self.directory, name)
		except Exception as error:
			self.failures[assembly] = "EMBL flat file: {}: {}".format(type(error).__name__, error)
			debug("ENA conversion failed for {}".format(assembly), exc=True)
			return {}
		finally:
			flat.unlink(missing_ok=True)
		if not got:
			self.failures[assembly] = "no sequence record in the ENA response"
		elif "faa" not in got:
			self.failures[assembly] = "the ENA record carries no annotated CDS with a translation"
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
