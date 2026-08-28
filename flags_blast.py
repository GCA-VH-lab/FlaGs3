import os
import re
import time
import shutil
import subprocess
import tempfile
from typing import List, NamedTuple, Optional

from Bio import Entrez, SeqIO

ACCESSION_RE = re.compile(r"^[ANYXW]P_\d+(\.\d+)?$", re.I)
INSDC_RE = re.compile(r"^[A-Z]{3}\d{5,7}(\.\d+)?$")
RESIDUES = set("ACDEFGHIKLMNPQRSTVWYBZXUO*-")

DB_ALIASES = {
	"remote": {"refseq_select": "refseq_select_prot",
			   "refseq_protein": "refseq_protein",
			   "genbank": "nr", "nr": "nr", "swissprot": "swissprot"},
	"local": {"refseq_select": "refseq_select_protein",
			  "refseq_protein": "refseq_protein",
			  "genbank": "nr", "nr": "nr", "swissprot": "swissprot"},
}


class BlastHit(NamedTuple):
	accession: str
	evalue: float
	bitscore: float
	description: str


class BlastQuery(NamedTuple):
	name: str
	sequence: str
	accession: Optional[str]


def _debug(message, exc=False):
	try:
		from flags_log import debug
		debug(message, exc=exc)
	except ImportError:
		pass


def read_query(path: str) -> BlastQuery:
	with open(path) as fh:
		lines = [ln.strip() for ln in fh if ln.strip()]
	if not lines:
		raise ValueError("--blast_input file is empty: {}".format(path))

	if len(lines) == 1 and not lines[0].startswith(">"):
		token = lines[0].split()[0]
		if ACCESSION_RE.match(token) or INSDC_RE.match(token):
			return BlastQuery(name=token, sequence="", accession=token)

	name = "query"
	if lines[0].startswith(">"):
		name = lines[0][1:].split()[0] or "query"
		lines = lines[1:]
	sequence = "".join(lines).replace(" ", "").upper()
	if not sequence:
		raise ValueError("--blast_input has a FASTA header but no sequence")
	bad = sorted(set(sequence) - RESIDUES)
	if bad:
		raise ValueError(
			"--blast_input is neither a RefSeq accession nor a protein sequence "
			"(unexpected characters: {})".format(" ".join(bad)))
	return BlastQuery(name=name, sequence=sequence, accession=None)


def fetch_sequence(accession: str) -> str:
	"""Resolve an accession to residues so both modes share one code path."""
	handle = Entrez.efetch(db="protein", id=accession, rettype="fasta",
						   retmode="text")
	try:
		record = SeqIO.read(handle, "fasta")
	finally:
		handle.close()
	return str(record.seq)


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
POLL_SECONDS = 60
MAX_WAIT_SECONDS = 3600


def _qblast_info(text: str, key: str):
	m = re.search(r"^\s*{}\s*=\s*(\S+)".format(key), text, re.M)
	return m.group(1) if m else None


def _human(seconds: float) -> str:
	seconds = int(seconds)
	if seconds < 60:
		return "{}s".format(seconds)
	return "{}m{:02d}s".format(seconds // 60, seconds % 60)


class BlastSearcher:

	def __init__(self, mode: str = "remote", database: str = "refseq_select",
				 evalue: float = 1e-5, max_hits: int = 50, threads: int = 0,
				 report=None, email: str = "", max_wait: float = MAX_WAIT_SECONDS):
		self.report = report or (lambda msg: None)
		self.email = email
		self.max_wait = max_wait
		self.mode = mode
		self.database = DB_ALIASES.get(mode, {}).get(database, database)
		self.evalue = evalue
		self.max_hits = max_hits
		self.threads = threads

	def search(self, query: BlastQuery) -> List[BlastHit]:
		sequence = query.sequence
		if not sequence:
			sequence = fetch_sequence(query.accession)
			_debug("blast: fetched {} residues for {}".format(
				len(sequence), query.accession))
		fasta = ">{}\n{}\n".format(query.name, sequence)
		if self.mode == "local":
			hits = self._local(fasta)
		else:
			hits = self._remote(fasta)
		return self._dedupe(hits, query.accession)

	def _remote(self, fasta: str) -> List[BlastHit]:
		import Bio.Blast as Blast
		_debug("blast: remote qblast against {} (evalue {}, {} hits)".format(
			self.database, self.evalue, self.max_hits))
		rid, rtoe = self._submit(fasta)
		note = "NCBI job {}".format(rid)
		if rtoe:
			note += ", NCBI estimates {}".format(_human(rtoe))
		self.report(note)
		self.report("follow it at {}?CMD=Get&RID={}".format(BLAST_URL, rid))
		stream = self._collect(rid, rtoe)
		try:
			record = Blast.read(stream)
		finally:
			stream.close()
		hits = []
		for hit in record:
			accession = self._accession_of(hit.target.id, hit.target.description)
			if not accession:
				continue
			hsp = hit[0] if len(hit) else None
			hits.append(BlastHit(
				accession=accession,
				evalue=float(hsp.annotations.get("evalue", 0.0)) if hsp else 0.0,
				bitscore=float(hsp.annotations.get("bit score", 0.0)) if hsp else 0.0,
				description=(hit.target.description or "").strip()))
		return hits

	def _post(self, params: dict):
		from urllib.parse import urlencode
		from urllib.request import Request, urlopen
		if self.email:
			params.setdefault("email", self.email)
		params.setdefault("tool", "FlaGs3")
		req = Request(BLAST_URL, data=urlencode(params).encode(),
					  headers={"User-Agent": "BiopythonClient"})
		return urlopen(req, timeout=120)

	def _submit(self, fasta: str):
		params = {"CMD": "Put", "PROGRAM": "blastp",
				  "DATABASE": self.database, "QUERY": fasta,
				  "EXPECT": str(self.evalue),
				  "HITLIST_SIZE": str(self.max_hits)}
		_debug("blast: Put {}".format(
			{k: v for k, v in params.items() if k != "QUERY"}))
		with self._post(params) as resp:
			text = resp.read().decode("utf-8", errors="replace")
		rid = _qblast_info(text, "RID")
		if not rid:
			raise RuntimeError(
				"NCBI did not return a job id; the database name {!r} may be "
				"unknown to QBLAST".format(self.database))
		try:
			rtoe = float(_qblast_info(text, "RTOE") or 0)
		except ValueError:
			rtoe = 0.0
		return rid, rtoe

	def _collect(self, rid: str, rtoe: float):
		from io import BytesIO
		started = time.time()
		time.sleep(min(max(rtoe, 5), POLL_SECONDS))
		while True:
			with self._post({"CMD": "Get", "RID": rid,
							 "FORMAT_OBJECT": "SearchInfo"}) as resp:
				text = resp.read().decode("utf-8", errors="replace")
			status = _qblast_info(text, "Status") or "UNKNOWN"
			waited = time.time() - started
			_debug("blast: {} status={} after {}".format(rid, status, _human(waited)))
			if status == "READY":
				if _qblast_info(text, "ThereAreHits") == "no":
					return BytesIO(b"")
				break
			if status == "FAILED":
				raise RuntimeError(
					"NCBI reports job {} failed. Check {}?CMD=Get&RID={}".format(
						rid, BLAST_URL, rid))
			if status == "UNKNOWN":
				raise RuntimeError(
					"NCBI no longer knows job {}; it expired or was rejected".format(rid))
			if waited > self.max_wait:
				raise TimeoutError(
					"still queued at NCBI after {}. The job may yet finish -- see "
					"{}?CMD=Get&RID={}".format(_human(waited), BLAST_URL, rid))
			self.report("still waiting on NCBI, {} elapsed".format(_human(waited)))
			time.sleep(POLL_SECONDS)
		with self._post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}) as resp:
			return BytesIO(resp.read())

	def _local(self, fasta: str) -> List[BlastHit]:
		binary = shutil.which("blastp")
		if binary is None:
			raise FileNotFoundError(
				"--blast_mode local needs blastp on PATH (NCBI BLAST+)")
		tmp = tempfile.mkdtemp(prefix="flags_blast_")
		try:
			query_path = os.path.join(tmp, "query.fasta")
			with open(query_path, "w") as out:
				out.write(fasta)
			cmd = [binary, "-query", query_path, "-db", self.database,
				   "-outfmt", "6 sacc evalue bitscore stitle",
				   "-evalue", str(self.evalue),
				   "-max_target_seqs", str(self.max_hits)]
			if self.threads:
				cmd += ["-num_threads", str(self.threads)]
			_debug("blast: running {}".format(" ".join(cmd)))
			result = subprocess.run(cmd, capture_output=True, text=True)
			_debug("blast: exit {}".format(result.returncode))
			if result.returncode != 0:
				raise RuntimeError("blastp failed ({}): {}".format(
					result.returncode, (result.stderr or "").strip()[:300]))
			hits = []
			for line in result.stdout.splitlines():
				parts = line.split("\t")
				if len(parts) < 3:
					continue
				accession = parts[0].strip()
				if not accession:
					continue
				try:
					evalue, bitscore = float(parts[1]), float(parts[2])
				except ValueError:
					continue
				hits.append(BlastHit(accession, evalue, bitscore,
									 parts[3].strip() if len(parts) > 3 else ""))
			return hits
		finally:
			shutil.rmtree(tmp, ignore_errors=True)

	@staticmethod
	def _accession_of(target_id: str, description: str) -> Optional[str]:
		"""Pull an accession out of 'ref|WP_000000001.1|' or a bare id."""
		for token in re.split(r"[|\s]+", "{} {}".format(target_id or "",
														description or "")):
			token = token.strip()
			if ACCESSION_RE.match(token) or INSDC_RE.match(token):
				return token
		return None

	@staticmethod
	def _dedupe(hits: List[BlastHit], self_accession: Optional[str]) -> List[BlastHit]:
		"""Keep best-first order, one row per accession, ignoring the version."""
		seen = set()
		out = []
		for hit in hits:
			key = hit.accession.split(".")[0]
			if key in seen:
				continue
			seen.add(key)
			out.append(hit)
		return out


def write_report(hits: List[BlastHit], query: BlastQuery, path: str):
	with open(path, "w") as out:
		out.write("#query\t{}\n".format(query.accession or query.name))
		out.write("#accession\tevalue\tbitscore\tdescription\n")
		for hit in hits:
			out.write("{}\t{:.3g}\t{:.1f}\t{}\n".format(
				hit.accession, hit.evalue, hit.bitscore, hit.description))
