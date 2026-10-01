import re
import subprocess
import time
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Optional
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import Entrez, SeqIO

from flags3.log import debug, record_command
from flags3.tools import Tool, brief

ACCESSION = re.compile(r"^[ANYXW]P_\d+(\.\d+)?$", re.I)
INSDC = re.compile(r"^[A-Z]{3}\d{5,7}(\.\d+)?$")
RESIDUES = set("ACDEFGHIKLMNPQRSTVWYBZXUO*-")
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
POLL_SECONDS = 60
DB_ALIASES = {
	"remote": {"refseq_select": "refseq_select_prot", "refseq_protein": "refseq_protein", "genbank": "nr", "nr": "nr", "swissprot": "swissprot"},
	"local": {"refseq_select": "refseq_select_protein", "refseq_protein": "refseq_protein", "genbank": "nr", "nr": "nr", "swissprot": "swissprot"},
}


class BlastError(RuntimeError):
	pass


@dataclass(frozen=True)
class Query:
	name: str
	sequence: str
	accession: Optional[str]

	@property
	def label(self) -> str:
		return self.accession or self.name


@dataclass(frozen=True)
class BlastHit:
	accession: str
	evalue: float
	bitscore: float
	description: str


def parse_query(lines, source: str) -> Query:
	lines = [ln.strip() for ln in lines if ln.strip()]
	if not lines:
		raise BlastError("{} is empty".format(source))
	if len(lines) == 1 and not lines[0].startswith(">"):
		token = lines[0].split()[0]
		if ACCESSION.match(token) or INSDC.match(token):
			return Query(token, "", token)
	name = "query"
	if lines[0].startswith(">"):
		name = lines[0][1:].split()[0] or "query"
		lines = lines[1:]
	sequence = "".join(lines).replace(" ", "").upper()
	if not sequence:
		raise BlastError("{} has a FASTA header but no sequence".format(source))
	bad = sorted(set(sequence) - RESIDUES)
	if bad:
		raise BlastError("{} is neither a protein accession nor a protein sequence (unexpected: {})".format(source, " ".join(bad)))
	return Query(name, sequence, None)


def accession_in(*texts) -> Optional[str]:
	for token in re.split(r"[|\s]+", " ".join(t or "" for t in texts)):
		if ACCESSION.match(token) or INSDC.match(token):
			return token
	return None


def fetch_sequence(accession: str) -> str:
	handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
	try:
		return str(SeqIO.read(handle, "fasta").seq)
	finally:
		handle.close()


def _info(text: str, key: str) -> Optional[str]:
	match = re.search(r"^\s*{}\s*=\s*(\S+)".format(key), text, re.M)
	return match.group(1) if match else None


def _human(seconds: float) -> str:
	seconds = int(seconds)
	return "{}s".format(seconds) if seconds < 60 else "{}m{:02d}s".format(seconds // 60, seconds % 60)


class RemoteBlast:
	def __init__(self, database: str, evalue: float, max_hits: int, email: str, max_wait: float, report=print):
		self.database = DB_ALIASES["remote"].get(database, database)
		self.evalue = evalue
		self.max_hits = max_hits
		self.email = email
		self.max_wait = max_wait
		self.report = report

	def search(self, fasta: str) -> list[BlastHit]:
		import Bio.Blast as Blast
		rid, rtoe = self._submit(fasta)
		self.report("NCBI job {}{}; follow it at {}?CMD=Get&RID={}".format(
			rid, ", NCBI estimates {}".format(_human(rtoe)) if rtoe else "", BLAST_URL, rid))
		stream = self._collect(rid, rtoe)
		if not stream.getvalue():
			return []
		record = Blast.read(stream)
		hits = []
		for hit in record:
			accession = accession_in(hit.target.id, hit.target.description)
			if not accession:
				continue
			hsp = hit[0] if len(hit) else None
			hits.append(BlastHit(accession, float(hsp.annotations.get("evalue", 0.0)) if hsp else 0.0,
				float(hsp.annotations.get("bit score", 0.0)) if hsp else 0.0, (hit.target.description or "").strip()))
		return hits

	def _post(self, params: dict):
		if self.email:
			params.setdefault("email", self.email)
		params.setdefault("tool", "flags3")
		request = Request(BLAST_URL, data=urlencode(params).encode(), headers={"User-Agent": "flags3"})
		return urlopen(request, timeout=120)

	def _submit(self, fasta: str) -> tuple[str, float]:
		params = {"CMD": "Put", "PROGRAM": "blastp", "DATABASE": self.database, "QUERY": fasta,
			"EXPECT": str(self.evalue), "HITLIST_SIZE": str(self.max_hits)}
		with self._post(params) as response:
			text = response.read().decode("utf-8", errors="replace")
		rid = _info(text, "RID")
		if not rid:
			raise BlastError("NCBI did not return a job id; the database name {!r} may be unknown to QBLAST".format(self.database))
		try:
			rtoe = float(_info(text, "RTOE") or 0)
		except ValueError:
			rtoe = 0.0
		return rid, rtoe

	def _collect(self, rid: str, rtoe: float) -> BytesIO:
		started = time.time()
		time.sleep(min(max(rtoe, 5), POLL_SECONDS))
		while True:
			with self._post({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}) as response:
				text = response.read().decode("utf-8", errors="replace")
			status = _info(text, "Status") or "UNKNOWN"
			waited = time.time() - started
			debug("blast: {} status={} after {}".format(rid, status, _human(waited)))
			if status == "READY":
				if _info(text, "ThereAreHits") == "no":
					return BytesIO(b"")
				break
			if status == "FAILED":
				raise BlastError("NCBI reports job {} failed; see {}?CMD=Get&RID={}".format(rid, BLAST_URL, rid))
			if status == "UNKNOWN":
				raise BlastError("NCBI no longer knows job {}; it expired or was rejected".format(rid))
			if waited > self.max_wait:
				raise BlastError("still queued at NCBI after {}; the job may yet finish at {}?CMD=Get&RID={}".format(_human(waited), BLAST_URL, rid))
			self.report("still waiting on NCBI, {} elapsed".format(_human(waited)))
			time.sleep(POLL_SECONDS)
		with self._post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}) as response:
			return BytesIO(response.read())


class LocalBlast:
	def __init__(self, tool: Tool, database: str, evalue: float, max_hits: int, threads: int, work: Path):
		self.tool = tool
		self.database = DB_ALIASES["local"].get(database, database)
		self.evalue = evalue
		self.max_hits = max_hits
		self.threads = threads
		self.work = work

	def search(self, fasta: str) -> list[BlastHit]:
		found, where = self.tool.locate()
		if not found:
			raise BlastError("blastp: {}; flags3 install core sets it up".format(where))
		self.work.mkdir(parents=True, exist_ok=True)
		query = self.work / "query.fasta"
		query.write_text(fasta)
		argv = self.tool.argv(db=self.database, evalue=self.evalue, hits=self.max_hits, **{"in": query})
		if self.threads:
			argv += ["-num_threads", str(self.threads)]
		debug("blast: " + " ".join(argv))
		done = subprocess.run(argv, capture_output=True, text=True, env=self.tool.environment())
		record_command(argv, done.returncode, stderr=done.stderr)
		if done.returncode != 0:
			raise BlastError("blastp exited {}: {}".format(done.returncode, brief(done.stderr)))
		hits = []
		for line in done.stdout.splitlines():
			parts = line.split("\t")
			if len(parts) < 3 or not parts[0].strip():
				continue
			try:
				hits.append(BlastHit(parts[0].strip(), float(parts[1]), float(parts[2]), parts[3].strip() if len(parts) > 3 else ""))
			except ValueError:
				continue
		return hits


def dedupe(hits: list[BlastHit]) -> list[BlastHit]:
	seen, out = set(), []
	for hit in hits:
		key = hit.accession.split(".")[0]
		if key not in seen:
			seen.add(key)
			out.append(hit)
	return out
