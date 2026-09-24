from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from Bio import Entrez

from flags3 import blast
from flags3.inputs import InputList
from flags3.log import note
from flags3.schema import Row
from flags3.stage import Stage
from flags3.tools import Tools

ACCESSIONS = "accessions.txt"


@dataclass(frozen=True)
class HitRow(Row):
	FILE: ClassVar[str] = "hits.tsv"
	query: str
	accession: str
	evalue: float
	bitscore: float
	description: str


def queries_of(run, config) -> list[blast.Query]:
	queries = []
	source = config.path("blast_input")
	if source:
		queries.append(blast.parse_query(source.read_text().splitlines(), str(source)))
	inputs = InputList()
	for name in config.text("inputs").split(","):
		if name:
			inputs.read(run.input_dir / name)
	for text, path, line in inputs.blast:
		queries.append(blast.parse_query([text], "line {} of {}".format(line, path)))
	return queries


class Blast(Stage):
	name = "blast"

	def wanted(self, config) -> bool:
		return bool(config.text("blast_input")) or config.flag("blast_inline")

	def run(self, run, config, out: Path) -> None:
		queries = queries_of(run, config)
		if not queries:
			raise blast.BlastError("no BLAST query: give -bi, or mark a line of the input list with a tab and BLAST")
		if config.integer("blast_hits", 50) < 2:
			raise blast.BlastError("--blast_hits must be at least 2")
		mode = config.text("blast_mode", "remote")
		email = config.text("user_email")
		if mode == "remote" and not email:
			raise blast.BlastError("remote BLAST needs -u/--user_email")
		Entrez.email = email
		Entrez.tool = "flags3"
		if config.text("api_key"):
			Entrez.api_key = config.text("api_key")
		searcher = self._searcher(mode, config, out)
		rows, accessions, seen = [], [], set()
		for query in queries:
			sequence = query.sequence or blast.fetch_sequence(query.accession)
			note("BlastP ({}) for {} against {}".format(mode, query.label, config.text("blast_db", "refseq_select")))
			if mode == "remote":
				print(">> waiting on NCBI's queue; this often takes several minutes.", flush=True)
			hits = blast.dedupe(searcher.search(">{}\n{}\n".format(query.name, sequence)))
			if not hits:
				print("Warning: BlastP returned no hits for {}. Try a larger --blast_evalue or a fuller --blast_db.".format(query.label))
			for hit in hits:
				rows.append(HitRow(query.label, hit.accession, hit.evalue, hit.bitscore, hit.description or "-"))
				key = hit.accession.split(".")[0]
				if key not in seen:
					seen.add(key)
					accessions.append(hit.accession)
			note("BlastP: {} hits for {}".format(len(hits), query.label))
		HitRow.write(out / HitRow.FILE, rows)
		with open(out / ACCESSIONS, "w") as handle:
			handle.write("# BlastP hits from {}\n# Reuse with -i to repeat the run without searching again.\n".format(
				", ".join(q.label for q in queries)))
			for accession in accessions:
				handle.write(accession + "\n")
		if not accessions:
			raise blast.BlastError("BlastP returned no hits for any query")

	@staticmethod
	def _searcher(mode: str, config, out: Path):
		database = config.text("blast_db", "refseq_select")
		evalue = config.number("blast_evalue", 1e-5)
		hits = config.integer("blast_hits", 50)
		if mode == "local":
			tools = Tools.load(config.path("tools"))
			return blast.LocalBlast(tools["blastp"], database, evalue, hits, config.workers(), out / "raw")
		return blast.RemoteBlast(database, evalue, hits, config.text("user_email"),
			config.number("blast_wait", 60) * 60, report=lambda m: print(">> " + m, flush=True))
