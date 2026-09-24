from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from flags3 import domains, fasta, home
from flags3.log import note
from flags3.schema import MISSING, Annotation, Row
from flags3.stage import Stage


@dataclass(frozen=True)
class Domain(Row):
	FILE: ClassVar[str] = "domains.tsv"
	protein: str
	database: str
	domain: str
	group: str
	pfam: str
	clan: str
	start: int
	end: int
	evalue: float
	interpro: str
	interpro_name: str
	interpro_type: str
	characterization: str
	informativeness: str
	interpretation: str


class Domains(Stage):
	name = "domains"
	requires = ("extract",)
	optional = True

	def wanted(self, config) -> bool:
		return config.flag("domains")

	def run(self, run, config, out: Path) -> None:
		specs = [s for s in config.text("hmmdb").split(",") if s] or [str(home.PFAM_HMM)]
		coverage = domains.parse_coverage([s for s in config.text("hmm_coverage").split(",") if s])
		sources = [domains.HmmSource.parse(spec, coverage) for spec in specs]
		for source in sources:
			source.files()
		sequences = {name: seq for name, _, seq in fasta.read(run.stage_file("extract", "proteins.faa"))}
		note("scanning {} proteins against {}".format(len(sequences), ", ".join(s.name for s in sources)))
		scanner = domains.DomainScanner(sources, evalue=config.number("ethreshold", 1e-3), cpus=config.workers())
		hits = scanner.scan(sequences)
		for name, count in scanner.counts.items():
			note("{}: {} domain hits".format(name, count))
		clans = domains.load_clans(config.path("clans")) if config.path("clans") else {}
		interpro = domains.InterPro()
		table = domains.find_interpro(config.path("interpro"))
		if table:
			note("InterPro: {} Pfam entries from {}".format(interpro.load(table), table.name))
		hits.sort(key=lambda h: (h.protein, h.start, h.end))
		Domain.write(out / Domain.FILE, (self._row(h, clans, interpro) for h in hits))
		Annotation.write(out / Annotation.FILE, (
			Annotation(h.protein, "aa", h.start, h.end, "wedge",
				clans.get(h.accession) or clans.get(h.name) or h.group, h.name, h.database, h.evalue)
			for h in hits))
		note("{} domain hits on {} proteins".format(len(hits), len({h.protein for h in hits})))

	@staticmethod
	def _row(h: domains.DomainHit, clans: dict, interpro: domains.InterPro) -> Domain:
		meta = interpro.get(h.accession) if h.accession else {}
		pick = lambda key: meta.get(key) or MISSING
		return Domain(h.protein, h.database, h.name, h.group, h.accession or MISSING,
			clans.get(h.accession) or clans.get(h.name) or MISSING, h.start, h.end, h.evalue,
			pick("accession"), pick("name"), pick("type"), pick("characterization_status"),
			pick("informativeness_label"), pick("interpretation"))
