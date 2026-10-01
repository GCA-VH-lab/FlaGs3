import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from flags3 import fasta

PSEUDOGENE = "pseudogene*"


@dataclass
class Feature:
	contig: str
	start: int
	end: int
	strand: str
	accession: Optional[str]
	product: str
	biotype: str
	locus_tag: str
	is_rna: bool


@dataclass
class Contig:
	name: str
	length: int
	first: int
	last: int
	reach: list[int] = field(default_factory=list)


class GeneTable:
	def __init__(self, gff: Path):
		self.genes: list[Feature] = []
		self.contigs: dict[str, Contig] = {}
		self._by_accession: dict[str, int] = {}
		self._parse(gff)
		self._index()

	def find(self, accessions: set[str]) -> Optional[int]:
		hits = [self._by_accession[a] for a in accessions if a in self._by_accession]
		return min(hits) if hits else None

	def contig_of(self, index: int) -> Contig:
		return self.contigs[self.genes[index].contig]

	def _index(self):
		self.genes.sort(key=lambda g: (g.contig, g.start))
		lengths = self._lengths
		for i, g in enumerate(self.genes):
			contig = self.contigs.get(g.contig)
			if contig is None:
				contig = Contig(g.contig, lengths.get(g.contig, 0), i, i + 1)
				self.contigs[g.contig] = contig
			contig.last = i + 1
			contig.reach.append(max(contig.reach[-1] if contig.reach else 0, g.end))
			self._by_accession.setdefault(g.accession, i)
		for contig in self.contigs.values():
			if not contig.length:
				contig.length = contig.reach[-1]

	def _parse(self, gff: Path):
		self._lengths = {}
		with fasta.open_text(gff) as handle:
			for raw in handle:
				if raw.startswith("#"):
					self._header(raw)
					continue
				col = raw.rstrip("\n").split("\t")
				if len(col) < 9:
					continue
				feature, attrs = col[2], _attrs(col[8])
				if feature == "region":
					self._lengths.setdefault(col[0], _int(col[4]))
				elif feature.endswith("gene"):
					self.genes.append(Feature(col[0], int(col[3]), int(col[4]), col[6], None, "",
						attrs.get("gene_biotype", ""), attrs.get("locus_tag", ""), False))
				elif feature == "CDS":
					self._cds(col, attrs)
				elif feature.endswith("RNA"):
					self._rna(col, feature, attrs)
		for g in self.genes:
			if g.biotype == "pseudogene":
				g.accession = PSEUDOGENE
			elif not g.accession:
				g.accession = (g.biotype or "noProtein") + "*"

	def _header(self, raw: str):
		if raw.startswith("##sequence-region"):
			parts = raw.split()
			if len(parts) >= 4:
				self._lengths[parts[1]] = _int(parts[3])

	def _cds(self, col, attrs):
		locus_tag = attrs.get("locus_tag", "")
		accession = (attrs.get("protein_id") or _strip_prefix(attrs.get("ID"))
			or locus_tag or attrs.get("Name"))
		product = attrs.get("product", "")
		pseudo = attrs.get("pseudo", "").lower() == "true"
		last = self.genes[-1] if self.genes else None
		if last is not None and last.accession is None:
			if pseudo or last.biotype == "pseudogene":
				last.biotype = "pseudogene"
			else:
				last.accession = accession
			last.product = product
		else:
			self.genes.append(Feature(col[0], int(col[3]), int(col[4]), col[6], accession, product,
				"pseudogene" if pseudo else "protein_coding", locus_tag, False))

	def _rna(self, col, feature, attrs):
		locus_tag = attrs.get("locus_tag", "")
		accession = (attrs.get("Name") or attrs.get("transcript_id")
			or _strip_prefix(attrs.get("ID")) or locus_tag)
		product = attrs.get("product", "")
		last = self.genes[-1] if self.genes else None
		if last is not None and last.accession is None:
			last.accession = accession
			last.product = product
			last.is_rna = True
			last.locus_tag = last.locus_tag or locus_tag
		else:
			self.genes.append(Feature(col[0], int(col[3]), int(col[4]), col[6], accession, product,
				feature, locus_tag, True))


class ProteinFasta:
	def __init__(self, path: Path):
		self.sequences: dict[str, str] = {}
		self.organisms: dict[str, str] = {}
		for name, description, seq in fasta.read(path):
			self.sequences[name] = seq
			match = re.search(r"\[([^\]]+)\]\s*$", description)
			if match:
				self.organisms[name] = match.group(1)

	def get(self, accession: str) -> Optional[str]:
		return self.sequences.get(accession)


class RnaFasta:
	def __init__(self, path: Path):
		self.sequences: dict[str, str] = {}
		for name, description, seq in fasta.read(path):
			self.sequences[name] = seq
			match = re.search(r"\[locus_tag=([^\]]+)\]", description)
			if match:
				self.sequences[match.group(1)] = seq

	def get(self, gene: Feature) -> Optional[str]:
		return self.sequences.get(gene.locus_tag) or self.sequences.get(gene.accession)


class GenomeFasta:
	def __init__(self, path: Path):
		self.contigs = {name: seq for name, _, seq in fasta.read(path)}

	def slice(self, gene: Feature) -> Optional[str]:
		seq = self.contigs.get(gene.contig)
		if seq is None:
			return None
		sub = seq[gene.start - 1:gene.end]
		if gene.strand == "-":
			sub = fasta.reverse_complement(sub)
		return sub or None


def _attrs(text: str) -> dict[str, str]:
	out: dict[str, str] = {}
	for item in text.split(";"):
		key, sep, value = item.partition("=")
		if sep and key not in out:
			out[key] = value
	return out


def _strip_prefix(value: Optional[str]) -> Optional[str]:
	if not value:
		return None
	for prefix in ("cds-", "rna-"):
		if value.startswith(prefix):
			return value[len(prefix):]
	return value


def _int(text: str) -> int:
	try:
		return int(text)
	except ValueError:
		return 0
