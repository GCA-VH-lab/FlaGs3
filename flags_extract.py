import gzip
import re
from typing import Dict, List, NamedTuple, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


class FlankingGene(NamedTuple):
	accession: str
	strand: str
	start: int
	end: int
	product: str
	offset: int
	query: str
	is_rna: bool = False
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
	up_available: int
	down_available: int
	up_reached: int
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
		self.sequences: Dict[str, str] = {}
		self.query_sequences: Dict[str, str] = {}
		self.row_sequences: Dict[str, str] = {}
		self.rna_sequences: Dict[str, str] = {}
		self.rna_products: Dict[str, str] = {}
		self.species: Dict[str, str] = {}
		self.row_label: Dict[str, str] = {}
		self.ranges: Dict[str, RangeInfo] = {}
		self.drawn_accessions: Dict[str, set] = {}
		self.window_genes: Dict[str, list] = {}
		self.window_sequences: Dict[str, str] = {}
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
				out[key] = value
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
