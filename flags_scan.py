import gzip
import os
import threading
from typing import Dict, List, NamedTuple, Optional, Tuple


class ScanWindow(NamedTuple):
	contig: str
	start: int         # analysis range, 1-based inclusive
	end: int
	slice_start: int   # what is actually handed to the tool, margin included
	slice_end: int


def merge_windows(spans, margin: int = 0,
				  contig_lengths: Optional[Dict[str, int]] = None
				  ) -> List[ScanWindow]:
	contig_lengths = contig_lengths or {}
	by_contig: Dict[str, List[Tuple[int, int]]] = {}
	for contig, start, end in spans:
		by_contig.setdefault(contig, []).append((start, end))
	out: List[ScanWindow] = []
	for contig in sorted(by_contig):
		limit = contig_lengths.get(contig)
		merged: List[List[int]] = []
		for start, end in sorted(by_contig[contig]):
			lo = max(1, start - margin)
			hi = end + margin
			if limit:
				hi = min(hi, limit)
			if merged and lo <= merged[-1][3] + 1:
				merged[-1][0] = min(merged[-1][0], start)
				merged[-1][1] = max(merged[-1][1], end)
				merged[-1][2] = min(merged[-1][2], lo)
				merged[-1][3] = max(merged[-1][3], hi)
			else:
				merged.append([start, end, lo, hi])
		for start, end, lo, hi in merged:
			out.append(ScanWindow(contig, start, end, lo, hi))
	return out


MAX_BATCH_BASES = 200_000_000   # one tool invocation's worth of sequence


def _records(genome_path, windows):
	from Bio import SeqIO
	wanted = {}
	for w in windows or []:
		wanted.setdefault(w.contig, []).append(w)
	opener = gzip.open if genome_path.endswith(".gz") else open
	with opener(genome_path, "rt", encoding="utf-8", errors="replace") as fin:
		for record in SeqIO.parse(fin, "fasta"):
			if windows:
				for w in wanted.get(record.id, []):
					lo = max(1, w.slice_start)
					hi = min(len(record.seq), w.slice_end)
					if hi >= lo:
						yield str(record.seq)[lo - 1:hi], lo - 1, w._replace(
							slice_start=lo, slice_end=hi)
			else:
				span = len(record.seq)
				if span:
					yield str(record.seq), 0, ScanWindow(
						record.id, 1, span, 1, span)


def describe_window(assembly, window, length, queries=None):
	"""The record id has to stay a short token, since the tools echo it back and
	it is what the offset table is keyed on. Everything a reader needs goes in
	the description after it."""
	primary = (queries or [None])[0]
	fields = ["assembly={}".format(assembly),
			  "contig={}".format(window.contig),
			  "len={}".format(length),
			  "genomic={}-{}".format(window.slice_start, window.slice_end),
			  "analysed={}-{}".format(window.start, window.end)]
	if queries:
		fields.insert(2, "query={}".format(",".join(queries)))
	if primary and primary in (WINDOW_QUERY_POS or {}):
		q_start, q_end, strand = WINDOW_QUERY_POS[primary]
		if strand == "-":
			lo, hi = q_end - window.slice_end, q_end - window.slice_start
		else:
			lo, hi = window.slice_start - q_start, window.slice_end - q_start
		fields += ["start={}".format(lo), "end={}".format(hi),
				   "query_at={}".format(abs(q_start - window.slice_start) + 1)]
	return " ".join(fields)


WINDOW_QUERY_POS = {}     # query -> (start, end, strand), filled in by the caller
WINDOW_QUERIES = {}       # (assembly, contig, start, end) -> [query, ...]


def _queries_for(assembly, window):
	direct = WINDOW_QUERIES.get((assembly, window.contig, window.start, window.end))
	if direct:
		return direct
	found = []
	for (asm, contig, lo, hi), names in WINDOW_QUERIES.items():
		if asm == assembly and contig == window.contig \
				and lo >= window.start and hi <= window.end:
			found.extend(names)
	return sorted(set(found))


def write_batches(jobs, out_dir: str, max_bases: int = MAX_BATCH_BASES):
	"""jobs: [(assembly, genome_path, windows)]. Yields (fasta_path, offsets),
	splitting so one file never holds more than max_bases of sequence."""
	os.makedirs(out_dir, exist_ok=True)
	batch, offsets, bases, index, number = None, {}, 0, 0, 0
	path = ""

	def open_batch(n):
		p = os.path.join(out_dir, "batch{:03d}.fna".format(n))
		return p, open(p, "w")

	for assembly, genome_path, windows in jobs:
		for seq, offset, window in _records(genome_path, windows):
			if batch is None:
				path, batch = open_batch(number)
			name = "w{}".format(index)
			index += 1
			batch.write(">{} {}\n{}\n".format(
				name, describe_window(assembly, window, len(seq),
									  _queries_for(assembly, window)), seq))
			offsets[name] = (assembly, offset, window)
			bases += len(seq)
			if bases >= max_bases:
				batch.close()
				yield path, offsets
				batch, offsets, bases = None, {}, 0
				number += 1
	if batch is not None:
		batch.close()
		yield path, offsets
	elif not offsets and index == 0:
		return


_STORE = {}
_STORE_LOCK = threading.Lock()


def shared_batches(key: str, jobs, out_dir: str,
				   max_bases: int = MAX_BATCH_BASES):
	"""Cut the windows once and let every tool asking for the same span reuse
	the files. Sismis and geNomad both take nucleotide FASTA, and at batch scale
	a second copy is hundreds of gigabytes."""
	with _STORE_LOCK:
		if key in _STORE:
			return _STORE[key]
		result = list(write_batches(jobs, os.path.join(out_dir, key), max_bases))
		_STORE[key] = result
		return result


def place_batched(record_name: str, start: int, end: int, offsets):
	placed = offsets.get(record_name)
	if placed is None and record_name:
		placed = offsets.get(str(record_name).split()[0])
	if placed is None:
		return None
	assembly, offset, window = placed
	lo, hi = start + offset, end + offset
	coverage = "full" if window.start <= lo and hi <= window.end else "partial"
	return assembly, window.contig, lo, hi, coverage

