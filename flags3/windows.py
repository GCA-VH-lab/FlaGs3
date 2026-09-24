from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Iterator

from flags3 import fasta
from flags3.schema import MISSING, Row, Window, split_row_id

MAX_BATCH_BASES = 200_000_000


@dataclass(frozen=True)
class Cut(Row):
	FILE: ClassVar[str] = "cuts.tsv"
	record: str
	batch: str
	assembly: str
	contig: str
	cut_lo: int
	cut_hi: int
	scan_lo: int
	scan_hi: int
	rows: str

	def place(self, start: int, end: int) -> tuple[int, int, str]:
		lo, hi = start + self.cut_lo - 1, end + self.cut_lo - 1
		coverage = "full" if self.scan_lo <= lo and hi <= self.scan_hi else "partial"
		return lo, hi, coverage


def merge(windows: list[Window]) -> dict[str, list[Cut]]:
	by_contig: dict[tuple[str, str], list[Window]] = {}
	for w in windows:
		by_contig.setdefault((w.assembly, w.contig), []).append(w)
	out: dict[str, list[Cut]] = {}
	for (assembly, contig), group in sorted(by_contig.items()):
		runs: list[list] = []
		for w in sorted(group, key=lambda w: (w.cut_lo, w.cut_hi)):
			if runs and w.cut_lo <= runs[-1][1] + 1:
				run = runs[-1]
				run[1] = max(run[1], w.cut_hi)
				run[2] = min(run[2], w.scan_lo)
				run[3] = max(run[3], w.scan_hi)
				run[4].append(w.row_id)
			else:
				runs.append([w.cut_lo, w.cut_hi, w.scan_lo, w.scan_hi, [w.row_id]])
		for cut_lo, cut_hi, scan_lo, scan_hi, rows in runs:
			out.setdefault(assembly, []).append(Cut("", "", assembly, contig, cut_lo, cut_hi, scan_lo, scan_hi, ",".join(rows)))
	return out


def describe(cut: Cut, length: int) -> str:
	queries = ",".join(split_row_id(r)[0] for r in cut.rows.split(","))
	return "assembly={} contig={} query={} len={} genomic={}-{} analysed={}-{}".format(
		cut.assembly, cut.contig, queries, length, cut.cut_lo, cut.cut_hi, cut.scan_lo, cut.scan_hi)


class Batches:
	def __init__(self, directory: Path, max_bases: int = MAX_BATCH_BASES):
		self.directory = directory
		self.max_bases = max_bases
		self.cuts: list[Cut] = []
		self.paths: list[Path] = []

	def write(self, genomes: dict[str, Path], cuts: dict[str, list[Cut]]) -> list[Path]:
		self.directory.mkdir(parents=True, exist_ok=True)
		handle, bases, index = None, 0, 0
		for assembly in sorted(cuts):
			wanted: dict[str, list[Cut]] = {}
			for cut in cuts[assembly]:
				wanted.setdefault(cut.contig, []).append(cut)
			for contig, _, seq in fasta.read(genomes[assembly]):
				for cut in wanted.get(contig, []):
					lo, hi = max(1, cut.cut_lo), min(len(seq), cut.cut_hi)
					if hi < lo:
						continue
					if handle is None:
						handle = self._open(len(self.paths))
					record = "w{}".format(index)
					index += 1
					placed = Cut(record, self.paths[-1].name, assembly, contig, lo, hi, cut.scan_lo, cut.scan_hi, cut.rows)
					piece = seq[lo - 1:hi]
					handle.write(">{} {}\n".format(record, describe(placed, len(piece))))
					for at in range(0, len(piece), 60):
						handle.write(piece[at:at + 60] + "\n")
					self.cuts.append(placed)
					bases += len(piece)
					if bases >= self.max_bases:
						handle.close()
						handle, bases = None, 0
		if handle is not None:
			handle.close()
		Cut.write(self.directory / Cut.FILE, self.cuts)
		return list(self.paths)

	def _open(self, number: int):
		path = self.directory / "batch{:03d}.fna".format(number)
		self.paths.append(path)
		return open(path, "w")

	def by_record(self) -> dict[str, Cut]:
		return {c.record: c for c in self.cuts}

	def per_assembly(self) -> dict[str, tuple[int, int]]:
		stats: dict[str, list[int]] = {}
		for c in self.cuts:
			entry = stats.setdefault(c.assembly, [0, 0])
			entry[0] += 1
			entry[1] += c.cut_hi - c.cut_lo + 1
		return {a: (n, b) for a, (n, b) in stats.items()}


def overlapping_rows(windows: list[Window], assembly: str, contig: str, start: int, end: int) -> list[str]:
	return [w.row_id for w in windows if w.assembly == assembly and w.contig == contig and start <= w.hi and end >= w.lo]
