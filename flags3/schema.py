import csv
from dataclasses import dataclass, fields
from pathlib import Path
from typing import ClassVar, Iterable, Iterator, Optional, Type, TypeVar

MISSING = "-"
KINDS = frozenset(("fill", "outline", "band", "wedge", "segment", "triangle"))
SPACES = frozenset(("aa", "bp", MISSING))
STRANDS = frozenset(("+", "-"))

R = TypeVar("R", bound="Row")


class SchemaError(ValueError):
	pass


def bp_subject(assembly: str, contig: str) -> str:
	return "{}|{}".format(assembly, contig)


def split_bp_subject(subject: str) -> tuple[str, str]:
	assembly, sep, contig = subject.rpartition("|")
	if not sep:
		raise SchemaError("not a bp subject: {}".format(subject))
	return assembly, contig


def row_id(query: str, assembly: str) -> str:
	return "{}|{}".format(query, assembly)


def split_row_id(row: str) -> tuple[str, str]:
	query, sep, assembly = row.rpartition("|")
	if not sep:
		raise SchemaError("not a row id: {}".format(row))
	return query, assembly


class Row:
	FILE: ClassVar[str] = ""

	@classmethod
	def columns(cls) -> tuple[str, ...]:
		return tuple(f.name for f in fields(cls))

	@classmethod
	def from_strings(cls: Type[R], values: dict[str, str]) -> R:
		kwargs = {}
		for f in fields(cls):
			raw = values.get(f.name, MISSING)
			kwargs[f.name] = _parse(raw, f.type)
		return cls(**kwargs)

	def to_strings(self) -> list[str]:
		return [_format(getattr(self, f.name)) for f in fields(self)]

	@classmethod
	def read(cls: Type[R], path: Path) -> list[R]:
		return list(cls.iterate(path))

	@classmethod
	def iterate(cls: Type[R], path: Path) -> Iterator[R]:
		with open(path, newline="", encoding="utf-8") as handle:
			reader = csv.DictReader(handle, delimiter="\t")
			missing = set(cls.columns()) - set(reader.fieldnames or ())
			if missing:
				raise SchemaError("{} lacks columns: {}".format(
					path, ", ".join(sorted(missing))))
			for line in reader:
				yield cls.from_strings(line)

	@classmethod
	def write(cls, path: Path, rows: Iterable["Row"]) -> int:
		count = 0
		with open(path, "w", newline="", encoding="utf-8") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(cls.columns())
			for row in rows:
				writer.writerow(row.to_strings())
				count += 1
		return count


def _parse(raw: str, kind):
	if raw == MISSING or raw == "":
		if kind is str:
			return MISSING
		return None
	if kind is str:
		return raw
	if kind is bool or kind == Optional[bool]:
		return raw.lower() in ("true", "1", "yes")
	if kind is int or kind == Optional[int]:
		return int(raw)
	if kind is float or kind == Optional[float]:
		return float(raw)
	raise SchemaError("no parser for column type {}".format(kind))


def _format(value) -> str:
	if value is None:
		return MISSING
	if isinstance(value, bool):
		return "true" if value else "false"
	if isinstance(value, float):
		return format(value, ".4g")
	return str(value)


@dataclass(frozen=True)
class Gene(Row):
	FILE: ClassVar[str] = "genes.tsv"
	row_id: str
	assembly: str
	contig: str
	accession: str
	start: int
	end: int
	strand: str
	product: str
	is_rna: bool
	offset: int

	def __post_init__(self):
		if self.strand not in STRANDS:
			raise SchemaError("bad strand {!r} on {}".format(self.strand, self.accession))
		if self.start < 1 or self.end < self.start:
			raise SchemaError("bad span {}..{} on {}".format(self.start, self.end, self.accession))


@dataclass(frozen=True)
class Window(Row):
	FILE: ClassVar[str] = "windows.tsv"
	row_id: str
	assembly: str
	contig: str
	contig_length: int
	lo: int
	hi: int
	scan_lo: int
	scan_hi: int
	cut_lo: int
	cut_hi: int

	def __post_init__(self):
		for lo, hi in ((self.lo, self.hi), (self.scan_lo, self.scan_hi), (self.cut_lo, self.cut_hi)):
			if lo < 1 or hi < lo or hi > self.contig_length:
				raise SchemaError("window {} has a bad span {}..{} on a contig of {}".format(
					self.row_id, lo, hi, self.contig_length))

	@property
	def subject(self) -> str:
		return bp_subject(self.assembly, self.contig)


@dataclass(frozen=True)
class GenomeFiles(Row):
	FILE: ClassVar[str] = "genomes.tsv"
	assembly: str
	source: str
	gff: str
	faa: str
	rna: str
	genome: str

	def path(self, slot: str) -> Optional[Path]:
		value = getattr(self, slot)
		return None if value == MISSING else Path(value)

	@property
	def usable(self) -> bool:
		return self.gff != MISSING and self.faa != MISSING


@dataclass(frozen=True)
class QueryTarget(Row):
	FILE: ClassVar[str] = "queries.tsv"
	query: str
	assembly: str
	accessions: str
	status: str

	@property
	def acceptable(self) -> set[str]:
		if self.accessions == MISSING:
			return {self.query}
		return set(self.accessions.split(",")) | {self.query}


@dataclass(frozen=True)
class Failure(Row):
	FILE: ClassVar[str] = "failures.tsv"
	subject: str
	reason: str


@dataclass(frozen=True)
class RowInfo(Row):
	FILE: ClassVar[str] = "rows.tsv"
	row_id: str
	query: str
	assembly: str
	accession: str
	contig: str
	strand: str
	species: str


@dataclass(frozen=True)
class Unmatched(Row):
	FILE: ClassVar[str] = "unmatched.tsv"
	query: str
	assembly: str
	reason: str


@dataclass(frozen=True)
class RangeReport(Row):
	FILE: ClassVar[str] = "range_report.tsv"
	row_id: str
	contig_length: int
	query_start: int
	query_end: int
	up_available: int
	down_available: int
	up_reached: int
	down_reached: int
	genes_up: int
	genes_down: int


@dataclass(frozen=True)
class Family(Row):
	FILE: ClassVar[str] = "families.tsv"
	family: int
	label: str
	size: int
	occurrences: int
	members: str

	@property
	def accessions(self) -> list[str]:
		return self.members.split(",")


@dataclass(frozen=True)
class ClusterHit(Row):
	FILE: ClassVar[str] = "hits.tsv"
	accession: str
	family: int
	hits: str


@dataclass(frozen=True)
class Annotation(Row):
	FILE: ClassVar[str] = "annotations.tsv"
	subject: str
	space: str
	start: Optional[int]
	end: Optional[int]
	kind: str
	category: str
	label: str
	tool: str
	score: Optional[float]

	def __post_init__(self):
		if self.kind not in KINDS:
			raise SchemaError("unknown glyph kind {!r}".format(self.kind))
		if self.space not in SPACES:
			raise SchemaError("unknown coordinate space {!r}".format(self.space))
		if self.space == MISSING:
			if self.start is not None or self.end is not None:
				raise SchemaError("coordinates given without a space on {}".format(self.subject))
		else:
			if self.start is None or self.end is None:
				raise SchemaError("space {} needs coordinates on {}".format(self.space, self.subject))
			if self.start < 1 or self.end < self.start:
				raise SchemaError("bad span {}..{} on {}".format(self.start, self.end, self.subject))

	@property
	def whole(self) -> bool:
		return self.space == MISSING
