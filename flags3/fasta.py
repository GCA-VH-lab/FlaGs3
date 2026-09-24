import gzip
from pathlib import Path
from typing import Iterable, Iterator

COMPLEMENT = str.maketrans("ACGTUNRYKMSWBDHVacgtunrykmswbdhv",
	"TGCAANYRMKSWVHDBtgcaanyrmkswvhdb")


def open_text(path: Path):
	path = Path(path)
	if path.suffix == ".gz":
		return gzip.open(path, "rt", encoding="utf-8", errors="replace")
	return open(path, "rt", encoding="utf-8", errors="replace")


def read(path: Path) -> Iterator[tuple[str, str, str]]:
	name, description, chunks = None, "", []
	with open_text(path) as handle:
		for line in handle:
			if line.startswith(">"):
				if name is not None:
					yield name, description, "".join(chunks)
				header = line[1:].strip()
				name, _, description = header.partition(" ")
				chunks = []
			else:
				chunks.append(line.strip())
	if name is not None:
		yield name, description, "".join(chunks)


def write(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> int:
	count = 0
	with open(path, "w", encoding="utf-8") as out:
		for name, seq in records:
			out.write(">{}\n".format(name))
			if width <= 0:
				out.write(seq + "\n")
			else:
				for i in range(0, len(seq), width):
					out.write(seq[i:i + width] + "\n")
			count += 1
	return count


def reverse_complement(seq: str) -> str:
	return seq.translate(COMPLEMENT)[::-1]
