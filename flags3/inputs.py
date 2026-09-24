from pathlib import Path

BLAST_MARK = "blast"


class InputList:
	def __init__(self):
		self.entries: list[tuple[str, str | None]] = []
		self.paired: list[tuple[str, str]] = []
		self.unpaired: list[str] = []
		self.blast: list[tuple[str, str, int]] = []
		self.warnings: list[str] = []

	def read(self, path: Path) -> "InputList":
		seen = set(self.unpaired)
		with open(path, encoding="utf-8") as handle:
			for n, line in enumerate(handle, 1):
				line = line.split("#", 1)[0].strip()
				if not line:
					continue
				if "\t" in line:
					fields = [c for c in line.split("\t") if c]
					if len(fields) < 2:
						self.warnings.append("line {} of {} is malformed, skipping it: {!r}".format(n, path, line))
					elif fields[1].strip().lower() == BLAST_MARK:
						self.blast.append((fields[0], str(path), n))
					else:
						self.paired.append((fields[0], fields[1]))
						self.entries.append((fields[0], fields[1]))
				elif line not in seen:
					seen.add(line)
					self.unpaired.append(line)
					self.entries.append((line, None))
		return self

	@property
	def queries(self) -> list[str]:
		out, seen = [], set()
		for q in [p for p, _ in self.paired] + self.unpaired:
			if q not in seen:
				seen.add(q)
				out.append(q)
		return out
