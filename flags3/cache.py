import time
from pathlib import Path

from flags3.run import RunError
from flags3.stages.fetch import GenomeCache

DAY = 86400.0


class CacheEntry:
	def __init__(self, base: str, files: dict[str, Path]):
		self.base = base
		self.files = files
		stats = [p.stat() for p in files.values() if p.is_file()]
		self.size = sum(s.st_size for s in stats)
		self.last_used = max((max(s.st_atime, s.st_mtime) for s in stats), default=0.0)

	def age_days(self, now: float) -> float:
		return (now - self.last_used) / DAY

	def remove(self) -> None:
		for path in self.files.values():
			path.unlink(missing_ok=True)


class CacheReport:
	def __init__(self, directory: Path):
		self.directory = Path(directory)
		if not self.directory.is_dir():
			raise FileNotFoundError("{} is not a directory".format(self.directory))
		cache = GenomeCache(self.directory)
		self.entries = [CacheEntry(base, files) for base, files in sorted(cache.index.items())]
		self.now = time.time()

	@property
	def size(self) -> int:
		return sum(e.size for e in self.entries)

	def older_than(self, days: int) -> list[CacheEntry]:
		return [e for e in self.entries if e.age_days(self.now) > days]

	def summary(self) -> str:
		if not self.entries:
			return "{}: empty".format(self.directory)
		ages = [e.age_days(self.now) for e in self.entries]
		return "{}: {}, {}, last used between {:.0f} and {:.0f} days ago".format(
			self.directory, _count(len(self.entries)), _human(self.size), min(ages), max(ages))

	def prune(self, days: int, dry_run: bool = False) -> list[CacheEntry]:
		stale = self.older_than(days)
		if not dry_run:
			for entry in stale:
				entry.remove()
		return stale


def _count(n: int) -> str:
	return "{} genome{}".format(n, "" if n == 1 else "s")


def _human(size: int) -> str:
	value = float(size)
	for unit in ("B", "KB", "MB", "GB", "TB"):
		if value < 1024 or unit == "TB":
			return "{:.1f} {}".format(value, unit) if unit != "B" else "{:.0f} B".format(value)
		value /= 1024
	return "{:.1f} TB".format(value)


def main(args) -> int:
	try:
		report = CacheReport(Path(args.directory))
	except FileNotFoundError as error:
		raise RunError(str(error))
	print(report.summary())
	if args.prune is None:
		return 0
	stale = report.prune(args.prune, args.dry_run)
	verb = "would delete" if args.dry_run else "deleted"
	for entry in stale:
		print("{} {} ({}, last used {:.0f} days ago)".format(verb, entry.base, _human(entry.size), entry.age_days(report.now)))
	print("{} {}, {}, unused for more than {} days".format(verb, _count(len(stale)), _human(sum(e.size for e in stale)), args.prune))
	return 0
