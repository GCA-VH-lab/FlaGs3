import subprocess
import sys
import time
from pathlib import Path

from flags3.run import KeyValueFile

STATUS = "batch_status.tsv"
LOG = "batch.log"


class Batch:
	def __init__(self, lists: Path, out: Path, options: list[str], retry: bool):
		self.lists = sorted(p for p in lists.iterdir() if p.is_file() and not p.name.startswith("."))
		self.out = out
		self.options = options
		self.retry = retry
		self.out.mkdir(parents=True, exist_ok=True)
		self.status_path = out / STATUS
		self.log_path = out / LOG

	def done(self, name: str) -> bool:
		report = self.out / name / "report" / "status.tsv"
		return report.is_file() and KeyValueFile(report).load().get("status") == "ok"

	def record(self, name: str, state: str, seconds: float) -> None:
		new = not self.status_path.exists()
		with open(self.status_path, "a", encoding="utf-8") as out:
			if new:
				out.write("list\tstatus\tseconds\tfinished\n")
			out.write("{}\t{}\t{:.0f}\t{}\n".format(name, state, seconds, time.strftime("%Y-%m-%d %H:%M:%S")))

	def run(self) -> int:
		failures = 0
		print("{} lists in {}; runs go to {}".format(len(self.lists), self.lists[0].parent if self.lists else "-", self.out))
		for i, path in enumerate(self.lists, 1):
			name = path.stem
			if self.done(name):
				print("[{}/{}] {}: done, skipped".format(i, len(self.lists), name))
				continue
			if (self.out / name).exists() and not self.retry:
				print("[{}/{}] {}: incomplete run directory exists; pass --retry to redo it".format(i, len(self.lists), name))
				failures += 1
				continue
			argv = [sys.executable, "-m", "flags3", "run", "-i", str(path), "-o", str(self.out / name), "-nt"] + self.options
			print("[{}/{}] {}: starting".format(i, len(self.lists), name), flush=True)
			started = time.time()
			with open(self.log_path, "a", encoding="utf-8") as log:
				log.write("\n==== {} {}\n{}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), name, " ".join(argv)))
				log.flush()
				if (self.out / name).exists():
					subprocess.run(["rm", "-rf", str(self.out / name)])
				done = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
			seconds = time.time() - started
			state = "ok" if done.returncode == 0 and self.done(name) else "failed"
			failures += state == "failed"
			self.record(name, state, seconds)
			print("[{}/{}] {}: {} in {:.0f} s".format(i, len(self.lists), name, state, seconds), flush=True)
		print("batch finished: {} failed; status in {}".format(failures, self.status_path))
		return 1 if failures else 0
