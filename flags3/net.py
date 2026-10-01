import gzip
import os
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from flags3 import VERSION
from flags3.log import debug

USER_AGENT = "flags3/" + VERSION


class RateLimiter:
	def __init__(self, rate: float):
		self.min_interval = 1.0 / rate if rate > 0 else 0.0
		self._lock = threading.Lock()
		self._next_slot = time.monotonic()

	def wait(self) -> None:
		if self.min_interval <= 0:
			return
		with self._lock:
			now = time.monotonic()
			start = max(now, self._next_slot)
			self._next_slot = start + self.min_interval
		delay = start - now
		if delay > 0:
			time.sleep(delay)


def session(workers: int) -> requests.Session:
	s = requests.Session()
	s.headers["User-Agent"] = USER_AGENT
	retry = Retry(total=5, backoff_factor=0.5, respect_retry_after_header=True,
		status_forcelist=[429, 500, 502, 503, 504], allowed_methods=frozenset(["GET"]))
	size = max(workers, 1) * 4
	s.mount("https://", HTTPAdapter(max_retries=retry, pool_connections=size, pool_maxsize=size))
	return s


class Downloader:
	def __init__(self, rate: float, workers: int):
		self.limiter = RateLimiter(rate)
		self.session = session(workers)
		self.failures: dict[str, str] = {}

	def get(self, url: str, timeout: int = 30) -> requests.Response:
		self.limiter.wait()
		return self.session.get(url, timeout=timeout)

	def stream(self, url: str, local: Path) -> bool:
		partial = local.with_name(local.name + ".part")
		return self._requests(url, local, partial)

	def _requests(self, url: str, local: Path, partial: Path) -> bool:
		for attempt in range(3):
			try:
				self.limiter.wait()
				with self.session.get(url, stream=True, timeout=120) as r:
					if r.status_code != 200:
						self.failures[url] = "HTTP {}".format(r.status_code)
						debug("HTTP {} {}".format(r.status_code, url))
						return False
					size = 0
					with open(partial, "wb") as out:
						for chunk in r.iter_content(chunk_size=1 << 16):
							out.write(chunk)
							size += len(chunk)
				if size > 0:
					os.replace(partial, local)
					self.failures.pop(url, None)
					debug("got {:,d} B {}".format(size, local.name))
					return True
				self.failures[url] = "empty response body"
			except Exception as error:
				self.failures[url] = "{}: {}".format(type(error).__name__, error)
				debug("attempt {}/3 failed for {}: {!r}".format(attempt + 1, url, error))
				time.sleep(0.5)
		partial.unlink(missing_ok=True)
		return False


def intact(path: Path) -> bool:
	try:
		if path.stat().st_size == 0:
			return False
		if path.suffix == ".gz":
			with gzip.open(path, "rb") as handle:
				while handle.read(1 << 20):
					pass
		return True
	except (OSError, EOFError, gzip.BadGzipFile):
		return False
