import os
import sys
import threading
import time
import traceback

_T0 = time.monotonic()
DEBUG = False

STDERR_PREFIX = "[stderr] "


def set_debug(on: bool):
	global DEBUG
	DEBUG = bool(on)


def debug(message: str, exc: bool = False):
	if not DEBUG:
		return
	sys.stderr.write("[debug {:7.2f}s] {}\n".format(time.monotonic() - _T0, message))
	if exc:
		traceback.print_exc(file=sys.stderr)
	sys.stderr.flush()


class _Transcript:
	def __init__(self):
		self.fh = None
		self.pending = []
		self.partial = {}
		self.lock = threading.Lock()
		self.path = None

	def write(self, tag, text):
		if not text:
			return
		with self.lock:
			text = self.partial.pop(tag, "") + text
			lines = text.split("\n")
			self.partial[tag] = lines.pop()
			for line in lines:
				self._line(tag, line)

	def _line(self, tag, line):
		out = (STDERR_PREFIX + line if tag == "stderr" else line) + "\n"
		if self.fh is None:
			self.pending.append(out)
		else:
			self.fh.write(out)
			self.fh.flush()

	def attach(self, path):
		with self.lock:
			if self.fh is not None:
				return
			parent = os.path.dirname(path)
			if parent:
				os.makedirs(parent, exist_ok=True)
			self.fh = open(path, "w", encoding="utf-8", errors="replace")
			self.path = path
			for line in self.pending:
				self.fh.write(line)
			self.pending = []
			self.fh.flush()

	def close(self):
		with self.lock:
			for tag in sorted(self.partial):
				if self.partial[tag]:
					self._line(tag, self.partial[tag])
			self.partial = {}
			if self.fh is not None:
				self.fh.close()
				self.fh = None


_transcript = _Transcript()


class _Tee:
	def __init__(self, stream, tag):
		self._stream = stream
		self._tag = tag

	def write(self, text):
		self._stream.write(text)
		_transcript.write(self._tag, text)
		return len(text)

	def flush(self):
		self._stream.flush()

	def isatty(self):
		return self._stream.isatty()

	def fileno(self):
		return self._stream.fileno()

	def __getattr__(self, name):
		return getattr(self._stream, name)


_buffering = False


def start_buffering():
	global _buffering
	if _buffering:
		return
	_buffering = True
	sys.stdout = _Tee(sys.stdout, "stdout")
	sys.stderr = _Tee(sys.stderr, "stderr")


def attach(path: str):
	if not _buffering:
		start_buffering()
	_transcript.attach(path)


def record(text: str, stream: str = "stdout"):
	if not text:
		return
	_transcript.write(stream, text if text.endswith("\n") else text + "\n")


def record_command(cmd, returncode, stdout=None, stderr=None):
	name = cmd[0] if cmd else "?"
	record("--- {} (exit {}) ---".format(os.path.basename(str(name)), returncode))
	record(" ".join(str(c) for c in cmd))
	record(stdout)
	record(stderr, "stderr")


def transcript_path():
	return _transcript.path


def close():
	_transcript.close()
