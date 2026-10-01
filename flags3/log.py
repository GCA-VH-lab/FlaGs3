import os
import sys
import threading
import time
import traceback

STDERR_PREFIX = "[stderr] "


class Transcript:
	def __init__(self):
		self.handle = None
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
		if self.handle is None:
			self.pending.append(out)
		else:
			self.handle.write(out)
			self.handle.flush()

	def attach(self, path):
		with self.lock:
			if self.handle is not None:
				return
			os.makedirs(os.path.dirname(str(path)) or ".", exist_ok=True)
			self.handle = open(path, "w", encoding="utf-8", errors="replace")
			self.path = path
			for line in self.pending:
				self.handle.write(line)
			self.pending = []
			self.handle.flush()

	def close(self):
		with self.lock:
			for tag in sorted(self.partial):
				if self.partial[tag]:
					self._line(tag, self.partial[tag])
			self.partial = {}
			if self.handle is not None:
				self.handle.close()
				self.handle = None


class Tee:
	def __init__(self, stream, tag, transcript):
		self._stream = stream
		self._tag = tag
		self._transcript = transcript

	def write(self, text):
		self._stream.write(text)
		self._transcript.write(self._tag, text)
		return len(text)

	def flush(self):
		self._stream.flush()

	def isatty(self):
		return self._stream.isatty()

	def fileno(self):
		return self._stream.fileno()

	def __getattr__(self, name):
		return getattr(self._stream, name)


class Console:
	def __init__(self):
		self.transcript = Transcript()
		self.debug_on = False
		self.verbose = False
		self.started = time.monotonic()
		self.capturing = False

	def capture(self):
		if self.capturing:
			return
		self.capturing = True
		sys.stdout = Tee(sys.stdout, "stdout", self.transcript)
		sys.stderr = Tee(sys.stderr, "stderr", self.transcript)

	def attach(self, path):
		self.capture()
		self.transcript.attach(path)

	def close(self):
		self.transcript.close()

	def debug(self, message, exc=False):
		if not self.debug_on:
			return
		sys.stderr.write("[debug {:7.2f}s] {}\n".format(
			time.monotonic() - self.started, message))
		if exc:
			traceback.print_exc(file=sys.stderr)
		sys.stderr.flush()

	def note(self, message):
		if self.verbose:
			print(">> {}".format(message), flush=True)

	def record(self, text, stream="stdout"):
		if not text:
			return
		self.transcript.write(stream, text if text.endswith("\n") else text + "\n")

	def record_command(self, cmd, returncode, stdout=None, stderr=None):
		name = cmd[0] if cmd else "?"
		self.record("--- {} (exit {}) ---".format(os.path.basename(str(name)), returncode))
		self.record(" ".join(str(c) for c in cmd))
		self.record(stdout)
		self.record(stderr, "stderr")


console = Console()
debug = console.debug
note = console.note
record = console.record
record_command = console.record_command
