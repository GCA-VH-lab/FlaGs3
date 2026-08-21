import sys
import time
import traceback

_T0 = time.monotonic()
DEBUG = False


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
