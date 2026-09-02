import os
import re
import shlex
import shutil
import subprocess
import tempfile
import threading
from typing import Dict, List, Tuple

FeatureRegion = Tuple[str, int, int]

_SUBMIT_LOCK = threading.Lock()   # pybiolib sign-in is not thread-safe; see Architecture.md
_WARMUP_LOCK = threading.Lock()
_warmed_up = False


def _debug(message, exc=False):
	try:
		from flags_log import debug
		debug(message, exc=exc)
	except ImportError:
		pass


def write_report(features: Dict[str, List[FeatureRegion]], path: str):
	with open(path, "w") as out:
		out.write("#protein\tkind\tstart\tend\n")
		for protein in sorted(features):
			for kind, start, end in sorted(features[protein], key=lambda r: r[1]):
				out.write("{}\t{}\t{}\t{}\n".format(protein, kind, start, end))


def warm_up() -> bool:

	global _warmed_up
	with _WARMUP_LOCK:
		if _warmed_up:
			return True
		try:
			from biolib.biolib_api_client import BiolibApiClient
		except ImportError:
			return False
		try:
			BiolibApiClient.get()
		except Exception:
			pass   
		_warmed_up = True
		return True


def _runs(topology: str, wanted: str, kind: str) -> List[FeatureRegion]:
	regions = []
	start = None
	for i, code in enumerate(topology):
		if code in wanted:
			if start is None:
				start = i
		elif start is not None:
			regions.append((kind, start + 1, i))
			start = None
	if start is not None:
		regions.append((kind, start + 1, len(topology)))
	return regions


def parse_deeptmhmm_3line(text: str, want_signal: bool = True) -> Dict[str, List[FeatureRegion]]:
	features: Dict[str, List[FeatureRegion]] = {}
	lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
	i = 0
	while i + 2 < len(lines) + 1 and i < len(lines):
		if not lines[i].startswith(">"):
			i += 1
			continue
		name = lines[i][1:].split("|")[0].split()[0].strip()
		topology = lines[i + 2] if i + 2 < len(lines) else ""
		regions = _runs(topology, "MB", "tm")
		if want_signal:
			regions += _runs(topology, "S", "signal")
		if regions:
			features[name] = regions
		i += 3
	return features


def parse_signalp_regions(text: str) -> Dict[str, List[FeatureRegion]]:
	features: Dict[str, List[FeatureRegion]] = {}
	lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
	i = 0
	while i < len(lines):
		if not lines[i].startswith(">"):
			i += 1
			continue
		name = lines[i][1:].split("|")[0].split()[0].strip()
		labels = lines[i + 2] if i + 2 < len(lines) else ""
		regions = _runs(labels, "STLP", "signal")
		if regions:
			features[name] = [regions[0]]
		i += 3
	return features


class _LocalScanner:

	def __init__(self, tool: str, needs=None):
		import flags_tools
		self.tool = tool
		self.command, self.directory = flags_tools.get(tool)
		self.directory = os.path.expanduser(self.directory or "")
		self.needs = needs

	def _find(self, name):
		if os.path.sep in name:
			expanded = os.path.expanduser(name)
			return expanded if os.path.isfile(expanded) else ""
		if self.directory:
			local = os.path.join(self.directory, name)
			if os.path.isfile(local):
				return local
		return shutil.which(name) or ""

	def available(self):
		cmd = shlex.split(self.command)
		if not cmd:
			return "no command configured for {}".format(self.tool)
		if self.directory and not os.path.isdir(self.directory):
			return "{} is not a directory".format(self.directory)
		if not self._find(cmd[0]):
			return "{} not found in {}".format(
				cmd[0], self.directory or "PATH" if os.path.sep not in cmd[0]
				else "the path given")
		for extra in (self.needs or []):
			if os.path.sep in extra:
				continue
			if self.directory and os.path.isfile(os.path.join(self.directory, extra)):
				continue
			if os.path.isfile(extra):
				continue
			return "{} not found in {}".format(extra, self.directory or "the working directory")
		return ""

	def _run(self, sequences: Dict[str, str], result_suffix: str) -> str:
		tmp = tempfile.mkdtemp(prefix="flags_local_")
		fasta = os.path.join(tmp, "query.fasta")
		out_dir = os.path.join(tmp, "out")
		with open(fasta, "w") as out:
			for name, seq in sequences.items():
				out.write(">{}\n{}\n".format(name, seq))
		import flags_tools
		cmd, _ = flags_tools.command(self.tool, fasta=fasta, out=out_dir)
		if self.directory:
			candidate = os.path.join(self.directory, cmd[0])
			if os.path.isfile(candidate):
				cmd[0] = candidate
		_debug("features: running {} (cwd={})".format(
			" ".join(cmd), self.directory or os.getcwd()))
		result = subprocess.run(cmd, cwd=self.directory or None,
								capture_output=True, text=True)
		_debug("features: exit {}".format(result.returncode))
		if result.returncode != 0:
			raise RuntimeError("{} exited {}: {}".format(
				cmd[0], result.returncode, (result.stderr or "").strip()[:300]))
		for root in (out_dir, tmp, self.directory or tmp):
			for base, _, names in os.walk(root):
				for name in names:
					if name.endswith(result_suffix):
						with open(os.path.join(base, name)) as fh:
							return fh.read()
		raise FileNotFoundError(
			"{} produced no {} under {}".format(cmd[0], result_suffix, tmp))


class _BioLibScanner:

	def __init__(self):
		import biolib 
		self._biolib = biolib

	def _write_fasta(self, sequences: Dict[str, str], path: str):
		with open(path, "w") as out:
			for name, seq in sequences.items():
				out.write(">{}\n{}\n".format(name, seq))

	def _run_batched(self, app_slug, sequences, result_suffix, parse,
					 args_template="--fasta {fasta}", batch_size=25):
		items = list(sequences.items())
		merged: Dict[str, List[FeatureRegion]] = {}
		for i in range(0, len(items), batch_size):
			chunk = dict(items[i:i + batch_size])
			text = self._run(app_slug, chunk, result_suffix, args_template=args_template)
			for name, regions in parse(text).items():
				merged.setdefault(name, []).extend(regions)
		return merged

	def _run(self, app_slug: str, sequences: Dict[str, str], result_suffix: str,
			 args_template: str = "--fasta {fasta}") -> str:
		tmp = tempfile.mkdtemp(prefix="flags_biolib_")
		fasta_name = "query.fasta"
		self._write_fasta(sequences, os.path.join(tmp, fasta_name))
		with _SUBMIT_LOCK:
			prev_cwd = os.getcwd()
			os.chdir(tmp)   # fasta must reach app.cli() as a bare relative name; see Architecture.md
			try:
				app = self._biolib.load(app_slug)
				job = app.cli(args=args_template.format(fasta=fasta_name))
			finally:
				os.chdir(prev_cwd)
		if hasattr(job, "wait"):
			job.wait()  

		out_dir = os.path.join(tmp, "out")
		saved = []
		problems = []
		try:
			job.save_files(out_dir)
			for root, _, names in os.walk(out_dir):
				for name in names:
					full = os.path.join(root, name)
					saved.append(full)
					if full.replace("\\", "/").endswith(result_suffix):
						with open(full) as fh:
							return fh.read()
		except Exception as e:
			problems.append("save_files: {!r}".format(e))
		try:
			listed = job.list_output_files()
			paths = [f if isinstance(f, str) else getattr(f, "path", str(f)) for f in listed]
			match = next((p for p in paths if p.replace("\\", "/").endswith(result_suffix)), None)
			if match is not None:
				out_file = job.get_output_file(match)
				data = (out_file.get_data() if hasattr(out_file, "get_data")
						else out_file.get_file_handle().read())
				return data.decode() if isinstance(data, bytes) else data
			saved = saved or paths
		except Exception as e:
			problems.append("list_output_files: {!r}".format(e))
		stdout = ""
		try:
			stdout = job.get_stdout().decode(errors="replace")[-500:]
		except Exception as e:
			problems.append("get_stdout: {!r}".format(e))
		raise FileNotFoundError(
			"{} produced no {}.\n  files seen: {}\n  retrieval errors: {}\n"
			"  stdout tail: {}".format(
				app_slug, result_suffix, saved or "(none)",
				"; ".join(problems) or "(none)", stdout or "(none)"))


class TMScanner(_BioLibScanner):
	APP = "DTU/DeepTMHMM"
	RESULT = "predicted_topologies.3line"

	def __init__(self, local=False):
		self.local = _LocalScanner("deeptmhmm", needs=["predict.py"]) if local else None
		if self.local is None:
			_BioLibScanner.__init__(self)

	def scan(self, sequences: Dict[str, str], want_signal: bool = True
			 ) -> Dict[str, List[FeatureRegion]]:
		if not sequences:
			return {}
		parse = lambda text: parse_deeptmhmm_3line(text, want_signal=want_signal)
		if self.local:
			missing = self.local.available()
			if missing:
				raise FileNotFoundError(
					"local DeepTMHMM unusable: {}. Run deeptmhmm_installer.sh to "
					"install it and fill in the deeptmhmm row of tools_table.tsv, "
					"or drop -lth to use the BioLib cloud.".format(missing))
			return parse(self.local._run(sequences, self.RESULT))
		return self._run_batched(self.APP, sequences, self.RESULT, parse)


class SignalPScanner(_BioLibScanner):

	APP = "DTU/SignalP-6"

	RESULT = "prediction_results.txt"
	ARGS = "--fastafile {fasta} --output_dir output --organism other --format txt --mode fast"

	def __init__(self, local=False):
		self.local = _LocalScanner("signalp") if local else None
		if self.local is None:
			_BioLibScanner.__init__(self)

	def scan(self, sequences: Dict[str, str]) -> Dict[str, List[FeatureRegion]]:
		if not sequences:
			return {}
		if self.local:
			missing = self.local.available()
			if missing:
				raise FileNotFoundError(
					"local SignalP unusable: {}. Run signalp_installer.sh to "
					"install it and fill in the signalp row of tools_table.tsv, "
					"or drop -lsp to use the BioLib cloud.".format(missing))
			return self._parse_prediction_results(
				self.local._run(sequences, self.RESULT))
		return self._run_batched(self.APP, sequences, self.RESULT,
								 self._parse_prediction_results, args_template=self.ARGS)

	@staticmethod
	def _parse_prediction_results(text: str) -> Dict[str, List[FeatureRegion]]:
		features: Dict[str, List[FeatureRegion]] = {}
		cs_re = re.compile(r"CS pos:\s*(\d+)")
		for line in text.splitlines():
			if line.startswith("#") or not line.strip():
				continue
			col = line.split("\t")
			if len(col) < 2:
				continue
			name = col[0].split()[0].strip()
			pred = col[1].strip().upper()
			if pred in ("", "OTHER"):
				continue   
			end = None
			m = cs_re.search(line)
			if m:
				end = int(m.group(1))
			if end is None:
				continue
			features.setdefault(name, []).append(("signal", 1, end))
		return features