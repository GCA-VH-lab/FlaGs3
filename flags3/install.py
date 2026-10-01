import gzip
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional

from flags3 import home
from flags3.tools import Tools

MICROMAMBA = "https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-{}"
CHANNELS = ("conda-forge", "bioconda")
PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
DEFENSE_MODELS = "https://api.github.com/repos/mdmparis/defense-finder-models/releases/latest"


class InstallError(RuntimeError):
	pass


def say(message: str) -> None:
	print("  " + message, flush=True)


def download(url: str, target: Path) -> Path:
	target.parent.mkdir(parents=True, exist_ok=True)
	partial = target.with_name(target.name + ".part")
	say("downloading {}".format(url))
	if shutil.which("curl"):
		done = subprocess.run(["curl", "-fL", "--retry", "3", "-C", "-", "--progress-bar", "-o", str(partial), url])
	elif shutil.which("wget"):
		done = subprocess.run(["wget", "-c", "-q", "--show-progress", "-O", str(partial), url])
	else:
		done = None
	if done is not None:
		if done.returncode != 0 or not partial.is_file() or partial.stat().st_size == 0:
			raise InstallError("download failed ({} exited {}): {}".format("curl" if shutil.which("curl") else "wget", done.returncode, url))
		os.replace(partial, target)
		return target
	request = urllib.request.Request(url, headers={"User-Agent": "flags3-install"})
	with urllib.request.urlopen(request, timeout=120) as response, open(partial, "wb") as out:
		total = int(response.headers.get("Content-Length") or 0)
		done, shown = 0, -1
		while True:
			chunk = response.read(1 << 20)
			if not chunk:
				break
			out.write(chunk)
			done += len(chunk)
			mark = done * 100 // total if total else done >> 24
			if mark != shown:
				shown = mark
				progress = "{:3d}%".format(done * 100 // total) if total else "{:,d} MB".format(done >> 20)
				sys.stdout.write("\r    {} of {}".format(progress, _size(total) if total else "?"))
				sys.stdout.flush()
		sys.stdout.write("\r    {} downloaded{}\n".format(_size(done), " " * 20))
	os.replace(partial, target)
	return target


def _size(n: int) -> str:
	return "{:.1f} GB".format(n / (1 << 30)) if n >= (1 << 30) else "{:.0f} MB".format(n / (1 << 20))


def fetch_json(url: str):
	request = urllib.request.Request(url, headers={"User-Agent": "flags3-install", "Accept": "application/vnd.github+json"})
	with urllib.request.urlopen(request, timeout=60) as response:
		return json.load(response)


def unpack(archive: Path, destination: Path) -> Path:
	destination.mkdir(parents=True, exist_ok=True)
	if archive.is_dir():
		shutil.copytree(archive, destination, dirs_exist_ok=True)
		return destination
	if tarfile.is_tarfile(archive):
		with tarfile.open(archive) as tar:
			tar.extractall(destination, filter="data")
	else:
		shutil.unpack_archive(str(archive), str(destination))
	return destination


def gunzip(source: Path, target: Path) -> Path:
	with gzip.open(source, "rb") as src, open(target, "wb") as dst:
		shutil.copyfileobj(src, dst, 1 << 20)
	return target


def run(argv: list, prefix: Optional[Path] = None, **kwargs) -> subprocess.CompletedProcess:
	say("running " + " ".join(str(a) for a in argv))
	if prefix is not None and "env" not in kwargs:
		env = dict(os.environ)
		env.pop("DYLD_LIBRARY_PATH", None)
		env["PATH"] = os.pathsep.join([str(prefix / "bin"), env.get("PATH", "")])
		kwargs["env"] = env
	done = subprocess.run([str(a) for a in argv], **kwargs)
	if done.returncode != 0:
		raise InstallError("{} exited {}".format(argv[0], done.returncode))
	return done


class Micromamba:
	def __init__(self):
		self.binary = home.TOOLS / "micromamba"
		self.root = home.HOME / "mamba"
		self.existing = shutil.which("micromamba")

	@staticmethod
	def asset() -> str:
		system, machine = platform.system(), platform.machine().lower()
		if system == "Linux":
			return "linux-aarch64" if machine in ("aarch64", "arm64") else "linux-64"
		if system == "Darwin":
			return "osx-arm64" if machine in ("arm64", "aarch64") else "osx-64"
		raise InstallError("no micromamba build for {}".format(system))

	def ensure(self) -> Path:
		if self.existing:
			self.binary = Path(self.existing)
			return self.binary
		if self.binary.is_file() and os.access(self.binary, os.X_OK):
			return self.binary
		download(MICROMAMBA.format(self.asset()), self.binary)
		self.binary.chmod(self.binary.stat().st_mode | stat.S_IEXEC)
		return self.binary

	def env(self) -> dict:
		env = dict(os.environ)
		env["MAMBA_ROOT_PREFIX"] = str(self.root)
		env.pop("DYLD_LIBRARY_PATH", None)
		return env

	def create(self, prefix: Path, packages: list[str]) -> Path:
		self.ensure()
		if (prefix / "conda-meta").is_dir():
			say("environment {} exists; updating".format(prefix))
			verb = "install"
		else:
			verb = "create"
		argv = [self.binary, verb, "-y", "-p", prefix]
		for channel in CHANNELS:
			argv += ["-c", channel]
		try:
			run(argv + packages, env=self.env())
		except InstallError:
			if platform.system() != "Darwin" or platform.machine().lower() not in ("arm64", "aarch64"):
				raise
			say("native osx-arm64 solve failed; retrying as osx-64 (runs under Rosetta 2)")
			run(argv + ["--platform", "osx-64"] + packages, env=self.env())
		return prefix

	def bin(self, prefix: Path, name: str) -> Path:
		path = prefix / "bin" / name
		if not path.is_file():
			raise InstallError("{} is missing from {}".format(name, prefix))
		return path


class Component:
	name = ""
	summary = ""
	takes_path = False
	licensed = False
	rows: dict[str, str] = {}

	def __init__(self, mamba: Micromamba, tools: Tools):
		self.mamba = mamba
		self.tools = tools

	@property
	def prefix(self) -> Path:
		return home.ENVS / self.name

	def installed(self) -> Optional[str]:
		raise NotImplementedError

	def install(self, source: Optional[Path] = None) -> None:
		raise NotImplementedError

	def update(self) -> None:
		self.install()

	def set_row(self, name: str, command: str, directory: str = "", **options) -> None:
		tool = self.tools[name]
		tool.command = command
		tool.directory = directory
		tool.options.update({k: str(v) for k, v in options.items()})
		say("{} row: {}".format(name, command))


class Core(Component):
	name = "core"
	summary = "mafft, trimal, VeryFastTree, IQ-TREE and BLAST+ for --tree and --blast_mode local"
	packages = ["mafft", "trimal", "veryfasttree", "iqtree", "blast"]
	binaries = {"mafft": "mafft", "trimal": "trimal", "veryfasttree": "VeryFastTree", "iqtree": "iqtree", "blastp": "blastp"}

	def installed(self) -> Optional[str]:
		found = [self.locate(row, binary) for row, binary in self.binaries.items()]
		if not all(found):
			return None
		return "installed" if any(str(self.prefix) in str(p) for p in found) else "on PATH"

	def locate(self, row: str, binary: str) -> Optional[Path]:
		here = self.prefix / "bin" / binary
		if here.is_file():
			return here
		found, where = self.tools[row].locate()
		if found:
			return Path(where)
		on_path = shutil.which(binary)
		return Path(on_path) if on_path else None

	def install(self, source=None) -> None:
		missing = [pkg for pkg, (row, binary) in zip(self.packages, self.binaries.items())
			if self.locate(row, binary) is None]
		if missing:
			self.mamba.create(self.prefix, missing)
		else:
			say("all core tools already present; nothing to build")
		for row, binary in self.binaries.items():
			path = self.locate(row, binary)
			if path is None:
				raise InstallError("{} is missing after the install".format(binary))
			self.set_row(row, "{} {}".format(path, self.tools.template(row)))
			say("{} from {}".format(binary, "the flags3 environment" if str(self.prefix) in str(path) else path.parent))


class Pfam(Component):
	name = "pfam"
	summary = "Pfam-A HMMs and clans for --domains (about 1.5 GB unpacked)"
	takes_path = True

	def installed(self) -> Optional[str]:
		return "installed" if home.PFAM_HMM.with_suffix(".hmm.h3m").is_file() else None

	def install(self, source: Optional[Path] = None) -> None:
		target = home.PFAM_HMM.parent
		target.mkdir(parents=True, exist_ok=True)
		if source is None:
			download(PFAM_BASE + "/Pfam-A.hmm.gz", target / "Pfam-A.hmm.gz")
			download(PFAM_BASE + "/Pfam-A.clans.tsv.gz", home.PFAM_CLANS)
			source = target / "Pfam-A.hmm.gz"
		if str(source).endswith(".gz"):
			gunzip(source, home.PFAM_HMM)
		elif source != home.PFAM_HMM:
			shutil.copy(source, home.PFAM_HMM)
		self.press(home.PFAM_HMM)

	@staticmethod
	def press(path: Path) -> int:
		import pyhmmer
		for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
			path.with_name(path.name + suffix).unlink(missing_ok=True)
		with pyhmmer.plan7.HMMFile(str(path)) as handle:
			count = pyhmmer.hmmpress(handle, str(path))
		say("pressed {} HMMs".format(count))
		return count


class InterPro(Component):
	name = "interpro"
	summary = "InterPro metadata table joined onto the domain table; give the file, it is not downloadable"
	takes_path = True

	def installed(self) -> Optional[str]:
		for path in (home.INTERPRO, home.INTERPRO.with_suffix(".tsv.gz")):
			if path.is_file():
				return "installed"
		return None

	def install(self, source: Optional[Path] = None) -> None:
		if source is None:
			raise InstallError("interpro needs the table: flags3 install interpro /path/to/interpro_metadata_processed.tsv")
		target = home.INTERPRO.with_suffix(".tsv.gz") if str(source).endswith(".gz") else home.INTERPRO
		target.parent.mkdir(parents=True, exist_ok=True)
		shutil.copy(source, target)
		say("copied to {}".format(target))


class DefenceHmm(Component):
	name = "defence-hmm"
	summary = "DefenseFinder HMM profiles for --domains -db (no DefenseFinder install needed)"
	takes_path = True

	@property
	def directory(self) -> Path:
		return home.DB / "defensefinder"

	def installed(self) -> Optional[str]:
		profiles = self.directory / "profiles"
		return "installed" if profiles.is_dir() and any(profiles.glob("*.hmm")) else None

	def install(self, source: Optional[Path] = None) -> None:
		if source is None:
			release = fetch_json(DEFENSE_MODELS)
			asset = next((a for a in release.get("assets", []) if a["name"].endswith(".tar.gz")), None)
			if asset is None:
				raise InstallError("no .tar.gz asset in the latest defense-finder-models release")
			source = download(asset["browser_download_url"], home.TOOLS / "downloads" / asset["name"])
		unpacked = unpack(source, self.directory / "unpacked")
		profiles = next((p for p in sorted(unpacked.rglob("profiles")) if p.is_dir() and any(p.glob("*.hmm"))), None)
		if profiles is None:
			raise InstallError("no profiles/ directory with .hmm files inside {}".format(source))
		target = self.directory / "profiles"
		if target.exists():
			shutil.rmtree(target)
		shutil.move(str(profiles), str(target))
		shutil.rmtree(unpacked, ignore_errors=True)
		say("use with: -d -db defensefinder={} -hc defensefinder=0.7,0.5".format(target))


class CondaTool(Component):
	packages: list[str] = []
	binary = ""
	row = ""
	command = ""

	def installed(self) -> Optional[str]:
		return "installed" if (self.prefix / "bin" / self.binary).is_file() else None

	def install(self, source=None) -> None:
		self.mamba.create(self.prefix, self.packages)
		path = self.mamba.bin(self.prefix, self.binary)
		self.set_row(self.row, "{} {}".format(path, self.command))
		self.after(path)

	def after(self, binary: Path) -> None:
		pass


class Mmseqs(Component):
	name = "mmseqs"
	summary = "MMseqs2 (static release binary) for -cm mmseqs_cluster"
	RELEASE = "https://github.com/soedinglab/MMseqs2/releases/latest/download/{}"

	@property
	def directory(self) -> Path:
		return home.TOOLS / "mmseqs"

	@property
	def binary(self) -> Path:
		return self.directory / "bin" / "mmseqs"

	@staticmethod
	def asset() -> str:
		system, machine = platform.system(), platform.machine().lower()
		if system == "Darwin":
			return "mmseqs-osx-universal.tar.gz"
		if system == "Linux":
			return "mmseqs-linux-arm64.tar.gz" if machine in ("aarch64", "arm64") else "mmseqs-linux-avx2.tar.gz"
		raise InstallError("no MMseqs2 binary for {}".format(system))

	def installed(self) -> Optional[str]:
		if self.binary.is_file():
			return "installed"
		return "on PATH" if shutil.which("mmseqs") else None

	def install(self, source=None) -> None:
		path = Path(shutil.which("mmseqs")) if (shutil.which("mmseqs") and source is None and not self.binary.is_file()) else None
		if path is None:
			archive = source or download(self.RELEASE.format(self.asset()), home.TOOLS / "downloads" / self.asset())
			unpacked = unpack(archive, home.TOOLS / "mmseqs-unpacked")
			found = next((p for p in unpacked.rglob("bin/mmseqs")), None)
			if found is None:
				raise InstallError("no bin/mmseqs inside {}".format(archive))
			if self.directory.exists():
				shutil.rmtree(self.directory)
			shutil.move(str(found.parent.parent), str(self.directory))
			shutil.rmtree(unpacked, ignore_errors=True)
			self.binary.chmod(self.binary.stat().st_mode | stat.S_IEXEC)
			path = self.binary
		else:
			say("mmseqs found on PATH at {}".format(path))
		for row in ("mmseqs", "mmseqs_cluster", "mmseqs_cluster_exhaustive"):
			self.set_row(row, "{} {}".format(path, self.tools.template(row)))


class Genomad(CondaTool):
	name = "genomad"
	summary = "geNomad and its database (1.6 GB) for --genomad"
	packages = ["genomad"]
	binary = "genomad"
	row = "genomad"
	command = "end-to-end --cleanup --threads {threads} {in} {out} {db}"

	def after(self, binary: Path) -> None:
		db = home.DB / "genomad"
		db.mkdir(parents=True, exist_ok=True)
		run([binary, "download-database", db], prefix=self.prefix)
		inner = next((p for p in db.iterdir() if p.is_dir() and p.name.startswith("genomad_db")), db)
		self.set_row(self.row, "{} {}".format(binary, self.command), db=str(inner))


class DefenseFinder(CondaTool):
	name = "defensefinder"
	summary = "DefenseFinder for --defensefinder"
	packages = ["python=3.10", "hmmer", "pip"]
	binary = "defense-finder"
	row = "defensefinder"
	command = "run --db-type gembase -o {out} {faa}"

	def install(self, source=None) -> None:
		self.mamba.create(self.prefix, self.packages)
		run([self.prefix / "bin" / "pip", "install", "--upgrade", "mdmparis-defense-finder"], prefix=self.prefix)
		path = self.mamba.bin(self.prefix, self.binary)
		self.set_row(self.row, "{} {}".format(path, self.command))
		self.after(path)

	def after(self, binary: Path) -> None:
		run([binary, "update"], prefix=self.prefix)


class Padloc(CondaTool):
	name = "padloc"
	summary = "PadLoc for --padloc"
	packages = ["padloc"]
	binary = "padloc"
	row = "padloc"
	command = "--faa {faa} --gff {gff} --outdir {out} --cpu {threads}"

	def after(self, binary: Path) -> None:
		run([binary, "--db-update"], prefix=self.prefix)


class Sismis(CondaTool):
	name = "sismis"
	summary = "Sismis for --sismis"
	packages = ["python=3.11", "hmmer", "pip"]
	binary = "sismis"
	row = "sismis"
	command = "run -g {in} -o {out}"

	def install(self, source=None) -> None:
		self.mamba.create(self.prefix, self.packages)
		run([self.prefix / "bin" / "pip", "install", "sismis"], prefix=self.prefix)
		path = self.mamba.bin(self.prefix, self.binary)
		self.set_row(self.row, "{} {}".format(path, self.command))


class SignalP(Component):
	name = "signalp"
	summary = "SignalP 6 (licensed; give the package you downloaded from DTU)"
	takes_path = True
	licensed = True

	def installed(self) -> Optional[str]:
		return "installed" if (self.prefix / "bin" / "signalp6").is_file() else None

	def install(self, source: Optional[Path] = None) -> None:
		if source is None:
			raise InstallError("signalp needs the package: flags3 install signalp /path/to/signalp-6-package[.tar.gz]")
		unpacked = unpack(source, home.TOOLS / "signalp6-src")
		setup = next((p.parent for p in unpacked.rglob("setup.py")), None) or next((p.parent for p in unpacked.rglob("pyproject.toml")), None)
		if setup is None:
			raise InstallError("no setup.py or pyproject.toml inside {}".format(source))
		self.mamba.create(self.prefix, ["python=3.9", "pip"])
		pip = self.prefix / "bin" / "pip"
		run([pip, "install", "torch<2.0", "numpy<2"], prefix=self.prefix)
		run([pip, "install", str(setup)], prefix=self.prefix)
		self.copy_weights(unpacked)
		binary = self.mamba.bin(self.prefix, "signalp6")
		self.set_row("signalp", "{} --fastafile {{fasta}} --output_dir {{out}} --organism other --format txt --mode fast".format(binary))


	def copy_weights(self, unpacked: Path) -> None:
		source = next((p for p in unpacked.rglob("models") if p.is_dir() and any(p.glob("*.pt"))), None)
		if source is None:
			source = next((p for p in unpacked.rglob("model_weights") if p.is_dir() and any(p.glob("*.pt"))), None)
		if source is None:
			say("no models/ directory with .pt files in the package; signalp6 will report missing weights")
			return
		done = subprocess.run([self.prefix / "bin" / "python", "-c", "import signalp, os; print(os.path.dirname(signalp.__file__))"],
			capture_output=True, text=True)
		if done.returncode != 0:
			raise InstallError("cannot locate the installed signalp package: {}".format(done.stderr.strip()[-200:]))
		target = Path(done.stdout.strip()) / "model_weights"
		target.mkdir(parents=True, exist_ok=True)
		copied = 0
		for path in source.iterdir():
			if path.is_file():
				shutil.copy(path, target / path.name)
				copied += 1
		say("copied {} model files to {}".format(copied, target))


class DeepTmhmm2(Component):
	name = "deeptmhmm"
	summary = "DeepTMHMM2, local and licence-free (MIT), for --tmhmm -lth; ~2 GB of torch and model weights"
	REPO = "git+https://github.com/fteufel/DeepTMHMM2.git"
	SAMPLE = ">warmup\nMKKLLIAGAVLLGLASSAWAQQTASLPADKSQQ\n"

	def installed(self) -> Optional[str]:
		return "installed" if (self.prefix / "bin" / "dtm2").is_file() else None

	def install(self, source=None) -> None:
		self.mamba.create(self.prefix, ["python=3.12", "pip"])
		pip = self.prefix / "bin" / "pip"
		run([pip, "install", "--upgrade", str(source) if source else self.REPO], prefix=self.prefix)
		binary = self.mamba.bin(self.prefix, "dtm2")
		self.set_row("deeptmhmm", "{} {}".format(binary, self.tools.template("deeptmhmm")))
		say("fetching the model weights with a warm-up prediction")
		work = home.TOOLS / "deeptmhmm2-warmup"
		work.mkdir(parents=True, exist_ok=True)
		(work / "sample.fasta").write_text(self.SAMPLE)
		run([binary, work / "sample.fasta", work / "out", "--simplify-io"], prefix=self.prefix)
		shutil.rmtree(work, ignore_errors=True)


COMPONENTS = (Core, Pfam, InterPro, DefenceHmm, Mmseqs, Genomad, DefenseFinder, Padloc, Sismis, DeepTmhmm2, SignalP)
BY_NAME = {c.name: c for c in COMPONENTS}


class Installer:
	def __init__(self):
		home.HOME.mkdir(parents=True, exist_ok=True)
		self.mamba = Micromamba()
		self.tools = Tools.load(home.user_tools_table())
		self.components = {c.name: c(self.mamba, self.tools) for c in COMPONENTS}

	def status(self) -> list[tuple[str, str, str]]:
		return [(c.name, c.installed() or "-", c.summary) for c in self.components.values()]

	def run(self, names: list[str], sources: dict[str, Optional[Path]], update: bool = False, force: bool = False) -> None:
		for name in names:
			component = self.components.get(name)
			if component is None:
				raise InstallError("unknown component {!r}; known: {}".format(name, ", ".join(BY_NAME)))
			print("== {}".format(name), flush=True)
			state = component.installed()
			if state == "installed" and not force and name not in sources:
				if update:
					component.update()
				else:
					say("already {}; --update refreshes it, --force rebuilds it".format(state))
					continue
			else:
				component.install(sources.get(name))
			self.tools.write(home.USER_TOOLS_TABLE)
			say("done; tools table at {}".format(home.USER_TOOLS_TABLE))


def main(args) -> int:
	installer = Installer()
	names = list(args.components or [])
	if args.all:
		names = [c.name for c in COMPONENTS if not c.licensed and c.name != "interpro"]
	if not names:
		width = max(len(n) for n in BY_NAME)
		for name, state, summary in installer.status():
			print("{:{w}}  {:9}  {}".format(name, state, summary, w=width))
		return 0
	sources = {}
	if args.path:
		with_path = [n for n in names if BY_NAME[n].takes_path]
		if len(with_path) != 1:
			raise InstallError("a path applies to exactly one component that takes one")
		sources[with_path[0]] = Path(args.path).expanduser().resolve()
	installer.run(names, sources, update=args.update, force=args.force)
	return 0
