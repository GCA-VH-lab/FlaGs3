import csv
import os
import shutil
import shlex

TABLE_NAME = "tools_table.tsv"
LOCAL_NAME = "tools_table.local.tsv"

COLUMNS = ("name", "command", "directory", "scan_range")

DEFAULTS = {
	"mafft": ("mafft --auto --anysymbol --quiet --thread {threads} {in}", ""),
	"trimal": ("trimal -in {in} -out {out} -fasta {mode}", ""),
	"veryfasttree": ("VeryFastTree {in}", ""),
	"iqtree": ("iqtree -s {in} -m {model} --prefix {prefix} -T {threads} --quiet", ""),
	"blastp": ("blastp -query {in} -db {db} -outfmt \"6 sacc evalue bitscore stitle\" "
			   "-evalue {evalue} -max_target_seqs {hits}", ""),
	"sismis": ("sismis run -g {in} -o {out}", ""),
	"defensefinder": ("defense-finder run --db-type gembase -o {out} {faa}", ""),
	"padloc": ("padloc --faa {faa} --gff {gff} --outdir {out} --cpu {threads}", ""),
	"genomad": ("genomad end-to-end --cleanup --threads {threads} {in} {out} {db}", ""),
	"mmseqs": ("mmseqs easy-linclust {in} {out} {tmp} --min-seq-id {id} "
			   "-c {cov} --cov-mode 0 --threads {threads} -v 1", ""),
	"deeptmhmm": ("python3 predict.py --fasta {fasta} --output-dir {out}", ""),
	"signalp": ("signalp6 --fastafile {fasta} --output_dir {out} --organism other "
				"--format txt --mode fast", ""),
}

_loaded = None
_scan = {}


def table_path(explicit=None) -> str:
	if explicit:
		return explicit
	return os.path.join(os.path.dirname(os.path.abspath(__file__)), TABLE_NAME)


def local_path(path=None):
	return os.path.join(os.path.dirname(table_path(path)), LOCAL_NAME)


def load(path=None):
	"""Read the shipped defaults, then tools_table.tsv, then
	tools_table.local.tsv. The local file is what the installers write, so
	machine-specific paths never end up in the committed table."""
	global _loaded, _scan
	tools = {name: (cmd, wd) for name, (cmd, wd) in DEFAULTS.items()}
	scan = {}
	targets = [table_path(path)]
	if path is None:
		targets.append(local_path())
	for target in targets:
		_read_into(target, tools, scan, required=bool(path))
	_loaded = tools
	_scan = scan
	return tools


def _read_into(target, tools, scan, required=False):
	if os.path.isfile(target):
		with open(target, newline="", encoding="utf-8") as fh:
			for row in csv.DictReader(
					(ln for ln in fh if ln.strip() and not ln.startswith("##")),
					delimiter="\t"):
				clean = {(k or "").lstrip("#").strip(): (v or "").strip()
						 for k, v in row.items()}
				name = clean.get("name", "").lower()
				if not name:
					continue
				if name not in DEFAULTS:
					raise ValueError(
						"{}: unknown tool {!r}; known tools: {}".format(
							target, name, ", ".join(sorted(DEFAULTS))))
				known = tools.get(name, DEFAULTS[name])
				command = clean.get("command") or known[0]
				directory = clean.get("directory", "")
				tools[name] = (command,
							   os.path.expanduser(directory) if directory else known[1])
				span = clean.get("scan_range", "")
				if span:
					if span.lower() == "genome":
						scan[name] = 0
					else:
						try:
							scan[name] = int(span)
						except ValueError:
							raise ValueError(
								"{}: scan_range for {} must be a number of bases "
								"or 'genome', got {!r}".format(target, name, span))
	elif required:
		raise FileNotFoundError("tool table not found: {}".format(target))


def scan_range(name: str, default):
	if _loaded is None:
		load()
	return _scan.get(name, default)


def get(name: str):
	if _loaded is None:
		load()
	return _loaded.get(name, DEFAULTS.get(name, ("", "")))


def command(name: str, **values):
	template, directory = get(name)
	out = []
	for part in shlex.split(template):
		filled = part
		for key, value in values.items():
			filled = filled.replace("{%s}" % key, str(value))
		if not filled:
			continue
		if filled != part and " " in filled and part.startswith("{"):
			out.extend(shlex.split(filled))
		else:
			out.append(filled)
	return out, directory


def missing_values(name: str, values: dict):
	cmd, _ = get(name)
	empty = [k for k in values
			 if "{" + k + "}" in cmd and not str(values.get(k, "")).strip()]
	return empty


def locate(cmd):
	if not cmd:
		return False, "no command"
	program = str(cmd[0])
	if os.path.isabs(program):
		if os.path.isfile(program) and os.access(program, os.X_OK):
			return True, program
		return False, "{} is not an executable file".format(program)
	found = shutil.which(program)
	if found:
		return True, found
	return False, "{} not found on PATH".format(program)


def brief(text, lines: int = 3, limit: int = 300) -> str:
	kept = [ln.strip() for ln in str(text or "").splitlines() if ln.strip()]
	if not kept:
		return "no output"
	joined = " | ".join(kept[-lines:])
	return joined if len(joined) <= limit else "..." + joined[-limit:]


def env_for(cmd):
	if not cmd:
		return None
	program = str(cmd[0])
	if not os.path.isabs(program):
		return None
	bindir = os.path.dirname(program)
	if not bindir or not os.path.isdir(bindir):
		return None
	env = dict(os.environ)
	parts = [p for p in env.get("PATH", "").split(os.pathsep) if p]
	if bindir in parts:
		return None
	env["PATH"] = os.pathsep.join([bindir] + parts)
	return env


def write_default_table(path: str):
	with open(path, "w") as out:
		out.write("#" + "\t".join(COLUMNS) + "\n")
		for name in sorted(DEFAULTS):
			cmd, wd = DEFAULTS[name]
			out.write("{}\t{}\t{}\t\n".format(name, cmd, wd))
