import csv
import os
import shlex

TABLE_NAME = "tools_table.tsv"

COLUMNS = ("name", "command", "directory")

DEFAULTS = {
	"mafft": ("mafft --auto --anysymbol --quiet --thread {threads} {in}", ""),
	"trimal": ("trimal -in {in} -out {out} -fasta {mode}", ""),
	"veryfasttree": ("VeryFastTree {in}", ""),
	"iqtree": ("iqtree -s {in} -m {model} --prefix {prefix} -T {threads} --quiet", ""),
	"blastp": ("blastp -query {in} -db {db} -outfmt \"6 sacc evalue bitscore stitle\" "
			   "-evalue {evalue} -max_target_seqs {hits}", ""),
	"sismis": ("sismis run -g {in} -o {out}", ""),
	"deeptmhmm": ("python3 predict.py --fasta {fasta} --output-dir {out}", ""),
	"signalp": ("signalp6 --fastafile {fasta} --output_dir {out} --organism other "
				"--format txt --mode fast", ""),
}

_loaded = None


def table_path(explicit=None) -> str:
	if explicit:
		return explicit
	return os.path.join(os.path.dirname(os.path.abspath(__file__)), TABLE_NAME)


def load(path=None):
	global _loaded
	tools = {name: (cmd, wd) for name, (cmd, wd) in DEFAULTS.items()}
	target = table_path(path)
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
				command = clean.get("command") or DEFAULTS[name][0]
				tools[name] = (command, os.path.expanduser(clean.get("directory", "")))
	elif path:
		raise FileNotFoundError("tool table not found: {}".format(path))
	_loaded = tools
	return tools


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


def write_default_table(path: str):
	with open(path, "w") as out:
		out.write("#" + "\t".join(COLUMNS) + "\n")
		for name in sorted(DEFAULTS):
			cmd, wd = DEFAULTS[name]
			out.write("{}\t{}\t{}\n".format(name, cmd, wd))
