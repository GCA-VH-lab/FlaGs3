import csv
from pathlib import Path

from flags3.scan import Hit, ScanError, WindowScan
from flags3.tools import Tool

DB_MARKERS = ("genomad_marker_metadata.tsv", "version.txt")


def database(tool: Tool, configured) -> Path:
	raw = configured or tool.options.get("db")
	if not raw:
		raise ScanError("no geNomad database; flags3 install genomad downloads it, or give -gdb")
	path = Path(raw).expanduser()
	for candidate in (path, path / "genomad_db"):
		if any((candidate / m).is_file() for m in DB_MARKERS):
			return candidate
	raise ScanError("{} does not look like a geNomad database (none of {})".format(path, ", ".join(DB_MARKERS)))


def virus_name(taxonomy: str) -> str:
	ranks = [r.strip() for r in (taxonomy or "").split(";") if r.strip()]
	ranks = [r for r in ranks if r.lower() not in ("viruses", "unclassified")]
	return ranks[-1] if ranks else "virus"


def plasmid_name(conjugation: str, amr: str) -> str:
	parts = []
	if (conjugation or "").strip() not in ("", "0", "NA"):
		parts.append("conjugative")
	parts.append("plasmid")
	if (amr or "").strip() not in ("", "0", "NA"):
		parts.append("(AMR)")
	return " ".join(parts)


def parse_summary(path: Path, kind: str):
	with open(path, newline="", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t"):
			name = (row.get("seq_name") or "").strip()
			if not name:
				continue
			record = name.split("|")[0]
			coords = (row.get("coordinates") or "").strip()
			if coords and coords.upper() != "NA" and "-" in coords:
				lo, _, hi = coords.partition("-")
				start, end = int(lo), int(hi)
			else:
				start, end = 1, int(row.get("length") or 0)
			if kind == "virus":
				label, score = virus_name(row.get("taxonomy", "")), row.get("virus_score")
			else:
				label, score = plasmid_name(row.get("conjugation_genes", ""), row.get("amr_genes", "")), row.get("plasmid_score")
			yield Hit(record, start, end, kind, float(score or 0.0), label, extra=dict(row))


class Genomad(WindowScan):
	name = "genomad"
	tool_name = "genomad"
	report_file = "mobile_elements.tsv"
	passthrough_prefix = "genomad_"

	def wanted(self, config) -> bool:
		return config.flag("genomad")

	def scan_batch(self, tool: Tool, fasta: Path, out_dir: Path, config) -> list[Hit]:
		db = database(tool, config.text("genomad_db"))
		self.execute(tool, tool.argv(**{"in": fasta, "out": out_dir, "db": db, "threads": config.workers()}))
		hits = []
		for path in sorted(out_dir.rglob("*_summary.tsv")):
			if path.name.endswith("_virus_summary.tsv"):
				hits.extend(parse_summary(path, "virus"))
			elif path.name.endswith("_plasmid_summary.tsv"):
				hits.extend(parse_summary(path, "plasmid"))
		return hits
