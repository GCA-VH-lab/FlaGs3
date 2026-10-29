import csv
from pathlib import Path

from flags3 import home
from flags3.scan import Hit, ScanError, WindowScan
from flags3.tools import Tool

COLUMNS = ("qseqid", "sseqid", "pident", "length", "slen", "sstart", "send", "qstart", "qend", "evalue", "bitscore")


def database(tool: Tool, configured) -> Path:
	raw = configured or tool.options.get("db") or home.PLASMIDFINDER
	path = Path(raw).expanduser()
	for candidate in (path, path / "plasmidfinder_db"):
		if (candidate / "replicons.fsa").is_file():
			return candidate
	raise ScanError("no PlasmidFinder database at {}; flags3 install plasmidfinder sets it up, or give -pfdb".format(path))


def replicon_type(subject: str) -> tuple:
	header = subject.strip()
	return header.split("_", 1)[0], header


def parse_blast(path: Path, min_identity: float, min_coverage: float):
	with open(path, newline="", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t", fieldnames=COLUMNS):
			identity = float(row["pident"]) / 100.0
			length, slen = int(row["length"]), int(row["slen"])
			coverage = length / slen if slen else 0.0
			if identity < min_identity or coverage < min_coverage:
				continue
			name, header = replicon_type(row["sseqid"])
			start, end = sorted((int(row["qstart"]), int(row["qend"])))
			yield Hit(row["qseqid"], start, end, name, identity * 100.0, name, extra={
				"replicon": name, "sequence": header, "identity": "{:.2f}".format(identity * 100.0),
				"coverage": "{:.2f}".format(coverage * 100.0), "aligned": str(length), "evalue": row["evalue"]})


class PlasmidFinder(WindowScan):
	name = "plasmidfinder"
	tool_name = "plasmidfinder"
	report_file = "replicons.tsv"
	passthrough_prefix = "pf_"
	kind = "verdict"

	def wanted(self, config) -> bool:
		return config.flag("plasmidfinder")

	def scan_batch(self, tool: Tool, fasta: Path, out_dir: Path, config) -> list[Hit]:
		db = database(tool, config.text("plasmidfinder_db"))
		out_dir.mkdir(parents=True, exist_ok=True)
		table = out_dir / "blast.tsv"
		self.execute(tool, tool.argv(**{"in": fasta, "db": db / "replicons", "out": table,
			"threads": config.workers(), "evalue": tool.options.get("evalue", "1e-10")}))
		identity = float(tool.options.get("identity", "0.90"))
		coverage = float(tool.options.get("coverage", "0.60"))
		return list(parse_blast(table, identity, coverage))
