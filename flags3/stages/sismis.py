import csv
from pathlib import Path

from flags3.scan import Hit, ScanError, WindowScan
from flags3.tools import Tool

CONTIG_COLUMNS = ("sequence_id", "seq_id", "contig", "scaffold", "record_id")


class Sismis(WindowScan):
	name = "sismis"
	tool_name = "sismis"
	report_file = "secretion.tsv"
	passthrough_prefix = "sismis_"

	def wanted(self, config) -> bool:
		return config.flag("sismis")

	def scan_batch(self, tool: Tool, fasta: Path, out_dir: Path, config) -> list[Hit]:
		self.execute(tool, tool.argv(**{"in": fasta, "out": out_dir}))
		table = next((p for p in out_dir.glob("*.clusters.tsv")), None) if out_dir.is_dir() else None
		if table is None:
			return []
		return list(parse_clusters(table))


def parse_clusters(table: Path):
	with open(table, newline="", encoding="utf-8") as handle:
		reader = csv.DictReader(handle, delimiter="\t")
		columns = [c.lower() for c in reader.fieldnames or []]
		contig = next((c for c in CONTIG_COLUMNS if c in columns), None)
		missing = [c for c in ("start", "end", "type", "max_p") if c not in columns]
		if contig is None or missing:
			raise ScanError("{} lacks columns {}; seen: {}".format(table, missing or ["contig"], columns))
		for row in reader:
			row = {k.lower(): v for k, v in row.items()}
			yield Hit(row[contig], int(row["start"]), int(row["end"]), row["type"], float(row["max_p"]),
				extra={k: v for k, v in row.items() if k not in (contig,)})
