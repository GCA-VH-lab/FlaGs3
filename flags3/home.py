from pathlib import Path

HOME = Path.home() / ".flags3"
GENOMES = HOME / "genomes"
DB = HOME / "db"
TOOLS = HOME / "tools"
ENVS = HOME / "envs"
USER_TOOLS_TABLE = HOME / "tools_table.tsv"
PFAM_HMM = DB / "pfam" / "Pfam-A.hmm"
PFAM_CLANS = DB / "pfam" / "Pfam-A.clans.tsv.gz"
INTERPRO = DB / "interpro" / "interpro_metadata_processed.tsv"


def user_tools_table() -> Path | None:
	return USER_TOOLS_TABLE if USER_TOOLS_TABLE.is_file() else None
