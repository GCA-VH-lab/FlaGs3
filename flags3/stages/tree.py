from pathlib import Path

from flags3 import fasta, tree
from flags3.log import note
from flags3.schema import Row
from flags3.stage import Stage
from flags3.tools import Tools
from dataclasses import dataclass
from typing import ClassVar


@dataclass(frozen=True)
class Leaf(Row):
	FILE: ClassVar[str] = "leaves.tsv"
	position: int
	row_id: str


class Tree(Stage):
	name = "tree"
	requires = ("extract",)
	optional = True

	def wanted(self, config) -> bool:
		return config.flag("tree") or config.flag("iqtree") or config.flag("tree_order")

	def run(self, run, config, out: Path) -> None:
		tools = Tools.load(config.path("tools"))
		sequences = {name: seq for name, _, seq in fasta.read(run.stage_file("extract", "queries.faa"))}
		engine = "iqtree" if config.flag("iqtree") else "veryfasttree"
		note("building a tree of {} queries with {}".format(len(sequences), engine))
		builder = tree.TreeBuilder(tools, out / "raw", threads=config.workers(), engine=engine,
			trimal_mode=config.text("trimal_mode", "gt"), trimal_value=config.number("trimal_value", 0.1),
			trimal_extra=config.text("trimal_extra"))
		newick = builder.build(sequences)
		if not newick:
			raise tree.TreeError("{} produced no tree".format(engine))
		(out / "tree.nwk").write_text(newick + "\n")
		fasta.write(out / "alignment.aln", builder.raw_alignment.items())
		fasta.write(out / "trimmed.aln", builder.alignment.items())
		(out / "commands.txt").write_text("\n".join(builder.commands) + "\n")
		Leaf.write(out / Leaf.FILE, (Leaf(i + 1, name) for i, name in enumerate(tree.leaf_order(newick))))
		note("tree with {} leaves".format(len(sequences)))
