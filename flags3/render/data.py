from dataclasses import dataclass, field
from pathlib import Path

from flags3.schema import Annotation, Gene, RowInfo, Window
from flags3.stages.tree import Leaf


@dataclass
class RunData:
	genes: dict[str, list[Gene]]
	rows: dict[str, RowInfo]
	windows: dict[str, Window]
	order: list[str]
	tree_order: list[str]
	newick: str
	annotations: dict[str, list[Annotation]] = field(default_factory=dict)

	@classmethod
	def load(cls, run) -> "RunData":
		genes: dict[str, list[Gene]] = {}
		for g in Gene.read(run.stage_file("extract", Gene.FILE)):
			genes.setdefault(g.row_id, []).append(g)
		for row in genes.values():
			row.sort(key=lambda g: g.start)
		rows = {r.row_id: r for r in RowInfo.read(run.stage_file("extract", RowInfo.FILE))}
		windows = {w.row_id: w for w in Window.read(run.stage_file("extract", Window.FILE))}
		order = [r for r in rows if r in genes]
		tree_order, newick = [], ""
		if run.has("tree"):
			tree_order = [l.row_id for l in Leaf.read(run.stage_file("tree", Leaf.FILE)) if l.row_id in genes]
			newick = run.stage_file("tree", "tree.nwk").read_text().strip()
		annotations = {}
		for stage in run.stages():
			path = run.stage_file(stage, Annotation.FILE)
			if run.has(stage) and path.is_file():
				found = Annotation.read(path)
				if found:
					annotations[stage] = found
		return cls(genes, rows, windows, order, tree_order, newick, annotations)

	def labels(self) -> dict[str, str]:
		per_query: dict[str, int] = {}
		for r in self.rows.values():
			per_query[r.query] = per_query.get(r.query, 0) + 1
		out = {}
		for row_id, r in self.rows.items():
			name = row_id if per_query[r.query] > 1 else r.query
			out[row_id] = "{}  {}".format(name, r.species) if r.species != "-" else name
		return out
