from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from flags3 import fasta, features
from flags3.log import note
from flags3.schema import Annotation, Row
from flags3.stage import Stage
from flags3.tools import Tools

GLYPH = {"tm": ("segment", "TM"), "signal": ("triangle", "SP")}


@dataclass(frozen=True)
class Feature(Row):
	FILE: ClassVar[str] = "features.tsv"
	protein: str
	kind: str
	start: int
	end: int
	tool: str


class Features(Stage):
	name = "features"
	requires = ("extract",)
	optional = True
	background = True

	def wanted(self, config) -> bool:
		return config.flag("tmhmm") or config.flag("signalp")

	def run(self, run, config, out: Path) -> None:
		tools = Tools.load(config.path("tools"))
		sequences = {n: s for n, _, s in fasta.read(run.stage_file("extract", "proteins.faa"))}
		scanners: list[tuple[str, object]] = []
		problems = []
		if config.flag("tmhmm"):
			local = config.flag("local_tmhmm")
			app = config.text("tmhmm_app") or features.TmScanner.APP
			runner = (features.LocalRunner(tools["deeptmhmm"], out / "raw" / "deeptmhmm") if local
				else features.BioLibRunner(app, features.TmScanner.remote_args(app), features.TmScanner.BATCH))
			note("DeepTMHMM on {} proteins ({})".format(len(sequences), "local DeepTMHMM2" if local else "BioLib " + app))
			scanners.append(("deeptmhmm", features.TmScanner(runner, want_signal=not config.flag("signalp"))))
		if config.flag("signalp"):
			local = config.flag("local_signalp")
			runner = (features.LocalRunner(tools["signalp"], out / "raw" / "signalp") if local
				else features.BioLibRunner(features.SignalPScanner.APP, features.SignalPScanner.ARGS, features.SignalPScanner.BATCH))
			note("SignalP on {} proteins ({})".format(len(sequences), "local" if local else "BioLib"))
			scanners.append(("signalp", features.SignalPScanner(runner)))
		started = []
		for name, scanner in scanners:
			try:
				scanner.start(sequences)
				started.append((name, scanner))
			except features.FeatureError as error:
				problems.append("{}: {}".format(name, error))
		found: list[Feature] = []
		for name, scanner in started:
			try:
				found += [Feature(r.protein, r.kind, r.start, r.end, name) for r in scanner.finish()]
			except features.FeatureError as error:
				problems.append("{}: {}".format(name, error))
		for problem in problems:
			print("Warning: " + problem)
		found.sort(key=lambda f: (f.protein, f.start, f.end))
		Feature.write(out / Feature.FILE, found)
		Annotation.write(out / Annotation.FILE, (
			Annotation(f.protein, "aa", f.start, f.end, GLYPH[f.kind][0], f.kind, GLYPH[f.kind][1], f.tool, None)
			for f in found))
		note("{} regions on {} proteins".format(len(found), len({f.protein for f in found})))
		if problems and not found:
			raise RuntimeError("; ".join(problems))
