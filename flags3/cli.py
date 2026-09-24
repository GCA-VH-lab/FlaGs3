import argparse
import atexit
import shutil
import sys
import time
from pathlib import Path

from flags3 import VERSION, home
from flags3.log import console
from flags3.run import RunDir, RunError
from flags3.inputs import InputList
from flags3.install import InstallError
from flags3.stage import Runner, StageError
from flags3.stages import BY_NAME, PIPELINE

RUN_OPTIONS = (
	(("-g", "--gene"), dict(type=int, default=4, help="Number of flanking genes up/downstream. Default = 4")),
	(("-r", "--range"), dict(type=int, metavar="BP", help="Neighbourhood as every gene within this many bases of the query; overrides -g")),
	(("-sr", "--scan_range"), dict(type=int, metavar="BP", help="Span around the query handed to sequence-scanning tools; without it, the whole contig")),
	(("-sm", "--scan_margin"), dict(type=int, default=10000, metavar="BP", help="Extra sequence beyond the scan range. Default = 10000")),
	(("-bi", "--blast_input"), dict(metavar="FILE", help="One starting point for BlastP: a protein accession, or a sequence (FASTA or bare). Its homologues become queries, added to -i")),
	(("-bm", "--blast_mode"), dict(choices=("remote", "local"), default="remote", help="remote: NCBI QBLAST, no install, minutes per search. local: blastp from flags3 install core against a local database. Default = remote")),
	(("-bd", "--blast_db"), dict(default="refseq_select", help="refseq_select, refseq_protein, genbank, swissprot, or a local database path. Default = refseq_select")),
	(("-be", "--blast_evalue"), dict(type=float, default=1e-5, help="BlastP E-value cutoff. Default = 1e-5")),
	(("-bh", "--blast_hits"), dict(type=int, default=50, help="Hits carried forward as queries; each one is a genome and a row. Default = 50")),
	(("-bw", "--blast_wait"), dict(type=float, default=60, help="Minutes to wait on NCBI's queue before giving up. Default = 60")),
	(("-u", "--user_email"), dict(help="Email address, required by NCBI Entrez")),
	(("-api", "--api_key"), dict(help="NCBI API key; raises the request rate")),
	(("-m", "--max_assemblies"), dict(type=int, default=1, help="Max assemblies per bare protein query. Default = 1")),
	(("-nc", "--no_cross_db"), dict(action="store_true", help="RefSeq proteins resolve only to GCF_ assemblies, INSDC only to GCA_")),
	(("-rm", "--remap"), dict(action="store_true", help="Ask IPG about paired queries too, so an accession renamed in its assembly still matches")),
	(("-gd", "--genomes"), dict(metavar="DIR", default=str(home.GENOMES), help="Genome cache: looked in first, downloaded into, never emptied. Default = ~/.flags3/genomes")),
	(("--offline",), dict(action="store_true", help="Resolve queries only against the genome directory; no NCBI or MGnify access")),
	(("-lf", "--local_first"), dict(action="store_true", help="Look bare queries up in the genome directory's protein files before asking IPG, so your own assemblies and NCBI queries can share one list")),
	(("--no_cache",), dict(action="store_true", help="Neither read nor fill the genome cache: genomes are downloaded into <run>/genomes/ and live with the run")),
	(("-cm", "--cluster_method"), dict(default="jackhmmer", metavar="NAME", help="Row of the tools table that finds homologous flanking proteins. Default = jackhmmer")),
	(("-cr", "--cluster_rna"), dict(action="store_true", help="Cluster flanking RNA genes too, with the rna_method row (nhmmer); RNAs without a sequence are grouped by product name")),
	(("--rna_method",), dict(default="nhmmer", metavar="NAME", help="Row of the tools table used for RNA clustering. Default = nhmmer")),
	(("-t", "--tree"), dict(action="store_true", help="Build a phylogenetic tree of the query proteins (mafft, trimal, VeryFastTree)")),
	(("-iq", "--iqtree"), dict(action="store_true", help="Build the tree with IQ-TREE (ModelFinder, 1000 ultrafast bootstraps) instead of VeryFastTree; implies --tree")),
	(("-to", "--tree_order"), dict(action="store_true", help="Order rows in figures and tables by tree leaf order; implies --tree")),
	(("-tm", "--trimal_mode"), dict(default="gt", help="trimal column filter: gt, cons, st, or a preset such as gappyout, strict, automated1. Default = gt")),
	(("-tv", "--trimal_value"), dict(type=float, default=0.1, help="Value for a trimal mode that takes one. Default = 0.1")),
	(("-tx", "--trimal_extra"), dict(default="", help="Extra trimal arguments, verbatim")),
	(("-d", "--domains"), dict(action="store_true", help="Scan flanking proteins for domains (pyhmmer) and draw them as wedges")),
	(("-db", "--hmmdb"), dict(action="append", metavar="[NAME=]PATH", help="HMM database: a .hmm file or a directory of .hmm files (e.g. DefenseFinder profiles). Repeat for several. Default = ~/.flags3/db/pfam/Pfam-A.hmm")),
	(("-hc", "--hmm_coverage"), dict(action="append", metavar="[NAME=]Q[,H]", help="Minimum fraction of the protein (Q) and of the model (H) a hit must span; NAME= limits it to one database")),
	(("-e", "--ethreshold"), dict(type=float, default=1e-3, help="Domain E-value for models without a gathering threshold. Default = 1e-3")),
	(("-cl", "--clans"), dict(metavar="TSV", help="Pfam-A.clans.tsv(.gz): colour domains by clan instead of family")),
	(("-ip", "--interpro"), dict(metavar="TSV", help="InterPro metadata table joined onto the domain table. Default: ~/.flags3/db/interpro/ if present")),
	(("-ss", "--sismis"), dict(action="store_true", help="Scan the genomic windows for secretion systems with Sismis and draw them as bands")),
	(("-gn", "--genomad"), dict(action="store_true", help="Scan the genomic windows for proviruses and plasmids with geNomad and draw them as bands")),
	(("-gdb", "--genomad_db"), dict(metavar="DIR", help="geNomad database, overriding the db option of the genomad tool row")),
	(("-df", "--defensefinder"), dict(action="store_true", help="Call anti-phage defence systems on the neighbourhood genes with DefenseFinder")),
	(("-pl", "--padloc"), dict(action="store_true", help="Call anti-phage defence systems with PadLoc; with -df, a system both agree on is drawn once")),
	(("-th", "--tmhmm"), dict(action="store_true", help="Predict transmembrane regions with DeepTMHMM (BioLib cloud unless -lth)")),
	(("--tmhmm_app",), dict(default="DTU/DeepTMHMM2", metavar="NAME", help="BioLib app for -th. Default = DTU/DeepTMHMM2; DTU/DeepTMHMM for the 1.0 model")),
	(("-lth", "--local_tmhmm"), dict(action="store_true", help="Run DeepTMHMM from the tools table instead of BioLib; implies -th")),
	(("-sp", "--signalp"), dict(action="store_true", help="Predict signal peptides with SignalP 6 (BioLib cloud unless -lsp)")),
	(("-lsp", "--local_signalp"), dict(action="store_true", help="Run SignalP from the tools table instead of BioLib; implies -sp")),
	(("-f", "--figures"), dict(metavar="TSV", help="Figure table: which figures to draw and how. Default: the shipped visualisation_table.tsv")),
	(("-nf", "--no_figures"), dict(action="store_true", help="Write the tables but draw nothing; flags3 figures <run> draws later")),
	(("-fh", "--figure_height"), dict(type=int, default=16383, metavar="PX", help="Split a figure into parts above this height. Default = 16383")),
	(("-no", "--no_overlaps"), dict(action="store_true", help="Classic mode: leave the family number out of a gene too small to hold it")),
	(("-pdf", "--pdf"), dict(action="store_true", help="Also write a PDF beside every figure (needs cairosvg)")),
	(("-c", "--cpu"), dict(type=int, help="Parallel workers (default: auto)")),
	(("--tools",), dict(metavar="TSV", help="Tool table overriding the shipped defaults. Default: ~/.flags3/tools_table.tsv if present")),
)


def build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(prog="flags3", description="Flanking gene analysis")
	parser.add_argument("-v", "--version", action="version", version="FlaGs3 " + VERSION)
	sub = parser.add_subparsers(dest="command", required=True)

	run = sub.add_parser("run", help="Run the pipeline into a new directory")
	run.add_argument("-i", "--input_list", action="append", metavar="FILE", default=[], help="Protein accessions, one per line, optionally with a tab and an assembly; a tab and BLAST marks a line as a BlastP starting point. May be repeated")
	run.add_argument("-o", "--output", default="output")
	run.add_argument("-nt", "--no_timestamp", action="store_true")
	for flags, kwargs in RUN_OPTIONS:
		run.add_argument(*flags, **kwargs)
	_common(run)

	install = sub.add_parser("install", help="Install tools and databases under ~/.flags3; no arguments lists them")
	install.add_argument("components", nargs="*", metavar="COMPONENT")
	install.add_argument("path", nargs="?", metavar="PATH", help="Archive, directory or file for a component that takes one (pfam, interpro, defence-hmm, signalp, deeptmhmm)")
	install.add_argument("--all", action="store_true", help="Every downloadable component")
	install.add_argument("--update", action="store_true", help="Refresh databases and environments of installed components")
	install.add_argument("--force", action="store_true", help="Reinstall a component even if it is already installed")

	for name in BY_NAME:
		stage = sub.add_parser(name, help="Rerun the {} stage on an existing run directory".format(name))
		stage.add_argument("run_dir")
		for flags, kwargs in RUN_OPTIONS:
			kwargs = dict(kwargs)
			kwargs.pop("default", None)
			if kwargs.get("action") == "store_true":
				kwargs["default"] = None
				stage.add_argument(*flags, **kwargs)
				long = flags[-1].lstrip("-")
				stage.add_argument("--no-" + long, dest=long, action="store_false", default=None,
					help="Switch off {} for this rerun even though the run had it on".format(flags[-1]))
				continue
			stage.add_argument(*flags, **kwargs)
		_common(stage)
	return parser


def _common(parser):
	parser.add_argument("-vb", "--verbose", action="store_true")
	parser.add_argument("-dbg", "--debug", action="store_true")


def main(argv=None) -> int:
	args = build_parser().parse_args(argv)
	console.capture()
	atexit.register(console.close)
	console.verbose = getattr(args, "verbose", False) or getattr(args, "debug", False)
	console.debug_on = getattr(args, "debug", False)
	try:
		if args.command == "install":
			return _install(args)
		if args.command == "run":
			return _run(args)
		return _stage(args)
	except (RunError, StageError, InstallError) as error:
		sys.exit("Error: {}".format(error))
	except KeyboardInterrupt:
		sys.exit("\nInterrupted.")


def _install(args) -> int:
	from flags3 import install
	components = list(args.components)
	if args.path is None and len(components) > 1 and Path(components[-1]).exists() and components[-1] not in install.BY_NAME:
		args.path = components.pop()
	args.components = components
	return install.main(args)


def _run(args) -> int:
	target = Path(args.output)
	if not args.no_timestamp:
		target = target.with_name(target.name + time.strftime("_%Y%m%d_%H%M%S"))
	if not args.input_list and not args.blast_input:
		sys.exit("Error: give -i/--input_list, -bi/--blast_input, or both.")
	run = RunDir(target).create(VERSION, " ".join(sys.argv))
	console.attach(run.console_log)
	names, inline = [], False
	for source in args.input_list:
		copied = run.input_dir / Path(source).name
		shutil.copy(source, copied)
		names.append(copied.name)
		inline = inline or bool(InputList().read(copied).blast)
	config = run.config()
	config.set("inputs", names)
	config.set("blast_inline", inline)
	if args.no_cache:
		args.genomes = str(run.path / "genomes")
	if args.local_tmhmm:
		args.tmhmm = True
	if args.local_signalp:
		args.signalp = True
	if args.iqtree or args.tree_order:
		args.tree = True
	for flags, _ in RUN_OPTIONS:
		key = flags[-1].lstrip("-")
		config.set(key, _resolved(key, getattr(args, key)))
	config.save()
	stages = [stage() for stage in PIPELINE]
	config.set("fetch.slots", sorted({slot for stage in stages if stage.wanted(config) for slot in stage.needs}))
	config.save()
	results = Runner(run, config).execute_all(stages)
	print("run directory: {}".format(run.path))
	return 0 if all(results.values()) else 1


def _stage(args) -> int:
	run = RunDir(args.run_dir).require()
	config = run.config()
	for flags, _ in RUN_OPTIONS:
		key = flags[-1].lstrip("-")
		value = getattr(args, key)
		if value is not None:
			config.set(key, _resolved(key, value))
	return 0 if Runner(run, config).execute(BY_NAME[args.command]()) else 1


PATH_KEYS = ("genomes", "tools", "clans", "interpro", "genomad_db", "figures", "blast_input")


def _resolved(key, value):
	if key == "tools" and not value:
		found = home.user_tools_table()
		return str(found) if found else value
	if key in PATH_KEYS and value:
		return str(Path(value).expanduser().resolve())
	return value


if __name__ == "__main__":
	sys.exit(main())
