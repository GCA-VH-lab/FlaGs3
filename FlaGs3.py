import argparse
import atexit
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, NamedTuple

import requests
from urllib3.util.retry import Retry
from Bio import Entrez

import flags_log
from flags_log import debug, set_debug
import flags_blast
import flags_domains
import flags_scan
from flags_domains import DEFAULT_INTERPRO
from flags_report import (ReportWriter, VERSION, note_skipped, plural,
						  print_summary, write_run_info, write_window_tables)
import flags_cluster
from flags_extract import NeighborhoodExtractor
from flags_fetch import (AssemblyDownloader, GenomeFiles,
						 LocalGenomeResolver, MgnifyGenomeDownloader,
						 ProteinAssemblyMapper)


DEFAULT_HMMDB = "./pfam_db/Pfam-A.hmm"


import pyhmmer


class AccessionListReader: 
	def __init__(self, paths):
		self.paths = [p for p in (paths or [])]

	BLAST_MARK = "blast"

	def read(self):
		proteins_assembly, proteins_only, blast_queries = [], [], []
		for path in self.paths:
			with open(path, 'r') as f:
				for n, line in enumerate(f, 1):
					line = line.split('#', 1)[0].strip()
					if not line:
						continue
					if '\t' in line:
						fields = [c for c in line.split('\t') if c]
						if len(fields) >= 2:
							if fields[1].strip().lower() == self.BLAST_MARK:
								blast_queries.append((fields[0], path, n))
							else:
								proteins_assembly.append([fields[0], fields[1]])
						else:
							print("Warning: line {} of {} is malformed, skipping it: "
								  "{!r}".format(n, path, line))
					else:
						proteins_only.append(line)
		seen = set()
		proteins_only = [p for p in proteins_only
						 if not (p in seen or seen.add(p))]
		return proteins_assembly, proteins_only, blast_queries


class InstanceLock:

	def __init__(self, path: str):
		self.path = path
		self.acquired = False

	def acquire(self) -> Optional[str]:
		for attempt in (1, 2):
			try:
				fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
			except FileExistsError:
				holder = self._read()
				if attempt == 1 and self._is_stale(holder):
					debug("removing stale lock {} ({})".format(self.path, holder))
					try:
						os.unlink(self.path)
						continue
					except OSError:
						pass
				return holder
			except OSError as e:
				debug("could not create lock {}: {!r}".format(self.path, e))
				return None
			with os.fdopen(fd, "w") as out:
				out.write("pid {}\nhost {}\nstarted {}\n".format(
					os.getpid(), platform.node(),
					time.strftime("%Y-%m-%d %H:%M:%S")))
			self.acquired = True
			debug("acquired lock {}".format(self.path))
			return None
		return self._read()

	def release(self):
		if not self.acquired:
			return
		try:
			os.unlink(self.path)
			debug("released lock {}".format(self.path))
		except OSError as e:
			debug("could not remove lock {}: {!r}".format(self.path, e))
		self.acquired = False

	def _read(self) -> str:
		try:
			with open(self.path) as fh:
				return " / ".join(fh.read().split("\n")).strip(" /")
		except OSError:
			return "unknown"

	def _is_stale(self, holder: str) -> bool:
		match = re.search(r"pid (\d+)", holder)
		host = re.search(r"host (\S+)", holder)
		if not match or (host and host.group(1) != platform.node()):
			return False
		try:
			os.kill(int(match.group(1)), 0)
		except ProcessLookupError:
			return True
		except PermissionError:
			return False
		except OSError:
			return False
		return False


def build_parser():
	usage = ''' Description: Identify flanking genes and cluster them based on similarity. Requirement= Python3, BioPython, pyhmmer, requests. '''
	parser = argparse.ArgumentParser(description=usage)
	parser.add_argument("-i", "--input_list", action="append", metavar="FILE", help=" Protein Accession eg. WP_047256880.1, optionally tab-separated with an assembly Identifier eg. GCF_000001765.3 (NCBI RefSeq/GenBank) or MGYG000454827 (MGnify Genomes catalogue -- the protein accession must then be the exact locus tag used in that genome's annotation). One per line. May be given more than once. Optional if --blast_input or a legacy list is given. ")
	parser.add_argument("-bi", "--blast_input", help=" File holding ONE starting point: a RefSeq protein accession (WP_/NP_/YP_/XP_/AP_), or a protein sequence as FASTA or bare residues over any number of lines. BlastP finds its homologues and they become the queries, appended to anything in -i. ")
	parser.add_argument("-bm", "--blast_mode", choices=("remote", "local"), default="remote", help=" Where to run BlastP. 'remote' uses NCBI QBLAST -- no install, but a search takes minutes. 'local' runs blastp from NCBI BLAST+ against a local database, which is far faster but needs the binary and the database. Default = remote ")
	parser.add_argument("-bd", "--blast_db", default="refseq_select", help=" Database to search. Aliases: refseq_select (representative RefSeq proteins), refseq_protein (full RefSeq), genbank (nr), swissprot. Each resolves to the right name for the chosen mode. Any other value passes through unchanged, which is how a local database path is given. Default = refseq_select ")
	parser.add_argument("-be", "--blast_evalue", type=float, default=1e-5, help=" E-value cutoff for BlastP. Default = 1e-5 ")
	parser.add_argument("-bw", "--blast_wait", type=float, default=60, help=" Give up waiting on NCBI's queue after this many minutes. The job usually still finishes on NCBI's side and the printed link stays valid. Default = 60 ")
	parser.add_argument("-bh", "--blast_hits", type=int, default=50, help=" How many BlastP hits to carry forward as queries. Each one becomes a genome download and a row in the figure, so this is the main control on how big the run gets. No upper limit, but webFlaGs caps at 200 and remote QBLAST may return fewer than asked. Default = 50 ")
	parser.add_argument("-u", "--user_email", required=True, help=" User Email Address (required by NCBI Entrez). ")
	parser.add_argument("-api", "--api_key", help=" NCBI API Key. ")
	parser.add_argument("-g", "--gene", type=int, default=4, help=" Number of flanking genes up/downstream. Default = 4 ")
	parser.add_argument("-r", "--range", type=int, metavar="BP", help=" Take the neighbourhood as every gene within this many bases of the query gene, instead of a fixed number of genes. Measured outwards from the query gene's own start and end, and a gene straddling the edge is included. Overrides -g. Contigs shorter than the window are reported in <dir>_rangeReport.tsv rather than silently truncated. ")
	parser.add_argument("-sm", "--scan_margin", type=int, default=10000, metavar="BP", help=" Extra sequence handed to the scanning tools beyond the analysis range, so a system straddling the edge is still called whole rather than cut in half. Hits reaching into the margin are kept and marked 'partial' in the coverage column. Set 0 to scan exactly the analysis range. Default = 10000 ")
	parser.add_argument("-sr", "--scan_range", type=int, metavar="BP", help=" Genomic span around the query gene given to the tools that scan sequence rather than genes (Sismis, geNomad). Without it they get the whole genome, as before. Independent of -g/-r, so -g 5 -sr 50000 clusters 11 genes and scans 100 kb. A tool can override it in the scan_range column of tools_table.tsv. ")
	parser.add_argument("-m", "--max_assemblies", type=int, default=1, help=" Max assemblies per protein. Default = 1 ")
	parser.add_argument("-rm", "--remap", action="store_true", help=" Look up again, through IPG, any protein whose paired assembly in the input produced nothing. A protein given with an assembly otherwise never goes near IPG, so a withdrawn assembly, a failed download, or an accession that is simply not in that assembly leaves it unresolved with no second attempt. Off by default because it costs extra NCBI requests and downloads; worth it when the assembly column was compiled some time ago. ")
	parser.add_argument("-nc", "--no_cross_db", action="store_true", help=" Keep protein and genome in the same database: RefSeq proteins (WP_, NP_, ...) resolve only to GCF_ assemblies and INSDC proteins only to GCA_. Proteins whose only assemblies are in the other database are then reported as unresolved. Does not affect assemblies given explicitly in the input file. ")
	parser.add_argument("-e", "--ethreshold", type=float, default=1e-3, help=" Jackhmmer inclusion E-value threshold for clustering flanking genes. Default = 1e-3 ")
	parser.add_argument("-n", "--number", type=int, default=3, help=" Number of jackhmmer iterations for clustering. Default = 3 ")
	parser.add_argument("-cm", "--cluster_method", default="jackhmmer", metavar="NAME", help=" Which row of tools_table.tsv decides that two flanking proteins are homologous. jackhmmer is the reference and the default. mmseqs_cluster is far faster and less sensitive. Run with an unknown name to see what the table offers. ")
	parser.add_argument("-c", "--cpu", type=int, help=" Max parallel CPU workers (default: auto-detect). ")
	parser.add_argument("-tmp", "--temporary", default="./genomes", help=" Temporary directory for downloaded assemblies; deleted at the end. Default = ./genomes ")
	parser.add_argument("-o", "--output", default="output", help=" Directory for result files; its name is also the file prefix. A YYYYMMDD_HHMMSS stamp of the run start is appended, so repeated runs do not overwrite each other. Default = output ")
	parser.add_argument("-nt", "--no_timestamp", action="store_true", help=" Use -O/--output verbatim instead of appending a date-time stamp. Repeated runs then overwrite each other; useful for scripted pipelines that need a fixed path. ")
	parser.add_argument("-t", "--tree", action="store_true", help=" Also build a phylogenetic tree and write it as <dir>_tree.nwk. Figures that use it are controlled by the figure table. ")
	parser.add_argument("-tm", "--trimal_mode", default="gt", help=" trimal column-filter mode: gt, cons, st, or a preset such as gappyout, strict, strictplus, automated1, nogaps, noallgaps. Default = gt ")
	parser.add_argument("-tv", "--trimal_value", type=float, default=0.1, help=" Value for the trimal mode that takes one (gt, cons, st). Default = 0.1 ")
	parser.add_argument("-tx", "--trimal_extra", default="", help=" Extra trimal arguments, passed through verbatim, e.g. \"-w 3\". ")
	parser.add_argument("-iq", "--iqtree", action="store_true", help=" Build the tree with IQ-TREE (ModelFinder plus 1000 ultrafast bootstrap replicates) instead of VeryFastTree. Implies --tree. Much slower but gives model selection and branch support. Needs iqtree on PATH. ")
	parser.add_argument("-to", "--tree_order", action="store_true", help=" Order the neighbours output by tree leaf order (implies --tree). ")
	parser.add_argument("-d", "--domains", action="store_true", help=" Scan flanking proteins for domains and write <dir>_domains.tsv (requires --hmmdb). ")
	parser.add_argument("-db", "--hmmdb", action="append", metavar="[NAME=]PATH", help=" HMM database for domain scanning: a .hmm file, or a directory of .hmm files such as DefenseFinder's profiles/. Repeat for several. Prefix with NAME= to label it in the outputs, otherwise the file or directory name is used. Models carrying a gathering threshold are scored by it; the rest use --ethreshold. Default = ./pfam_db/Pfam-A.hmm ")
	parser.add_argument("-hc", "--hmm_coverage", action="append", metavar="[NAME=]Q[,H]", help=" Minimum fraction of the protein (Q) and of the model (H) an alignment must span, dropping partial hits. Give NAME= to apply it to one database, or omit NAME to apply it to all. Sensible for full-length protein models such as DefenseFinder (e.g. 0.7,0.5); leave off for domain databases like Pfam, where partial coverage is normal. ")
	parser.add_argument("-pdf", "--pdf", action="store_true", help=" Also write a PDF beside every figure. Needs one of cairosvg, svglib, rsvg-convert or inkscape; without one the figures are still written as SVG. ")
	parser.add_argument("-nf", "--no_figures", action="store_true", help=" Run the analysis and write the tables, but draw nothing. Figures can be produced later with flags_redraw.py. ")
	parser.add_argument("-fh", "--figure_height", type=int, default=16383, metavar="PX", help=" Split a figure into parts when it would be taller than this, writing <name>_part1.svg and so on. The default is the canvas limit of Illustrator and librsvg, past which a figure opens in a browser but nowhere else. Raise it if you only ever view figures in a browser. Figures carrying a tree panel are never split, since the tree spans every row. ")
	parser.add_argument("-no", "--no_overlaps", action="store_true", help=" Leave a family number out of a gene too small to hold it. By default the number is drawn anyway, since knowing which family a small gene belongs to usually matters more than the overlap. ")
	parser.add_argument("-f", "--figures", metavar="TSV", help=" Figure table controlling which figures are drawn and every parameter of how. Default: visualisation_table.tsv, written into the output directory for you to edit and re-apply with flags_redraw.py. ")
	parser.add_argument("-cl", "--clans", help=" Pfam-A.clans.tsv(.gz): colour domains by clan instead of family. ")
	parser.add_argument("-ip", "--interpro", default=DEFAULT_INTERPRO, help=" InterPro metadata table (.tsv or .tsv.gz) with 'accession' and 'pfam_members' columns. Adds the InterPro entry, its name and type, and short characterisation/informativeness summaries to _domains.tsv, joined on the Pfam accession. Looked for in the working directory and next to FlaGs3.py; if it is not there the domain table is written without those columns. Default = " + DEFAULT_INTERPRO + " ")
	parser.add_argument("-lth", "--local_tmhmm", action="store_true", help=" Predict transmembrane regions with DeepTMHMM on this machine rather than in the BioLib cloud. Implies --tmhmm. The command and its directory come from the deeptmhmm row of tools_table.tsv; run deeptmhmm_installer.sh to set that up. ")
	parser.add_argument("-lsp", "--local_signalp", action="store_true", help=" Predict signal peptides with SignalP on this machine rather than in the BioLib cloud. Implies --signalp. The command and its directory come from the signalp row of tools_table.tsv; run signalp_installer.sh to set that up. ")
	parser.add_argument("--tools", metavar="TSV", help=" Table of external tool commands (mafft, trimal, VeryFastTree, IQ-TREE, blastp, sismis, DeepTMHMM, SignalP). Default: the tools_table.tsv next to FlaGs3.py. ")
	parser.add_argument("-th", "--tmhmm", action="store_true", help=" Predict transmembrane regions (DeepTMHMM, via BioLib cloud) and draw them as double red dotted lines in the domain figure. Off by default; needs pybiolib and network. ")
	parser.add_argument("-sp", "--signalp", action="store_true", help=" Predict signal peptides (SignalP-6, via BioLib cloud) and draw them as black triangles in the domain figure. Off by default; needs pybiolib and network. ")
	parser.add_argument("-df", "--defensefinder", action="store_true", help=" Call anti-phage defence systems with DefenseFinder and draw each one as a band spanning its genes, labelled with the system name. Runs on the neighbourhood genes, so it follows -g/-r and ignores -sr. Needs defensefinder_installer.sh. ")
	parser.add_argument("-pl", "--padloc", action="store_true", help=" Call anti-phage defence systems with PadLoc, drawn the same way as DefenseFinder's. Both can run together; a system both tools agree on is drawn once and the called_by column names both. Needs padloc_installer.sh. ")
	parser.add_argument("-gn", "--genomad", action="store_true", help=" Scan for proviruses and plasmids with geNomad and draw them as bands like Sismis' secretion systems. Each band is labelled with what geNomad actually called it -- the virus taxon, or plasmid/conjugative plasmid -- not a generic 'mobile element'. Uses -sr like Sismis does, so it sees a window around each query rather than whole genomes. Needs genomad and its database; run genomad_installer.sh then genomad_loader.sh. ")
	parser.add_argument("-gdb", "--genomad_db", metavar="DIR", help=" geNomad database directory, overriding the db column of the genomad row in tools_table.tsv. ")
	parser.add_argument("-ss", "--sismis", action="store_true", help=" Scan each assembly's genomic FASTA for secretion systems using Sismis (github.com/lmc297/Sismis) and write <dir>_secretion.tsv, noting which query neighborhoods (if any) each hit overlaps. Downloads the genomic FASTA per assembly. Off by default; needs sismis installed (pip install sismis). ")
	parser.add_argument("-k", "--keep", action="store_true", help=" Keep the downloaded assemblies instead of deleting the temporary directory at the end. ")
	parser.add_argument("-ul", "--use_local", metavar="DIR", help=" Directory of local .gff/.faa genome files to search before falling back to NCBI. Files may be gzipped; a genome is a .gff and .faa sharing a basename. ")
	parser.add_argument("-cr", "--cluster_rna", action="store_true", help=" Cluster flanking RNA genes into families too (off by default). Uses RNA sequences (nhmmer) when available, otherwise groups by product name. ")
	parser.add_argument("-vb", "--verbose", action="store_true", help=" Print progress and a stage-by-stage funnel. ")
	parser.add_argument("-dbg", "--debug", action="store_true", help=" Print diagnostics to stderr: HTTP status codes and Retry-After headers, external command lines and exit codes, and full tracebacks for errors that are otherwise only summarised. Implies --verbose. ")
	parser.add_argument("-nl", "--no_lock", action="store_true", help=" Skip the lock that stops two FlaGs3 runs sharing one -tmp directory. Two runs writing the same genome files also double the request rate against NCBI and EBI, which gets you throttled, so only use this with separate -tmp directories. ")
	parser.add_argument("-v", "--version", action="version", version="FlaGs3 " + VERSION)
	return parser


class Scans(NamedTuple):
	module: object
	hits: list
	rows: dict
	statuses: dict
	features: dict
	scanner: object = None
	genomad: object = None
	genomad_hits: list = ()
	genomad_statuses: dict = {}
	genomad_scanner: object = None
	defence: object = None
	defence_hits: list = ()
	defence_statuses: dict = {}
	defence_replicons: int = 0


def scan_outcome(statuses):
	ok, failed, skipped = 0, [], []
	for assembly, status in statuses.items():
		if status.startswith("error"):
			failed.append((assembly, status))
		elif status.startswith("skipped"):
			skipped.append((assembly, status))
		else:
			ok += 1
	return ok, failed, skipped


def report_scan(label, statuses, hits, noun, verbose):
	ok, failed, skipped = scan_outcome(statuses)
	if failed and not ok:
		print("Warning: {} failed on every genome. First error: {}".format(
			label, failed[0][1]))
		print("         See <prefix>_{}_diagnostics.txt, and the "
			  "'--- {} (exit N) ---' block in <prefix>_console.log for the "
			  "tool's own message.".format(label.lower(), label.lower()))
	elif failed:
		print("Warning: {} failed on {} of {}; the rest were scanned. "
			  "See <prefix>_{}_diagnostics.txt.".format(
				  label, len(failed), plural(len(statuses), "genome"),
				  label.lower()))
	if verbose and ok:
		print(">> {}: scanned {}, found {}{}".format(
			label, plural(ok, "genome"), plural(len(hits), noun),
			"" if not skipped else
			" ({} skipped)".format(len(skipped))), flush=True)


def run_background_scans(args, extractor, downloaded, all_neighborhoods, timings):
	genomad_box = {"mod": None, "hits": [], "statuses": {}, "scanner": None}
	defence_box = {"mod": None, "hits": [], "statuses": {}, "replicons": 0}
	sismis_scanner = []
	sismis_mod = None
	sismis_hits: List = []
	sismis_rows: Dict[str, Tuple[str, str, int, int]] = {}
	sismis_statuses: Dict[str, str] = {}
	features = {}

	def _run_sismis():
		nonlocal sismis_mod
		t = time.perf_counter()
		rows, hits, statuses = {}, [], {}
		if not (args.sismis and all_neighborhoods):
			return hits, rows, statuses, 0.0
		try:
			import flags_secretion as mod
		except ImportError:
			print("Warning: --sismis needs the sismis package (pip install sismis); skipping secretion-system detection.")
			return hits, rows, statuses, time.perf_counter() - t
		for g in all_neighborhoods:
			assembly = g.query.rsplit("|", 1)[-1]
			if g.query in rows:
				_, contig, lo, hi = rows[g.query]
				rows[g.query] = (assembly, contig, min(lo, g.start), max(hi, g.end))
			else:
				rows[g.query] = (assembly, g.contig, g.start, g.end)
		span = flags_scan.flags_tools_span(args, "sismis")
		scanner = mod.SismisScanner(out_dir=os.path.join(args.output, "sismis"))
		scanner.windows_dir = os.path.join(args.output, "windows")
		scanner.window_key = span
		windows = flags_scan.scan_windows(args, extractor, mod, "sismis")
		jobs = []
		for assembly in sorted({asm for asm, _, _, _ in rows.values()}):
			genome_path = downloaded.get(assembly, GenomeFiles()).genome
			if not genome_path:
				statuses[assembly] = "skipped: no genomic FASTA downloaded"
				continue
			jobs.append((assembly, genome_path, windows.get(assembly, [])))
		if jobs:
			hits = scanner.scan(jobs, statuses)
		sismis_mod = mod
		sismis_scanner.append(scanner)
		return hits, rows, statuses, time.perf_counter() - t

	def _run_defence():
		t = time.perf_counter()
		wanted = [n for n, flag in (("defensefinder", args.defensefinder),
									("padloc", args.padloc)) if flag]
		if not (wanted and all_neighborhoods):
			return [], {}, 0, 0.0
		try:
			import flags_defence as mod
		except ImportError as e:
			print("Warning: defence system calling could not be loaded ({}); "
				  "skipping.".format(e))
			return [], {}, 0, time.perf_counter() - t
		wide = [g for row in extractor.window_genes.values() for g in row]
		replicons = mod.build_replicons(wide or all_neighborhoods)
		scanner = mod.DefenceScanner(out_dir=os.path.join(args.output, "defence"),
									 threads=args.cpu or 0)
		hits = []
		for tool in wanted:
			ok, where = scanner.available(tool)
			if not ok:
				print("Warning: --{} needs {}; skipping.".format(tool, where))
				note_skipped(args, tool, where)
				scanner.statuses[tool] = "skipped: {}".format(where)
				continue
			try:
				hits.extend(scanner.run(
					tool, replicons,
					extractor.window_sequences or extractor.sequences))
			except Exception as e:
				scanner.statuses[tool] = "error: {}".format(e)
		defence_box["mod"] = mod
		return (mod.merge_calls(hits), scanner.statuses, len(replicons),
				time.perf_counter() - t)

	def _run_genomad():
		t = time.perf_counter()
		statuses, hits = {}, []
		if not (args.genomad and all_neighborhoods):
			return hits, statuses, 0.0
		try:
			import flags_genomad as mod
		except ImportError as e:
			print("Warning: --genomad could not be loaded ({}); skipping.".format(e))
			return hits, statuses, time.perf_counter() - t
		import flags_tools
		database = args.genomad_db or flags_tools.get("genomad")[1]
		scanner = mod.GenomadScanner(
			out_dir=os.path.join(args.output, "genomad"),
			database=database or "", threads=args.cpu or 0,
			windows_dir=os.path.join(args.output, "windows"),
			window_key=flags_scan.flags_tools_span(args, "genomad"))
		ok, where = scanner.available()
		if not ok:
			print("Warning: --genomad needs geNomad and {}; skipping.".format(where))
			note_skipped(args, "genomad", where)
			return hits, statuses, time.perf_counter() - t
		windows = flags_scan.scan_windows(args, extractor, mod, "genomad")
		jobs = []
		for assembly in sorted({asm for asm, _, _, _ in flags_scan.row_spans(all_neighborhoods).values()}):
			genome_path = downloaded.get(assembly, GenomeFiles()).genome
			if not genome_path:
				statuses[assembly] = "skipped: no genomic FASTA downloaded"
				continue
			jobs.append((assembly, genome_path, windows.get(assembly, [])))
		if jobs:
			hits = scanner.scan(jobs, statuses)
			if args.verbose:
				print(">> geNomad: {} over {}".format(
					plural(getattr(scanner, "batches", 1), "invocation"),
					plural(len(jobs), "genome")), flush=True)
		genomad_box["mod"] = mod
		genomad_box["scanner"] = scanner
		return hits, statuses, time.perf_counter() - t

	def _run_tmhmm():
		t = time.perf_counter()
		if not args.tmhmm:
			return {}, 0.0
		try:
			import flags_features as feat_mod
			tm = feat_mod.TMScanner(local=args.local_tmhmm).scan(
				extractor.sequences, want_signal=not args.signalp)
		except ImportError:
			print("Warning: --tmhmm needs pybiolib (pip install pybiolib); skipping transmembrane prediction.")
			return {}, time.perf_counter() - t
		except Exception as e:
			print("Warning: DeepTMHMM did not finish, skipping transmembrane regions ({}).".format(e))
			note_skipped(args, "tmhmm", e)
			return {}, time.perf_counter() - t
		return tm, time.perf_counter() - t

	def _run_signalp():
		t = time.perf_counter()
		if not args.signalp:
			return {}, 0.0
		try:
			import flags_features as feat_mod
			sp = feat_mod.SignalPScanner(local=args.local_signalp).scan(
				extractor.sequences)
		except ImportError:
			print("Warning: --signalp needs pybiolib (pip install pybiolib); skipping signal-peptide prediction.")
			return {}, time.perf_counter() - t
		except Exception as e:
			print("Warning: SignalP did not finish, skipping signal peptides ({}).".format(e))
			note_skipped(args, "signalp", e)
			return {}, time.perf_counter() - t
		return sp, time.perf_counter() - t

	tm, sp = {}, {}
	active = [n for n, flag in (("sismis", args.sismis), ("genomad", args.genomad),
								 ("defence", args.defensefinder or args.padloc),
								 ("tmhmm", args.tmhmm),
								 ("signalp", args.signalp)) if flag]
	if active:
		if args.verbose:
			local = [n for n, on in (("tmhmm", args.local_tmhmm),
									 ("signalp", args.local_signalp))
					 if on and n in active]
			where = "locally: {}".format(", ".join(local)) if local else \
				"cloud/subprocess, not local CPU"
			print(">> running {} in the background ({})...".format(
				", ".join(active), where), flush=True)
		if ((args.tmhmm and not args.local_tmhmm)
				or (args.signalp and not args.local_signalp)):
			import flags_features as feat_mod
			feat_mod.warm_up()
		task = {"sismis": _run_sismis, "genomad": _run_genomad,
				"defence": _run_defence,
				"tmhmm": _run_tmhmm, "signalp": _run_signalp}
		with ThreadPoolExecutor(max_workers=len(task)) as pool:
			futures = {pool.submit(task[name]): name for name in active}
			for fut in as_completed(futures):
				name = futures[fut]
				if name == "sismis":
					sismis_hits, sismis_rows, sismis_statuses, elapsed = fut.result()
					timings["sismis_scan"] = elapsed
					if sismis_mod:
						report_scan("sismis", sismis_statuses, sismis_hits,
									"secretion system", args.verbose)
				elif name == "genomad":
					genomad_box["hits"], genomad_box["statuses"], elapsed = fut.result()
					timings["genomad_scan"] = elapsed
					report_scan("geNomad", genomad_box["statuses"],
								genomad_box["hits"], "mobile element", args.verbose)
				elif name == "defence":
					(defence_box["hits"], defence_box["statuses"],
					 defence_box["replicons"], elapsed) = fut.result()
					timings["defence_scan"] = elapsed
					bad = [t for t, v in defence_box["statuses"].items()
						   if v.startswith("error")]
					if bad and not defence_box["hits"]:
						print("Warning: {} produced nothing but an error. See "
							  "<prefix>_defence_diagnostics.txt and the console "
							  "log.".format(" and ".join(sorted(bad))))
					elif args.verbose:
						print(">> defence systems: {} across {}".format(
							plural(len(defence_box["hits"]), "system"),
							plural(defence_box["replicons"], "neighbourhood")), flush=True)
				elif name == "tmhmm":
					tm, elapsed = fut.result()
					timings["tmhmm_scan"] = elapsed
					if args.verbose:
						print(">> DeepTMHMM: features on {} proteins".format(len(tm)), flush=True)
				elif name == "signalp":
					sp, elapsed = fut.result()
					timings["signalp_scan"] = elapsed
					if args.verbose:
						print(">> SignalP: signal peptides on {} proteins".format(len(sp)), flush=True)
		for acc, regs in tm.items():
			features.setdefault(acc, []).extend(regs)
		for acc, regs in sp.items():
			features.setdefault(acc, []).extend(regs)

	return Scans(sismis_mod, sismis_hits, sismis_rows, sismis_statuses, features,
				 sismis_scanner[0] if sismis_scanner else None,
				 genomad_box["mod"], genomad_box["hits"],
				 genomad_box["statuses"], genomad_box["scanner"],
				 defence_box["mod"], defence_box["hits"],
				 defence_box["statuses"], defence_box["replicons"])


def main():
	flags_log.start_buffering()
	atexit.register(flags_log.close)
	parser = build_parser()
	args = parser.parse_args()

	args.input_list = list(args.input_list or [])
	all_lists = args.input_list
	if not all_lists and not args.blast_input:
		sys.exit("Error: give -i/--input_list, -bi/--blast_input, or both.")
	for path in all_lists:
		if not os.path.isfile(path):
			sys.exit("Error: input list not found: {}".format(path))
	if args.blast_input and not os.path.isfile(args.blast_input):
		sys.exit("Error: --blast_input file not found: {}".format(args.blast_input))
	if args.blast_input and args.blast_hits < 2:
		sys.exit("Error: --blast_hits must be at least 2.")
	if args.use_local and not os.path.isdir(args.use_local):
		sys.exit("Error: --use_local directory not found: {}".format(args.use_local))
	args.hmmdb = list(args.hmmdb or [])
	if args.domains and not args.hmmdb:
		args.hmmdb = [DEFAULT_HMMDB]
	for spec in args.hmmdb:
		path = spec.partition("=")[2] or spec
		if not os.path.exists(os.path.expanduser(path)):
			sys.exit("Error: --hmmdb path not found: {}".format(path))
	try:
		args.hmm_coverage = flags_domains.parse_coverage(args.hmm_coverage)
	except ValueError as e:
		sys.exit("Error: {}".format(e))
	if args.clans and not os.path.isfile(args.clans):
		sys.exit("Error: --clans file not found: {}".format(args.clans))
	if args.range is not None and args.range < 1:
		sys.exit("Error: --range must be a positive number of bases.")
	if args.scan_range is not None and args.scan_range < 1:
		sys.exit("Error: --scan_range must be a positive number of bases.")
	if args.scan_margin < 0:
		sys.exit("Error: --scan_margin cannot be negative.")
	if args.gene < 0:
		sys.exit("Error: --gene cannot be negative.")
	if False:
		try:
			pass
		except ValueError as e:
			sys.exit("Error: {}".format(e))

	if not args.no_timestamp:
		args.output = "{}_{}".format(os.path.normpath(args.output),
									 time.strftime("%Y%m%d_%H%M%S"))

	args.output = os.path.abspath(args.output)
	args.temporary = os.path.abspath(args.temporary)

	try:
		flags_log.attach(os.path.join(
			args.output,
			os.path.basename(os.path.normpath(args.output)) + "_console.log"))
	except OSError as e:
		print("Warning: could not open the console log ({}); the run continues "
			  "without one.".format(e))

	if args.local_tmhmm:
		args.tmhmm = True
	if args.local_signalp:
		args.signalp = True
	set_debug(args.debug)
	import flags_tools
	try:
		flags_tools.load(args.tools)
	except (OSError, ValueError) as e:
		sys.exit("Error: {}".format(e))
	if args.debug:
		args.verbose = True
		debug("FlaGs3 {} on {} / python {}".format(
			VERSION, platform.platform(), sys.version.split()[0]))

	lock = InstanceLock(args.temporary + ".lock")
	if not args.no_lock:
		holder = lock.acquire()
		if holder:
			sys.exit(
				"Error: another FlaGs3 run is using {}\n"
				"  holder: {}\n"
				"Two runs sharing one temporary directory overwrite each other's "
				"genome files and double the request rate against NCBI and EBI, which "
				"gets the host throttled.\n"
				"Give this run its own directory with -tmp, or pass --no_lock to "
				"override. If no such run exists, delete {}".format(
					args.temporary, holder, lock.path))
		atexit.register(lock.release)

	write_run_info(args, parser)

	timings = {}
	t_start = time.perf_counter()
	t0 = t_start

	proteins_assembly, proteins_only, inline_blast = (
		AccessionListReader(args.input_list).read() if all_lists
		else ([], [], []))
	proteins_only = list(proteins_only)
	timings["1_read_input"] = time.perf_counter() - t0; t0 = time.perf_counter()

	blast_hits, blast_queries, blast_mod = flags_blast.resolve_blast(
		args, proteins_assembly, proteins_only, inline_blast, timings, t0)
	t0 = time.perf_counter()

	all_queries = [p for p, _ in proteins_assembly] + list(proteins_only)
	if args.verbose:
		print(">> read {} queries ({} paired, {} protein-only)".format(
			len(all_queries), len(proteins_assembly), len(proteins_only)), flush=True)

	local = LocalGenomeResolver(args.use_local) if args.use_local else None
	local_files: Dict[str, GenomeFiles] = {}
	protein_to_assemblies = {}
	local_acceptable = {}
	pending_pairs = []
	pending_only = []

	if local:
		for protein, assembly in proteins_assembly:
			hit = local.resolve_pair(assembly)
			if hit:
				base, files = hit
				local_files[base] = files
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_pairs.append([protein, assembly])
		for protein in proteins_only:
			hit = local.resolve_protein(protein)
			if hit:
				base, files = hit
				local_files[base] = files
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_only.append(protein)
		if args.verbose:
			print(">> local: {} genomes indexed; resolved {} of {} queries locally".format(
				len(local.genomes), len(protein_to_assemblies),
				len(proteins_assembly) + len(proteins_only)), flush=True)
	else:
		pending_pairs = proteins_assembly
		pending_only = proteins_only

	mapper = ProteinAssemblyMapper(email=args.user_email, api_key=args.api_key,
								   max_assemblies=args.max_assemblies,
								   cross_db=not args.no_cross_db)
	ncbi_map = mapper.map(pending_only)
	for protein, asms in ncbi_map.items():
		protein_to_assemblies.setdefault(protein, []).extend(asms)
	for protein, assembly in pending_pairs:
		protein_to_assemblies.setdefault(protein, []).append(assembly)
	timings["2_ipg_mapping"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose and pending_only:
		print(">> NCBI IPG: mapped {} of {} remaining proteins to assemblies".format(
			len(ncbi_map), len(pending_only)), flush=True)
	if mapper.unreachable:
		print("Warning: could not reach NCBI to resolve {} remaining protein(s) "
			  "({}). They are reported as unresolved; the rest of the run is "
			  "unaffected.".format(len(pending_only), mapper.unreachable))
	if args.verbose and mapper.dropped_cross_db:
		lost = sum(1 for a in mapper.dropped_cross_db if not ncbi_map.get(a))
		print(">> --no_cross_db: excluded cross-database assemblies for {} proteins"
			  "{}".format(len(mapper.dropped_cross_db),
						  ", {} left with none".format(lost) if lost else ""), flush=True)

	assemblies = sorted({asm for asms in protein_to_assemblies.values()
						 for asm in asms if asm not in local_files})
	mgnify_assemblies = [a for a in assemblies if MgnifyGenomeDownloader.is_mgnify_accession(a)]
	ncbi_assemblies = [a for a in assemblies if a not in set(mgnify_assemblies)]
	dl_workers = args.cpu if args.cpu else min(max(len(assemblies), 1), 10)
	dl_rate = 10.0 if args.api_key else 5.0

	def progress(label):
		if not args.verbose:
			return None
		return lambda done, total: print(
			">> {} download: {}/{}".format(label, done, total), flush=True)

	if args.verbose and assemblies:
		print(">> downloading {}...".format(plural(len(assemblies), "genome")), flush=True)
	downloaded: Dict[str, GenomeFiles] = {}
	failures: Dict[str, str] = {}
	if ncbi_assemblies:
		dl = AssemblyDownloader(out_dir=args.temporary, workers=dl_workers, rate=dl_rate,
								want_rna=args.cluster_rna,
								want_genome=args.sismis or args.genomad)
		downloaded.update(dl.download_many(ncbi_assemblies, progress("NCBI")))
		failures.update(dl.failures)
		timings["3a_download_ncbi"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if mgnify_assemblies:
		mg = MgnifyGenomeDownloader(out_dir=args.temporary, workers=dl_workers, rate=dl_rate,
									want_genome=args.sismis or args.genomad or args.cluster_rna)
		downloaded.update(mg.download_many(mgnify_assemblies, progress("MGnify")))
		failures.update(mg.failures)
		timings["3b_download_mgnify"] = time.perf_counter() - t0; t0 = time.perf_counter()
	downloaded.update(local_files)
	if args.verbose:
		ready = sum(1 for f in downloaded.values() if f.gff and f.faa)
		print(">> genomes ready: {} of {} usable ({} NCBI, {} MGnify, {} local)".format(
			ready, len(downloaded), len(ncbi_assemblies), len(mgnify_assemblies),
			len(local_files)), flush=True)
		if failures:

			print(">> {}, first few:".format(plural(len(failures), "download error")), flush=True)
			for key, msg in list(failures.items())[:5]:
				print("     {}: {}".format(key.rsplit("/", 1)[-1] or key, msg), flush=True)

	extractor = NeighborhoodExtractor(flank=args.gene, range_bp=args.range,
									  scan_range=args.scan_range,
									  label_assembly=args.max_assemblies > 1)
	all_neighborhoods = []
	matched = set()
	plan = []
	for protein, asms in protein_to_assemblies.items():
		for asm in asms:
			files = downloaded.get(asm, GenomeFiles())
			if not (files.gff and files.faa):
				continue
			if asm in local_files:
				acceptable = local_acceptable.get(protein, {}).get(asm)
			else:
				acceptable = mapper.accessions_in.get(protein, {}).get(asm)
			plan.append((protein, asm, acceptable, files))

	grouped: Dict[str, List[int]] = {}
	for i, entry in enumerate(plan):
		grouped.setdefault(entry[1], []).append(i)
	extracted: Dict[int, list] = {}
	for asm, indexes in grouped.items():
		for i in indexes:
			protein, _, acceptable, files = plan[i]
			extracted[i] = extractor.extract(
				asm, files.gff, files.faa, protein, acceptable,
				rna_path=files.rna if args.cluster_rna else None,
				genome_path=files.genome if args.cluster_rna else None)
		extractor.forget(asm)
	for i, (protein, asm, _, _) in enumerate(plan):
		rows = extracted.get(i)
		if rows:
			matched.add(protein)
			all_neighborhoods.extend(rows)

	paired = {p for p, _ in proteins_assembly}
	stranded = sorted(p for p in paired if p not in matched)
	if stranded and not args.remap and args.verbose:
		print(">> {} whose paired assembly gave nothing; --remap looks them up "
			  "again through IPG".format(
				  plural(len(stranded), "protein")), flush=True)
	if stranded and args.remap:
		if args.verbose:
			print(">> retrying {} through IPG; their paired assembly gave "
				  "nothing".format(plural(len(stranded), "protein")), flush=True)
		retry = mapper.map(stranded)
		fresh = sorted({a for asms in retry.values() for a in asms}
					   if retry else [])
		fresh = [a for a in fresh if a not in downloaded and a not in local_files]
		if fresh:
			rdl = AssemblyDownloader(out_dir=args.temporary, workers=dl_workers,
									 rate=dl_rate, want_rna=args.cluster_rna,
									 want_genome=args.sismis or args.genomad)
			downloaded.update(rdl.download_many(
				[a for a in fresh if not MgnifyGenomeDownloader.is_mgnify_accession(a)],
				progress("IPG retry")))
			failures.update(rdl.failures)
		regained = 0
		for protein in stranded:
			for asm in retry.get(protein, []):
				files = downloaded.get(asm, GenomeFiles())
				if not (files.gff and files.faa):
					continue
				rows = extractor.extract(
					asm, files.gff, files.faa, protein,
					mapper.accessions_in.get(protein, {}).get(asm),
					rna_path=files.rna if args.cluster_rna else None,
					genome_path=files.genome if args.cluster_rna else None)
				extractor.forget(asm)
				if rows:
					matched.add(protein)
					all_neighborhoods.extend(rows)
					regained += 1
					protein_to_assemblies.setdefault(protein, []).append(asm)
					break
		timings["4b_ipg_retry"] = time.perf_counter() - t0; t0 = time.perf_counter()
		if args.verbose:
			print(">> IPG retry recovered {} of {}".format(
				regained, len(stranded)), flush=True)
	timings["4_extract_neighbors"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		print(">> extracted {} flanking-gene records; {} of {} queries matched".format(
			len(all_neighborhoods), len(matched), len(all_queries)), flush=True)
	if args.range and extractor.ranges:
		short = [r for r in extractor.ranges.values()
				 if r.up_available < args.range or r.down_available < args.range]
		if short:
			print("Note: {} of {} rows sit closer than {} bp to a contig end, so their "
				  "window is truncated; see <prefix>_rangeReport.tsv.".format(
					  len(short), len(extractor.ranges), args.range))

	scans = run_background_scans(args, extractor, downloaded, all_neighborhoods,
								 timings)
	sismis_mod, sismis_hits = scans.module, scans.hits
	sismis_rows, sismis_statuses = scans.rows, scans.statuses
	features = scans.features
	if not args.keep:
		shutil.rmtree(args.temporary, ignore_errors=True)

	if not all_neighborhoods:
		os.makedirs(args.output, exist_ok=True)
		prefix = os.path.basename(os.path.normpath(args.output))
		reporter = ReportWriter(all_neighborhoods, [], extractor.species,
								all_queries, protein_to_assemblies, matched)
		issues_path = os.path.join(args.output, prefix + "_accessionIssues.txt")
		reporter.accession_issues(issues_path)
		reporter.query_status(os.path.join(args.output, prefix + "_QueryStatus.txt"))
		sys.exit("No flanking neighbourhoods could be extracted for any query. "
				 "See {} for per-query details.".format(issues_path))

	t0 = time.perf_counter()
	cluster_input = extractor.sequences
	try:
		clusterer = flags_cluster.build(args.cluster_method, workers=args.cpu)
		families = clusterer.cluster(cluster_input)
	except flags_cluster.MethodError as e:
		sys.exit("Error: {}".format(e))

	rna_families = []
	if args.cluster_rna:
		rna = flags_cluster.build("nhmmer", workers=args.cpu)
		have_seq = set(extractor.rna_sequences)
		all_rna = set(extractor.rna_products)
		missing = all_rna - have_seq
		rna_families = rna.cluster(extractor.rna_sequences)
		if missing:
			fallback = {acc: extractor.rna_products[acc] for acc in missing}
			rna_families += rna.cluster_by_name(fallback)
			print("Warning: {} of {} flanking RNAs had no nucleotide sequence available, "
				  "so they were grouped by product name rather than by sequence."
				  .format(len(missing), len(all_rna)))
		families = families + rna_families
	timings["5_clustering"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		msg = ">> clustered {} into {}".format(
			plural(len(extractor.sequences), "flanking protein"),
			plural(len(families) - len(rna_families), "family", "families"))
		if args.cluster_rna:
			msg += "; {} into {}".format(
				plural(len(extractor.rna_products), "RNA"),
				plural(len(rna_families), "family", "families"))
		print(msg, flush=True)

	want_tree = args.tree or args.tree_order or args.iqtree
	newick = ""
	leaf_order = None
	tree_mod = None
	builder = None
	if want_tree:
		if args.verbose:
			print(">> building tree with {} ({} leaves)...".format(
				"IQ-TREE" if args.iqtree else "VeryFastTree",
				len(extractor.row_sequences)), flush=True)
		import flags_tree as tree_mod
		builder = tree_mod.TreeBuilder(
			threads=args.cpu or 0,
			engine="iqtree" if args.iqtree else "veryfasttree",
			trimal_mode=args.trimal_mode, trimal_value=args.trimal_value,
			trimal_extra=args.trimal_extra)
		newick, _ = builder.build(extractor.row_sequences)
		if newick:
			leaf_order = tree_mod.ladderized_leaf_order(newick)
	if want_tree:
		timings["6_tree"] = time.perf_counter() - t0
	t0 = time.perf_counter()
	if args.verbose and want_tree:
		print(">> tree: {}".format("built ({} leaves)".format(len(leaf_order))
			  if newick else "skipped (fewer than 3 query sequences)"), flush=True)

	os.makedirs(args.output, exist_ok=True)
	prefix = os.path.basename(os.path.normpath(args.output))
	def out_path(suffix):
		return os.path.join(args.output, prefix + suffix)
	order = leaf_order if args.tree_order else None
	tree_written = False
	if newick:
		tree_dir = os.path.join(args.output, "tree")
		os.makedirs(tree_dir, exist_ok=True)
		tree_path = lambda suffix: os.path.join(tree_dir, prefix + suffix)
		with open(tree_path("_tree.nwk"), "w") as out:
			out.write(newick + "\n")
		if builder:
			for suffix, aln in (("_alignment.aln", builder.raw_alignment),
								("_trimmed.aln", builder.alignment)):
				if not aln:
					continue
				with open(tree_path(suffix), "w") as out:
					for name, seq in aln.items():
						out.write(">{}\n{}\n".format(name, seq))
			if builder.commands:
				with open(tree_path("_commands.txt"), "w") as out:
					out.write("\n".join(builder.commands) + "\n")
				with open(out_path("_runinfo.txt"), "a") as out:
					out.write("\ntree commands\n")
					for c in builder.commands:
						out.write("  {}\n".format(c))
		tree_written = True
	timings["6b_tree_files"] = time.perf_counter() - t0; t0 = time.perf_counter()
	want_domain_fig = args.domains or features
	domains, clan_map, domain_table_written = flags_domains.scan_domains(
		args, extractor, families, all_neighborhoods, out_path, timings, t0,
		want_domain_fig)
	t0 = time.perf_counter()
	if features:
		try:
			import flags_features as feat_mod
			feat_mod.write_report(features, out_path("_features.tsv"))
		except Exception as e:
			print("Warning: could not write the feature table ({}).".format(e))
	t0 = time.perf_counter()

	if scans.defence:
		t0 = time.perf_counter()
		dmod = scans.defence
		dmatches = dmod.match_rows(scans.defence_hits, flags_scan.row_spans(all_neighborhoods))
		dmod.write_report(scans.defence_hits, dmatches, out_path("_defence.tsv"))
		dmod.write_diagnostics(scans.defence_statuses, scans.defence_replicons,
							   out_path("_defence_diagnostics.txt"))
		timings["9c_defence_report"] = time.perf_counter() - t0

	if args.genomad and scans.genomad:
		t0 = time.perf_counter()
		gmod = scans.genomad
		gmatches = gmod.match_rows(scans.genomad_hits, flags_scan.row_spans(all_neighborhoods))
		gmod.write_report(scans.genomad_hits, gmatches, out_path("_genomad.tsv"))
		gscanner = scans.genomad_scanner
		gmod.write_diagnostics(
			scans.genomad_statuses, out_path("_genomad_diagnostics.txt"),
			scanned=getattr(gscanner, "scanned_bases", None),
			windows=getattr(gscanner, "window_count", None),
			genome_size={a: sum(extractor.contig_lengths(a).values())
						 for a in scans.genomad_statuses if extractor.contig_lengths(a)})
		timings["9b_genomad_report"] = time.perf_counter() - t0

	if args.sismis and sismis_mod:
		t0 = time.perf_counter()
		matches = sismis_mod.match_rows(sismis_hits, sismis_rows)
		sismis_mod.write_report(sismis_hits, matches, out_path("_secretion.tsv"))
		scanner = scans.scanner
		sismis_mod.write_diagnostics(
			sismis_statuses, out_path("_sismis_diagnostics.txt"),
			scanned=getattr(scanner, "scanned_bases", None),
			windows=getattr(scanner, "window_count", None),
			genome_size={a: sum(extractor.contig_lengths(a).values())
						 for a in sismis_statuses if extractor.contig_lengths(a)})
		timings["9_sismis_report"] = time.perf_counter() - t0

	adjacency = dict(clusterer.adjacency)
	if args.cluster_rna:
		adjacency.update(rna.adjacency)
	reporter = ReportWriter(all_neighborhoods, families, extractor.species,
							all_queries, protein_to_assemblies, matched,
							order=order, adjacency=adjacency,
							sequences=extractor.sequences,
							row_sequences=extractor.row_sequences,
							ranges=extractor.ranges, requested=args.range)
	n_issues = reporter.write_all(out_path)
	if extractor.window_genes:
		write_window_tables(extractor, out_path)
	if blast_hits:
		blast_mod.write_report(blast_hits, blast_queries[0],
							   out_path("_blast_hits.tsv"))
		blast_mod.write_accessions(blast_hits, blast_queries,
								   out_path("_blast_accessions.txt"))
	if args.verbose:
		print(">> wrote data tables and reports ({} queries with issues)".format(n_issues),
			  flush=True)

	t0 = time.perf_counter()
	figures_written = []
	if not args.no_figures:
		redraw = os.path.join(os.path.dirname(os.path.abspath(__file__)),
							  "flags_redraw.py")
		cmd = [sys.executable, redraw, "--data", args.output, "--prefix", prefix,
			   "--max_height", str(args.figure_height)]
		if args.no_overlaps:
			cmd.append("--no_overlaps")
		if args.figures:
			cmd += ["--format", args.figures]
		if args.pdf:
			cmd.append("--pdf")
		if args.verbose:
			cmd.append("--verbose")
		debug("running: {}".format(" ".join(cmd)))
		try:
			done = subprocess.run(cmd, capture_output=True, text=True)
			debug("flags_redraw exit {}".format(done.returncode))
			flags_log.record("--- flags_redraw.py (exit {}) ---".format(done.returncode))
			flags_log.record(done.stdout)
			flags_log.record(done.stderr, "stderr")
			if done.returncode != 0:
				print("Warning: figures were not drawn ({}).".format(
					(done.stderr or done.stdout or "").strip()[:300]))
			else:
				for line in done.stdout.splitlines():
					name = line.strip()
					if name.endswith(".svg"):
						figures_written.append(name)
		except OSError as e:
			print("Warning: could not run {} ({}).".format(redraw, e))
	timings["7_visualize"] = time.perf_counter() - t0

	print_summary(args, prefix, extractor, families, rna_families, figures_written,
				  tree_written, want_tree, domain_table_written, features,
				  sismis_mod, blast_hits,
				  genomad_written=bool(scans.genomad),
				  defence_written=bool(scans.defence))
	if args.verbose:
		print("\n--- timing (seconds) ---")
		for stage in sorted(timings):
			print("  {:24s} {:8.2f}".format(stage, timings[stage]))
		print("  {:24s} {:8.2f}".format("TOTAL", time.perf_counter() - t_start))


if __name__ == '__main__':
	try:
		main()
	except FileNotFoundError as e:
		sys.exit("Error: file not found - {}".format(e))
	except KeyboardInterrupt:
		sys.exit("\nInterrupted.")
	except Exception as e:
		sys.exit("Error: {}".format(e))