import os
import platform
import shutil
import sys
import time

import flags_log
from collections import Counter

SECRET_ARGS = ("api_key",)

SECRET_FLAGS = ("--api_key",)

VERSION = "2.3.0"


def plural(n, word, plural_form=None):
	return "{} {}".format(n, word if n == 1 else (plural_form or word + "s"))


class ReportWriter: 
	def __init__(self, neighborhoods, families, species,
				 queries, protein_to_assemblies, matched,
				 order=None, adjacency=None, sequences=None, row_sequences=None,
				 ranges=None, requested=None):
		self.ranges = ranges or {}
		self.requested = requested
		self.neighborhoods = neighborhoods
		self.families = families
		self.species = species
		self.queries = queries
		self.protein_to_assemblies = protein_to_assemblies
		self.matched = matched
		self.adjacency = adjacency or {}
		self.sequences = sequences or {}
		self.row_sequences = row_sequences or {}
		rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
		queries = {g.accession for g in neighborhoods if g.offset == 0}
		self.occurrences = Counter(g.accession for g in neighborhoods)
		self.fam_of = family_numbers(families, rna_accessions, queries,
									 self.occurrences)
		self.by_query = {}
		for g in neighborhoods:
			self.by_query.setdefault(g.query, []).append(g)
		if order:
			rank = {row: i for i, row in enumerate(order)}
			self.by_query = {row: self.by_query[row] for row in
							 sorted(self.by_query, key=lambda r: rank.get(r, len(rank)))}
		self.products = {}
		for g in neighborhoods:
			self.products.setdefault(g.accession, g.product)

	def write_all(self, out_path):
		self.operon_tsv(out_path("_operon.tsv"))
		self.clusters_tsv(out_path("_clusters.tsv"))
		self.outdesc_txt(out_path("_outdesc.txt"))
		self.species_info(out_path("_speciesInfo.txt"))
		self.query_status(out_path("_QueryStatus.txt"))
		self.flankgene_report(out_path("_flankgene_Report.log"))
		self.fasta_outputs(out_path)
		if self.ranges:
			self.range_report(out_path("_rangeReport.tsv"))
		if self.adjacency:
			self.clusterhits_tsv(out_path("_clusterhits.tsv"))
		return self.accession_issues(out_path("_accessionIssues.txt"))

	@staticmethod
	def _split_row(row_id):
		query, _, assembly = row_id.partition("|")
		return query, assembly

	def operon_tsv(self, path):
		with open(path, "w") as out:
			out.write("#query\tassembly\tspecies\tfamily\tstrand\toffset\t"
					  "start\tend\tlength\tcontig\tis_rna\taccession\tproduct\n")
			for row_id in self.by_query:
				query, assembly = self._split_row(row_id)
				sp = self.species.get(row_id, "")
				for g in sorted(self.by_query[row_id], key=lambda x: x.offset):
					out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
						query, assembly or "-", sp, self.fam_of.get(g.accession, "-"),
						g.strand, g.offset, g.start, g.end, g.end - g.start + 1,
						g.contig or "-", "True" if g.is_rna else "False",
						g.accession, g.product))

	def clusters_tsv(self, path):
		with open(path, "w") as out:
			out.write("#family\tsize\tmembers\n")
			for fam in self.families:
				label = (self.fam_of.get(fam[0], "-")
						 if family_shared(fam, self.occurrences) else "-")
				out.write("{}\t{}\t{}\n".format(label, len(fam), ",".join(fam)))

	def outdesc_txt(self, path):
		blocks = [fam for fam in self.families
				  if family_shared(fam, self.occurrences)]
		blocks.sort(key=lambda fam: (-sum(self.occurrences.get(a, 0) for a in fam),
									 fam[0]))
		with open(path, "w") as out:
			for fam in blocks:
				label = self.fam_of.get(fam[0], "-")
				for acc in fam:
					out.write("{}({})\t{}\t{}\n".format(
						label, self.occurrences.get(acc, 0), acc,
						self.products.get(acc, "")))
				out.write("\n\n")

	def fasta_outputs(self, out_path):
		if not (self.sequences or self.row_sequences):
			return
		query_accs = {g.accession for g in self.neighborhoods if g.offset == 0}
		flanking = [a for a in sorted(self.sequences) if a not in query_accs]

		def label(acc):
			product = self.products.get(acc, "")
			return "{}|{}".format(acc, product) if product else acc

		rows_out = [(row, self.row_sequences[row]) for row in sorted(self.row_sequences)]
		flank_out = [(label(a), self.sequences[a]) for a in sorted(flanking)]
		self._write_fasta(out_path("_tree.fasta"), rows_out)
		self._write_fasta(out_path("_flankgene.fasta"), flank_out)
		self._write_fasta(out_path("_all.fasta"), rows_out + flank_out)

	def range_report(self, path):
		want = self.requested
		with open(path, "w") as out:
			out.write("#query\tassembly\tcontig\tcontig_length\tquery_start\t"
					  "query_end\tquery_strand\trequested_bp\tup_available\t"
					  "down_available\tup_reached\tdown_reached\ttruncated\t"
					  "genes_up\tgenes_down\tgenes_total\tscan_start\tscan_end\t"
					  "scan_span\n")
			for row_id in self.by_query:
				info = self.ranges.get(row_id)
				if info is None:
					continue
				sides = []
				if want:
					if info.up_available < want:
						sides.append("up")
					if info.down_available < want:
						sides.append("down")
				out.write("\t".join(str(v) for v in (
					info.query, info.assembly or "-", info.contig or "-",
					info.contig_length, info.q_start, info.q_end, info.q_strand,
					want if want else "-",
					info.up_available, info.down_available,
					info.up_reached, info.down_reached,
					",".join(sides) if sides else "-",
					info.genes_up, info.genes_down,
					info.genes_up + info.genes_down + 1,
					info.scan_start or "-", info.scan_end or "-",
					(info.scan_end - info.scan_start + 1) if info.scan_start else "-",
					)) + "\n")

	@staticmethod
	def _write_fasta(path, records):
		with open(path, "w") as out:
			for name, seq in records:
				out.write(">{}\n{}\n".format(name, seq))

	def clusterhits_tsv(self, path):
		with open(path, "w") as out:
			out.write("#accession\tfamily\tn_hits\thits\n")
			for acc in sorted(self.adjacency):
				hits = sorted(self.adjacency[acc])
				out.write("{}\t{}\t{}\t{}\n".format(
					acc, self.fam_of.get(acc, "-"), len(hits), ";".join(hits)))

	def species_info(self, path):
		with open(path, "w") as out:
			out.write("#query\tassembly\tspecies\n")
			for row_id in sorted(self.species):
				query, assembly = self._split_row(row_id)
				out.write("{}\t{}\t{}\n".format(query, assembly or "-",
											    self.species[row_id]))

	def query_status(self, path):
		with open(path, "w") as out:
			out.write("#query\tassemblies\tflanking_genes_found\n")
			for q in self.queries:
				asms = self.protein_to_assemblies.get(q, [])
				status = "Yes" if q in self.matched else "No"
				out.write("{}\t{}\t{}\n".format(q, ";".join(asms) if asms else "-", status))

	def flankgene_report(self, path):
		with open(path, "w") as out:
			for query in self.by_query:
				genes = sorted(self.by_query[query], key=lambda x: x.offset)
				chain = " ".join("{}({})".format(g.accession, self.fam_of.get(g.accession, "-"))
								 for g in genes)
				out.write("{}\t{}\n".format(query, chain))

	def accession_issues(self, path):
		lines = []
		for q in self.queries:
			if not self.protein_to_assemblies.get(q):
				lines.append("{}\tno assembly resolved (not found locally or via NCBI)".format(q))
			elif q not in self.matched:
				lines.append("{}\tassembly resolved but no flanking neighborhood extracted".format(q))
		with open(path, "w") as out:
			out.write("#query\tissue\n")
			for line in lines:
				out.write(line + "\n")
		return len(lines)


def family_shared(fam, occurrences=None):
	if occurrences:
		return sum(occurrences.get(acc, 0) for acc in fam) > 1
	return len(fam) > 1


def family_numbers(families, rna_accessions=None, query_accessions=None,
				   occurrences=None):
	rna_accessions = rna_accessions or set()
	query_accessions = query_accessions or set()
	number = {}
	prot_n, rna_n, query_n = 0, 0, 0
	shared = [fam for fam in families if family_shared(fam, occurrences)]
	if occurrences:
		shared.sort(key=lambda fam: (-sum(occurrences.get(a, 0) for a in fam),
									 fam[0]))
	for fam in shared:
		if fam[0] in rna_accessions:
			rna_n += 1
			label = "R{}".format(rna_n)
		elif query_accessions.intersection(fam):
			query_n += 1
			label = "Q{}".format(query_n)
		else:
			prot_n += 1
			label = str(prot_n)
		for acc in fam:
			number[acc] = label
	return number


def _redact_command_line(argv):
	out, skip = [], False
	for i, token in enumerate(argv):
		if skip:
			out.append("<given>")
			skip = False
			continue
		flag = token.split("=", 1)[0]
		if flag in SECRET_FLAGS:
			out.append(flag + "=<given>" if "=" in token else token)
			skip = "=" not in token
		else:
			out.append(token)
	return " ".join(out)


def note_skipped(args, tool, reason):
	path = os.path.join(args.output, os.path.basename(
		os.path.normpath(args.output)) + "_runinfo.txt")
	try:
		with open(path, "a") as out:
			out.write("\nskipped {}\n  {}\n".format(tool, reason))
	except OSError:
		pass


def write_run_info(args, parser):
	os.makedirs(args.output, exist_ok=True)
	prefix = os.path.basename(os.path.normpath(args.output))
	out_path = lambda suffix: os.path.join(args.output, prefix + suffix)

	defaults = {a.dest: a.default for a in parser._actions}
	given, default = [], []
	for dest in sorted(vars(args)):
		value = getattr(args, dest)
		shown = "<given>" if dest in SECRET_ARGS and value else value
		(given if value != defaults.get(dest) else default).append(
			"  {:18} {}".format(dest, shown))

	with open(out_path("_runinfo.txt"), "w") as out:
		out.write("FlaGs3 {}\n".format(VERSION))
		out.write("run started   {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S %Z")))
		out.write("host          {}\n".format(platform.node()))
		out.write("python        {}\n".format(sys.version.split()[0]))
		out.write("working dir   {}\n".format(os.getcwd()))
		out.write("\ncommand line\n  {}\n".format(_redact_command_line(sys.argv)))
		out.write("\noptions set explicitly\n")
		out.write("\n".join(given) + "\n" if given else "  (none)\n")
		out.write("\noptions left at default\n")
		out.write("\n".join(default) + "\n" if default else "  (none)\n")

	sources = [(p, "_input.txt" if i == 0 else "_input{}.txt".format(i + 1))
			   for i, p in enumerate(list(args.input_list or []))]
	if args.blast_input:
		sources.append((args.blast_input, "_blast_input.txt"))
	args.input_copies = [suffix for _, suffix in sources]
	for source, suffix in sources:
		try:
			shutil.copyfile(source, out_path(suffix))
		except OSError as e:
			print("Warning: could not copy {} into the output directory "
				  "({}).".format(source, e))


def print_summary(args, prefix, extractor, families, rna_families, figures_written,
				  tree_written, want_tree, domain_table_written, features, sismis_mod,
				  blast_hits, genomad_written=False,
				  defence_written=False):
	print("\n{} -> {}".format(
		plural(len(extractor.sequences), "flanking protein"),
		plural(len(families) - len(rna_families), "family", "families")))
	print("\noutputs in {}/".format(os.path.relpath(args.output)
									if args.output.startswith(os.getcwd() + os.sep)
									else args.output))
	for name in figures_written:
		print("  {}".format(name))
	if tree_written:
		print("  tree/ ({}_tree.nwk, alignments, commands)".format(prefix))
	elif want_tree:
		print("  (tree skipped: fewer than 3 query sequences)")
	if domain_table_written:
		print("  {}_domains.tsv".format(prefix))
	if features:
		print("  {}_features.tsv".format(prefix))
	if args.sismis and sismis_mod:
		print("  {}_secretion.tsv / {}_sismis_diagnostics.txt".format(prefix, prefix))
	if genomad_written:
		print("  {}_genomad.tsv / {}_genomad_diagnostics.txt".format(prefix, prefix))
	if defence_written:
		print("  {}_defence.tsv / {}_defence_diagnostics.txt".format(prefix, prefix))
	for suffix in ("_operon.tsv", "_clusters.tsv", "_outdesc.txt", "_speciesInfo.txt",
				   "_QueryStatus.txt", "_flankgene_Report.log", "_clusterhits.tsv",
				   "_accessionIssues.txt", "_tree.fasta", "_flankgene.fasta", "_all.fasta",
				   "_runinfo.txt"):
		print("  {}{}".format(prefix, suffix))
	if extractor.ranges:
		print("  {}_rangeReport.tsv".format(prefix))
	if extractor.window_genes:
		print("  {}_window.tsv / {}_window.fasta".format(prefix, prefix))
	if flags_log.transcript_path():
		print("  {}_console.log".format(prefix))
	for suffix in getattr(args, "input_copies", []):
		print("  {}{}".format(prefix, suffix))
	if blast_hits:
		print("  {}_blast_hits.tsv / {}_blast_accessions.txt".format(prefix, prefix))
		if args.blast_input:
			print("  {}_blast_input.txt".format(prefix))


def write_window_tables(extractor, out_path):
	with open(out_path("_window.tsv"), "w") as out:
		out.write("#query\tassembly\tcontig\taccession\tstart\tend\tstrand\t"
				  "offset\tin_neighbourhood\tproduct\n")
		for row in sorted(extractor.window_genes):
			info = extractor.ranges.get(row)
			assembly = row.rsplit("|", 1)[-1]
			inside = extractor.drawn_accessions.get(row, set())
			for g in sorted(extractor.window_genes[row], key=lambda g: g.start):
				out.write("\t".join(str(v) for v in (
					row.rsplit("|", 1)[0], assembly, g.contig, g.accession,
					g.start, g.end, g.strand, g.offset,
					"yes" if g.accession in inside else "no",
					g.product or "")) + "\n")
	with open(out_path("_window.fasta"), "w") as out:
		for row in sorted(extractor.window_genes):
			info = extractor.ranges.get(row)
			query = row.rsplit("|", 1)[0]
			assembly = row.rsplit("|", 1)[-1]
			inside = extractor.drawn_accessions.get(row, set())
			for g in sorted(extractor.window_genes[row], key=lambda g: g.start):
				seq = extractor.window_sequences.get(g.accession)
				if not seq:
					continue
				lo, hi = relative_span(g, info)
				out.write(">{} {}\n{}\n".format(
					g.accession,
					window_header(query, assembly, g, len(seq), lo, hi,
								  g.accession in inside),
					seq))


def relative_span(gene, info):
	if info is None:
		return gene.start, gene.end
	if info.q_strand == "-":
		return info.q_end - gene.end, info.q_end - gene.start
	return gene.start - info.q_start, gene.end - info.q_start


def window_header(query, assembly, gene, length, lo, hi, drawn):
	return ("query={} assembly={} contig={} len={} start={} end={} strand={} "
			"offset={} in_neighbourhood={} product={}").format(
				query, assembly, gene.contig, length, lo, hi, gene.strand,
				gene.offset, "yes" if drawn else "no",
				gene.product or "hypothetical protein")
