# FlaGs3

Predicting protein functional association by analysis of conservation of genomic
context (Flanking Genes).

See `User_Guide.md` for usage and `Architecture.md` for internals.

## What it does

Takes a set of proteins, finds the genes flanking each one in its genome, groups
those neighbours into families by sequence similarity, and draws the
neighbourhoods so conserved gene arrangements are visible at a glance.

### Getting the proteins in

* Protein accessions, alone or paired with a specific genome
* Accessions or raw sequences expanded to their homologues by BlastP, from a
  separate file or marked `BLAST` in the main list —
  remotely through NCBI or locally with BLAST+, against RefSeq, GenBank or
  SwissProt
* NCBI RefSeq and GenBank, the MGnify Genomes catalogue, or your own local
  genome files
* Several input lists combined in one run

### Analysis

* Flanking-gene extraction by gene count or by distance in bases
* Defence systems from DefenseFinder and PadLoc, proviruses and plasmids from
  geNomad, and secretion systems from Sismis, drawn as labelled bands
* All-vs-all clustering of neighbours into families with `pyhmmer` jackhmmer
* Optional clustering of RNA genes alongside proteins
* Phylogenetic trees from the query proteins: MAFFT alignment, gap-threshold
  trimming, then VeryFastTree or IQ-TREE with model selection and bootstrap
  support
* Domain scanning against any HMM set — Pfam, or directories of profiles such as
  DefenseFinder and PADLOC — optionally annotated with InterPro entries, names,
  types and characterisation summaries
* Transmembrane helices via DeepTMHMM and signal peptides via SignalP
* Secretion-system detection via Sismis, matched to neighbourhoods by
  coordinate overlap
* The sequence-scanning tools take their span from `-sr`, independently of how
  many genes are clustered, with `-sm` padding it so a system at a window edge is
  still called whole
* `--remap` looks a protein up again through IPG when the assembly it was paired
  with in the input turns out to hold nothing

### Figures

* Neighbourhood maps with genes coloured and numbered by family
* Trees with neighbourhoods aligned to their leaves, in any figure style
* Domain wedges, transmembrane and signal-peptide marks, secretion bands —
  drawn together or separately, each with its own legend
* Three styles: the current look, tree-aligned triangles, and the original
  FlaGs layout
* Every figure defined by a row in an editable table: which overlays, which
  colours and numbers, and every size and spacing value
* Figures redrawn from a finished run without repeating the analysis, so the
  look can be tuned without re-downloading anything
* Overlapping bands stacked in lanes, so two tools calling one locus are both
  visible, with defence systems numbered `D1`, `D2` and the legend split by which
  tool made the call
* Figures too tall for an editor written as `<name>_part1.svg` and so on, at a
  ceiling set by `-fh`
* PDF output alongside the SVGs, from any run or redraw

### Results

* Tables for neighbourhoods, families, family descriptions, clustering
  evidence, domains, features, secretion hits, species and per-query status
* Sequence files for the queries, the flanking proteins, and both together
* A record of every run: options used, versions, and a copy of the input
* Timestamped output directories, so runs never overwrite each other

### Running it

* Parallel downloads with rate limiting, and cloud tools run concurrently
* Local genome reuse, and a lock that stops two runs colliding
* Graceful degradation: a missing tool or an unreachable NCBI skips that
  feature rather than ending the run
* Verbose progress with per-stage timings, and a debug mode reporting HTTP
  status, external commands and full tracebacks
* A full console transcript saved beside the results, surviving crashes and
  Ctrl-C

## Version history

**2.2.0** — the scan window's genes and proteins are written out and its FASTA
headers say what they are, clustering no longer exhausts memory on a dense
subfamily, and every unreferenced function and import is gone

**2.1.0** — `--remap` looks a protein up again through IPG when its paired
assembly gives nothing, IPG resolution runs in chunks rather than one request,
and figures too tall for an editor are written as parts

**2.0.0** — defence system and mobile element calling, windowed scanning, and
the scale work needed to run all of it on thousands of genomes.

Two changes affect existing output. Sismis' passthrough columns in
`_secretion.tsv` are prefixed `sismis_`, because they collided with the
normalised `start`/`end`/`type` and the renderer was reading the wrong ones.
Family numbers, clusters and FASTA output are now ordered deterministically
rather than following whatever order the input arrived in.

* `-df` / `-pl` call defence systems with DefenseFinder and PadLoc, drawn as
  bands labelled `D1`, `D2` with the names in a legend split by which tool
  called them
* `-gn` finds proviruses and plasmids with geNomad, each labelled with what it
  was actually called rather than a generic label
* `-sm` pads a scan window so a system at its edge is called whole and marked
  partial rather than dropped
* Sismis and geNomad are invoked once for many genomes rather than once each;
  their cost is dominated by loading a database, so this is the difference
  between hours and weeks on a large input
* Parsed genomes are freed as soon as they are finished with, cutting retained
  memory 226x on a large input
* DefenseFinder and PadLoc see every gene within `-sr` of the query rather than
  only the drawn neighbourhood, so a system reaching past the figure is still
  called whole
* Sismis and geNomad invoked once for many genomes rather than once each, since
  their cost is dominated by loading a database
* Installers for MMseqs2, geNomad, DefenseFinder and PadLoc, all offered by
  `build.sh` in turn or selected with `--all`, `--none` and `--with LIST`; none
  needs a user-supplied package
* Machine-specific tool paths go to a git-ignored `tools_table.local.tsv`, and
  `tools_table.tsv` keeps the shipped defaults

**1.4.0** — optional MMseqs2 collapsing before clustering, and connected
components taken on a symmetrised graph so a protein can no longer land in two
families

**1.3.0** — neighbourhoods can be defined as a distance in bases rather than a
number of genes, with a per-row report of what each contig could actually
supply, and a separate span for the scanning tools

**1.2.0** — every line the run prints is also saved as `<prefix>_console.log`,
including the output of external tools that never reached the terminal

**1.1.0** — DeepTMHMM and SignalP can run locally instead of in the cloud, and
every external tool command is configurable in `tools_table.tsv`

**1.0.13** — BlastP hits saved as a reusable accession list, and a loader for the
DefenseFinder profiles

**1.0.12** — domain scanning takes any HMM set, including directories of profiles
such as DefenseFinder

**1.0.11** — requests identify the tool to NCBI as `flags3`

**1.0.10** — entries in the main input list can be marked `BLAST` for expansion

**1.0.9** — the tree figure draws again, with the same gene styling as the others

**1.0.8** — pseudogenes are recognised in every annotation style, and classic
mode fills unclustered genes white

**1.0.7** — remote BlastP shows the NCBI job id and progress while it waits,
not after

**1.0.5** — conserved genes are coloured as well as numbered

**1.0.4** — figures open correctly in Illustrator and other SVG 1.1 editors

**1.0.3** — a gene conserved across genomes is no longer treated as a singleton,
and families are numbered in order of how widespread they are

**1.0.2** — trimal for alignment trimming, tree files in a `tree/` subfolder,
and figure styling closer to the original FlaGs

**1.0.1** — PDF output alongside every figure

**1.0.0** — First FlaGs3 release
