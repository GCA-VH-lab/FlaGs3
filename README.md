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
* One accession or one raw sequence, expanded to its homologues by BlastP —
  remotely through NCBI or locally with BLAST+, against RefSeq, GenBank or
  SwissProt
* NCBI RefSeq and GenBank, the MGnify Genomes catalogue, or your own local
  genome files
* Several input lists combined in one run

### Analysis

* Flanking-gene extraction with a configurable window
* All-vs-all clustering of neighbours into families with `pyhmmer` jackhmmer
* Optional clustering of RNA genes alongside proteins
* Phylogenetic trees from the query proteins: MAFFT alignment, gap-threshold
  trimming, then VeryFastTree or IQ-TREE with model selection and bootstrap
  support
* Pfam domain scanning, optionally annotated with InterPro entries, names,
  types and characterisation summaries
* Transmembrane helices via DeepTMHMM and signal peptides via SignalP
* Secretion-system detection via Sismis, matched to neighbourhoods by
  coordinate overlap

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

## Version history

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
