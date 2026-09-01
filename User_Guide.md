# FlaGs3 — User Guide

Predicting protein functional association by analysis of conservation of genomic context (Flanking Genes).

---

## Installing

With conda:

```bash
bash build.sh
conda activate FlaGs3
```

`build.sh` creates the environment, verifies it, and offers the two large
optional extras (sismis and the Pfam-A database). It is safe to re-run — it
detects what is already installed and skips it.

Or install the core dependencies by hand:

```bash
pip install biopython "pyhmmer>=0.12,<0.13" requests
```

Optional features need extra pieces, each only required if you use the matching
flag:

| Flag | Needs |
|---|---|
| `-t`, `--tree`, `--tree_order` | `mafft` and `VeryFastTree` on `PATH` |
| `--blast_mode local` | `blastp` from NCBI BLAST+ on `PATH`, and a local protein database |
| `-tm`, `--trimal_mode` | `gt` | trimal column filter: `gt`, `cons`, `st`, or a preset (`gappyout`, `strict`, `strictplus`, `automated1`, `nogaps`, `noallgaps`). |
| `-tv`, `--trimal_value` | `0.1` | Value for the modes that take one. |
| `-tx`, `--trimal_extra` | — | Extra trimal arguments passed through verbatim. |
| `-iq`, `--iqtree` | `mafft` and `iqtree` on `PATH` |
| `-d`, `--domains` | an HMM database — run `pfamA_loader.sh` to fetch Pfam-A |
| `-th`, `--tmhmm`, `--signalp` | `pip install pybiolib` and a network connection |
| `-ss`, `--sismis` | `pip install sismis` |

If a tool is missing, FlaGs3 prints a warning, skips that feature, and finishes
the rest of the run.

MGnify accessions are resolved through MGnify's API v2. MGnify's own API v1 was
switched off in September 2026, so older FlaGs3 releases will not resolve `MGYG`
accessions at all.

### One run at a time per `-tmp` directory

FlaGs3 takes a lock on its temporary directory, so a second run started against
the same `-tmp` stops with a message naming the holder instead of starting. Two
runs sharing that directory overwrite each other's downloads and double the
request rate against NCBI and EBI, which gets your host throttled — and the
throttling outlives the run that caused it, so downloads stay slow for a while
afterwards. To run two analyses at once, give each its own `-tmp`.

If a run is killed outright, its lock is left behind; the next run notices the
owning process is gone and reclaims it, so a stale lock never needs clearing by
hand. `--no_lock` skips the check entirely.

`mafft`, `VeryFastTree` and `iqtree` are all in `environment.yml`, so the conda
route covers the tree flags without anything extra. The table above matters only
if you install by hand.

Keep `pyhmmer` in the 0.12 series. FlaGs3 is developed against it, and sismis
(via gecco) requires it — an unpinned install can leave the two in conflict.

---

## The input list

One query per line. Two accepted forms, which can be mixed in the same file:

```
WP_047256880.1                          # protein only — FlaGs3 finds the genome
WP_047256880.1    GCF_000001765.3       # protein + genome, tab-separated
MGYG000454827_00001   MGYG000454827     # MGnify genome
```

**Protein only.** FlaGs3 resolves the genome through NCBI. This works for RefSeq
and GenBank proteins; `-m` controls how many genomes a protein may expand to when
it appears in several.

**Protein + genome.** Skips the lookup. The genome may be:

- an NCBI assembly — `GCF_...` (RefSeq) or `GCA_...` (GenBank)
- an MGnify Genomes accession — `MGYG...`

For MGnify genomes the protein accession must be the **exact locus tag** from
that genome's annotation (`MGYG000454827_00001`, not an `MGYP...` protein ID).
MGnify has no equivalent of NCBI's lookup service, so a bare MGnify protein
accession cannot be resolved to its genome — always supply the pair.

To find real locus tags for a genome:

```bash
curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000454827/downloads/MGYG000454827.faa" \
  | grep '^>' | head | sed 's/^>//' | cut -d' ' -f1
```

---

## Running it

```bash
python3 FlaGs3.py -i input.txt -u you@example.com -o results
```

`-u` is required by NCBI on any Entrez request. `-o` names the output directory
*and* becomes the prefix on every file inside it.

That is the whole default pipeline: it resolves each query to a genome, pulls the
flanking genes, clusters them, and writes the figure and tables. Everything else
is optional. See [All options](#all-options) for the full list and
[Examples](#examples) for common combinations.

## All options

### Required

| Option | Description |
|---|---|
| `-i`, `--input_list FILE` | Query list, one per line. See [The input list](#the-input-list). |
| `-u`, `--user_email ADDR` | Your email. NCBI requires it on every Entrez request. Not used for anything else. |

### Where genomes come from

| Option | Default | Description |
|---|---|---|
| `-ul`, `--use_local DIR` | — | Search a directory of local genomes before going to NCBI. A genome is a `.gff` and `.faa` sharing a basename; `.fna` and RNA FASTAs are picked up if present. Files may be gzipped. Anything not found falls back to NCBI. |
| `-m`, `--max_assemblies N` | `1` | How many genomes one protein may expand to when it occurs in several. Each genome becomes its own row. Raise to compare strains. Above `1`, row labels in the figures become `protein\|genome` so the rows stay distinguishable. |
| `-nc`, `--no_cross_db` | off | Keep protein and genome in the same database: RefSeq proteins (`WP_`, `NP_`, `YP_`, ...) resolve only to `GCF_` assemblies, INSDC proteins only to `GCA_`. A protein whose only assemblies sit in the other database is then reported as unresolved rather than annotated against a mirrored genome. Assemblies you supply yourself in the input file are never filtered. |
| `-api`, `--api_key KEY` | — | NCBI API key. Also raises the download rate cap from 5/s to 10/s. |
| `-tmp`, `--temporary DIR` | `./genomes` | Where downloads are stored. Deleted at the end unless `-k`. |
| `-k`, `--keep` | off | Keep downloaded genomes instead of deleting them. Useful for reruns — the directory can be fed straight back in via `--use_local`. |

### Finding queries with BlastP

Used only with `--blast_input`. See [Starting from one protein](#starting-from-one-protein).

| Option | Default | Description |
|---|---|---|
| `-bi`, `--blast_input FILE` | — | File holding one RefSeq accession, or one protein sequence as FASTA or bare residues. BlastP finds its homologues and they become the queries. |
| `-bh`, `--blast_hits N` | `50` | **Cap on how many BlastP hits are carried forward as queries.** Allowed 2-200. Hits are taken best-first, so a smaller number keeps the closest homologues. This is the main control on how big the run gets: every hit becomes a genome to download and a row in the figure. |
| `-be`, `--blast_evalue E` | `1e-5` | E-value cutoff for the BlastP search. Loosen it (e.g. `1e-3`) if a short or divergent query returns too few hits; tighten it to drop marginal ones. |
| `-bd`, `--blast_db NAME` | `refseq_select` | `refseq_select` (representative RefSeq proteins, faster), `refseq_protein` (full RefSeq), `genbank` (nr), or `swissprot`. Any other value is passed through unchanged, which is how a local database name or path is given. |
| `-bm`, `--blast_mode` | `remote` | `remote` uses NCBI QBLAST — nothing to install, but a search takes minutes. `local` runs `blastp` from NCBI BLAST+ against a local database: far faster, but you need the binary and the database. |

### What gets analysed

| Option | Default | Description |
|---|---|---|
| `-g`, `--gene N` | `4` | Flanking genes to take each side of the query. |
| `-e`, `--ethreshold X` | `1e-3` | Inclusion E-value for clustering. Lower is stricter, giving more and smaller families. |
| `-n`, `--number N` | `3` | Jackhmmer iterations. More iterations find remoter homology but blur family boundaries. |
| `-cr`, `--cluster_rna` | off | Also cluster flanking RNA genes into families. Uses nhmmer on RNA sequences where available, otherwise groups by product name. |

### Output and progress

| Option | Default | Description |
|---|---|---|
| `-o`, `--output DIR` | `output` | Result directory. A `_YYYYMMDD_HHMMSS` stamp of the run start is appended so repeated runs do not overwrite each other, and the stamped name is also the prefix on every file inside, so `-o myrun` produces `myrun_20260810_093134/myrun_20260810_093134_neighbors.svg`. |
| `-nt`, `--no_timestamp` | off | Use `-o` verbatim, without the stamp. Repeated runs then overwrite each other; use it when a pipeline needs a fixed path. |
| `-vb`, `--verbose` | off | Per-stage progress and a timing breakdown. Worth using on any long run. |
| `-dbg`, `--debug` | off | Diagnostics to stderr, each line stamped with seconds since start: per-file download timings split into limiter wait, time-to-first-byte and body transfer, plus HTTP status codes, external command lines with exit codes, and full tracebacks. Implies `--verbose`. Start here when downloads are slow or a tool silently produces nothing. |
| `-nl`, `--no_lock` | off | Skip the lock that stops two runs sharing one `-tmp` directory. Only safe if each run has its own `-tmp`. |
| `-v`, `--version` | — | Print the version and exit. |
| `-h`, `--help` | — | Print all options and exit. |

### Optional figures and annotation

Each needs an extra dependency. If it is missing, FlaGs3 warns, skips that
feature, and completes the rest of the run.

| Option | Needs | Description |
|---|---|---|
| `-t`, `--tree` | mafft, VeryFastTree | Also build a phylogenetic tree with the neighbourhoods aligned to its leaves (`_tree.svg`, `_tree.nwk`, `_tree.aln`). Does not change the main figure. |
| `-tm`, `--trimal_mode` | `gt` | trimal column filter: `gt`, `cons`, `st`, or a preset (`gappyout`, `strict`, `strictplus`, `automated1`, `nogaps`, `noallgaps`). |
| `-tv`, `--trimal_value` | `0.1` | Value for the modes that take one. |
| `-tx`, `--trimal_extra` | — | Extra trimal arguments passed through verbatim. |
| `-iq`, `--iqtree` | mafft, iqtree | Build the tree with IQ-TREE instead of VeryFastTree: ModelFinder picks the substitution model and 1000 ultrafast bootstrap replicates give branch support. Implies `--tree`. Far slower, so use it for the final figure rather than while exploring. |
| `--tree_order` | mafft, VeryFastTree | Order rows by tree leaf order, in the main figure and in `_operon.tsv`. Implies `--tree`. |
| `-d`, `--domains` | `-db` | Scan flanking proteins for domains and write `_domains.tsv`. Figures using them are controlled by the figure table. |
| `-db`, `--hmmdb [NAME=]PATH` | `./pfam_db/Pfam-A.hmm` | HMM database for `--domains`: a `.hmm` file, or a directory of `.hmm` files such as DefenseFinder's `profiles/`. Repeat for several. `NAME=` labels it in the outputs; otherwise the file or directory name is used. Models carrying a gathering threshold are scored by it, the rest by `-e`. |
| `-hc`, `--hmm_coverage [NAME=]Q[,H]` | — | Minimum fraction of the protein (Q) and of the model (H) an alignment must span. `NAME=` applies it to one database, omitting it applies to all. Use for full-length protein models such as DefenseFinder (`0.7,0.5`); leave off for Pfam, where partial coverage is normal. |
| `-ip`, `--interpro FILE` | `interpro_metadata_processed.tsv` | InterPro metadata table (`.tsv` or `.tsv.gz`). Adds the InterPro entry, name, type and short characterisation/informativeness summaries to `_domains.tsv`, joined on the Pfam accession. Needs `accession` and `pfam_members` columns. Found automatically in the working directory or next to `FlaGs3.py`; if it isn't there the domain table is simply written without those columns. A path you pass yourself must exist. |
| `--clans FILE` | — | `Pfam-A.clans.tsv.gz`. Colours domains by clan rather than family, which groups related domains together. |
| `-th`, `--tmhmm` | pybiolib, network | Predict transmembrane regions with DeepTMHMM, drawn as double red dotted lines on the domain figure. Uploads your sequences to the BioLib cloud. |
| `-sp`, `--signalp` | pybiolib, network | Predict signal peptides with SignalP-6, drawn as black triangles on the domain figure. Also uploads sequences. |
| `-ss`, `--sismis` | sismis | Scan each genome for secretion systems and write `_secretion.tsv` plus `_secretion.svg`, noting which neighbourhoods each hit overlaps. Downloads the genomic FASTA per genome. |

`--tmhmm`, `--signalp` and `--sismis` run concurrently with each other, since
each spends a lot of time waiting on a remote service or a subprocess.

### Performance

| Option | Default | Description |
|---|---|---|
| `-c`, `--cpu N` | auto | Worker cap for clustering, domain scanning and downloads. |

Download rate is capped independently of `-c` at 5 requests/second, or 10 with
`-api`. This is deliberate: raising the worker count on a fast machine would
otherwise raise the request rate and get the run throttled or rejected.

## Examples

Default run — neighbours figure and data tables:

```bash
python3 FlaGs3.py -i input.txt -u you@example.com -o myrun
```

Wider neighbourhood across several strains, ordered by a tree:

```bash
python3 FlaGs3.py -i input.txt -u you@example.com -o myrun \
  -g 6 -m 5 --tree --tree_order -vb
```

Domain annotation with clan colouring:

```bash
python3 FlaGs3.py -i input.txt -u you@example.com -o myrun \
  --domains --hmmdb pfam_db/Pfam-A.hmm --clans pfam_db/Pfam-A.clans.tsv.gz
```

Everything on, with an API key for faster downloads:

```bash
python3 FlaGs3.py -i input.txt -u you@example.com -o myrun -api YOUR_KEY \
  --tree --tree_order --domains --hmmdb pfam_db/Pfam-A.hmm \
  --tmhmm --signalp --sismis --cluster_rna -vb
```

Reuse genomes from a previous run instead of downloading again:

```bash
python3 FlaGs3.py -i input.txt -u you@example.com -o run1 -k
python3 FlaGs3.py -i input.txt -u you@example.com -o run2 --use_local ./genomes
```

---

## Starting from one protein

Instead of writing out a list of homologues yourself, you can hand FlaGs3 a single
starting point and let BlastP find them, the same way webFlaGs does. Put **one** of
these in a file and pass it with `-bi`/`--blast_input`:

- a RefSeq protein accession (`WP_`, `NP_`, `YP_`, `XP_`, `AP_`), or
- a protein sequence, as FASTA or as bare residues over any number of lines.

```bash
echo "WP_047256880.1" > start.txt
python3 FlaGs3.py --blast_input start.txt -u you@example.com -o run1
```

You can also mark entries in the main input list: a line whose second column is
`BLAST` is expanded instead of used directly, so a single file can mix ordinary
accessions with things to expand.

```
WP_000028540.1	GCF_022493555.1
WP_061892803.1
WP_201476908.1	BLAST
MKKATLARQLVDGT	BLAST
```

Every marked entry gets its own search and all the hits are pooled. `--blast_input`
works exactly as before and can be combined with marked lines.

The hits become the queries. They are appended to whatever `-i` holds, so you can
mix a curated list with a BlastP expansion in one run; anything already in `-i` is
not added twice, matching on accession regardless of version, and the assembly you
paired it with is kept. `-i` is optional when `--blast_input` is given.

**Controlling how many hits you get.** `--blast_hits` caps how many hits become
queries (default 50, allowed 2-200) and `--blast_evalue` sets the cutoff (default
`1e-5`). Hits are taken best-first, so `--blast_hits 10` gives you the ten closest
homologues. The cap matters: each hit becomes a genome to download, a neighbourhood
to extract and a row in the figure, so 200 hits is a much longer run than 20.

```bash
python3 FlaGs3.py --blast_input start.txt --blast_hits 20 --blast_evalue 1e-10 \
    -u you@example.com -o run1
```

The hits that were used are written to `_blast_hits.tsv` with their E-values and
descriptions, and your starting file is copied to `_blast_input.txt`.

**Remote versus local.** The default `--blast_mode remote` needs nothing installed
but a QBLAST search takes minutes, and much longer when NCBI is busy. FlaGs3 prints
the job id and NCBI's own time estimate as soon as the job is accepted, then a
progress line each minute while it waits, so a slow search is distinguishable from
a hung one. The printed link opens the job on NCBI's site. `-bw`/`--blast_wait`
caps the wait in minutes (default 60); giving up does not cancel the job, and the
link stays valid. `--blast_mode local` runs `blastp` from NCBI
BLAST+ against a local database and is far faster, but you need the binary and the
database — `--blast_db` then takes the database name or path, and `-c/--cpu` sets
`-num_threads`. An accession is resolved to its sequence through NCBI first, so
both modes accept the same input.

## Redrawing the figures

Rendering is the cheapest stage and the one most subject to taste, so it can be
re-run on its own against a finished run — no downloads, no clustering:

```bash
python3 flags_redraw.py --data testout_20260819_125350
```

Every run calls `flags_redraw.py`, which copies the `visualisation_table.tsv` from
the FlaGs3 directory into the output directory if there isn't one there already.
That root copy is the standard set — edit it and every later run starts from your
version, not a built-in one. That table *is* the figure list, and it carries
every drawing parameter — nothing about a figure's look is buried in the code.
Edit it and re-apply:

```
#name mode  tree_width  features_allowed  family_numbers  font_size row_height  gene_height gene_gap  bases_per_pixel pad domain_height label_step  arrow_head  min_gene_width  band_opacity
neighbors versatile False cluster_rna TRUE  13  26  8 default 10.4  16  6 12  default default default
tree  triangles 1 cluster_rna TRUE  13  24  20  1 False 16  default default default default default
classic classic False cluster_rna TRUE  12  20  15  default 10.4  16  6 12  7 13  default
```

Each row produces `<prefix>_<name>.svg`. Any numeric cell may be `default`.

| Column | Meaning |
|---|---|
| `name` | Figure name, used as the filename suffix. |
| `mode` | `versatile` (the current look), `triangles` (fixed-width genes beside a tree), or `classic` (original FlaGs3: tighter rows, blunter arrows, numbers inside them, query protein solid black). |
| `tree_width` | `False` for no tree, otherwise a multiplier on the default panel width. Works in **every** mode, not just `triangles`. |
| `features_allowed` | Comma-separated: `cluster_rna`, `domains`, `tmhmm`, `signalp`, `sismis`, `monochrome`, `none`. `monochrome` turns off family colouring. Tokens compose: `sismis` only adds secretion bands, `domains` only adds domain wedges, and family colours and numbers stay unless you ask for `monochrome` or `family_numbers FALSE`. |
| `family_numbers` | `TRUE`/`FALSE` — show family number labels and family colouring. They still appear when domains or secretion bands are drawn; family numbers are prefixed `G` on figures that also draw domains, so they are not confused with domain numbers; only plain family numbers take it, since `Q` and `R` already say what the family is. Under domain wedges the gene is filled in a pastel tint of its family colour, opaque so it hides any secretion band beneath and stays clear of the wedges on top, with the outline in the same tint. RNA, pseudogene and query outlines keep their full-strength accent. `FALSE` leaves genes as bare outlines when wedges are drawn, so the domains alone carry the colour. Domain and secretion numbers are unaffected either way. |
| `font_size`, `row_height`, `gene_height`, `gene_gap`, `pad` | Type and layout sizes in px. |
| `bases_per_pixel` | Genomic scale. `False` gives fixed-width genes. |
| `domain_height`, `label_step`, `arrow_head`, `min_gene_width`, `band_opacity` | Domain wedge height, label spacing, classic arrow head length, smallest arrow width, secretion band opacity. Band height follows `gene_height`; the type label sits to the right of the row in black at 85% of `font_size`. |

A gene that is more than one thing shows both: the fill gives the family, and the
outline gives the rest — green for RNA, navy for a pseudogene, black for a query
protein. A gene that is only its family is outlined in its own colour, and an
unclustered one in mid grey. RNA genes, pseudogenes and queries also get a
double-weight outline. So a clustered RNA gene keeps its family colour instead of
losing it to the RNA styling.

Every figure carries a legend for whatever it actually draws: domains, secretion
systems, transmembrane helices and signal peptides, the query protein, RNA genes
and pseudogenes, and a distance scale bar under any figure that shows a tree.
Nothing is listed that is not in the picture.

`--pdf` works the same way here, so a finished run can be turned into PDFs
without redoing the analysis.

`--format` takes a table from anywhere, `-o/--output` writes elsewhere, and
`--write_table` drops a starting table into a run directory. Because the main run
shells out to the same script, a redraw reproduces the run's figures exactly.
`-nf`/`--no_figures` skips drawing entirely; `-f`/`--figures` points a run at your own table.

## Output files

`_runinfo.txt` and `_input.txt` are written before any work starts, so a run that
fails partway still records what it was asked to do. Your `--api_key` value is
masked in `_runinfo.txt`.

Every file is prefixed with the output directory's name, stamp included — a run
with `-o results` writes `results_20260810_093134/results_20260810_093134_operon.tsv`.
The tables below drop the stamp and write `results_...` for readability.

**Figures**

| File | Contents |
|---|---|
| `results_neighbors.svg` | the main diagram: one row per query, genes as arrows coloured by family |
| `results_tree.svg` | same rows aligned to a phylogenetic tree (`--tree`) |
| `results_tree.nwk` | the tree in Newick format |
| `tree/` | the alignment, the trimmed alignment, the Newick tree and the exact alignment/trimming/tree commands used |
| `results_domains.svg` | neighbourhoods with domains, TM regions and signal peptides drawn on |
| `results_secretion.svg` | neighbourhoods with predicted secretion systems marked (`--sismis`) |

**Tables**

| File | Contents |
|---|---|
| `results_operon.tsv` | one row per flanking gene: query, genome, family, strand, offset, coordinates, length, contig, product |
| `results_clusters.tsv` | each family and its members |
| `results_outdesc.txt` | families as readable blocks: `family(occurrences)`, accession, product description |
| `results_features.tsv` | transmembrane and signal-peptide regions per protein (`--tmhmm`/`--signalp`), so the figures can be redrawn later |
| `visualisation_table.tsv` | the figure list for this run; edit and re-apply with `flags_redraw.py` |
| `results_domains.tsv` | one row per domain hit: protein, family, domain, Pfam accession, clan, coordinates, E-value (`--domains`); with `--interpro`, also the InterPro entry, name, type, characterisation status, informativeness and a one-line interpretation |
| `results_jackhits.tsv` | per-protein jackhmmer inclusion lists — the audit trail behind the families |
| `results_speciesInfo.txt` | genome and organism per query row |
| `results_QueryStatus.txt` | which genomes each query resolved to, and whether it produced a row |
| `results_accessionIssues.txt` | queries that produced nothing, and why |
| `results_flankgene_Report.log` | each neighbourhood as a compact family chain |
| `results_secretion.tsv` | Sismis hits and which neighbourhoods they overlap (`--sismis`) |
| `results_sismis_diagnostics.txt` | per-genome Sismis status: scanned, nothing found, or skipped and why (`--sismis`) |
| `results_runinfo.txt` | how the run was invoked: version, host, command line, and every option split into those you set and those left at default |
| `results_input.txt` | a copy of the input list, so the results stay self-contained |
| `results_blast_hits.tsv` | BlastP hits used as queries: accession, E-value, bitscore, description (`--blast_input`) |
| `results_blast_input.txt` | a copy of the BlastP starting accession or sequence (`--blast_input`) |

**Sequences**

| File | Contents |
|---|---|
| `results_tree.fasta` | one query protein per row, named by row id — the tree input |
| `results_flankgene.fasta` | the flanking proteins, headers `accession\|product` |
| `results_all.fasta` | both of the above in one file |

### Reading the figure

Genes are arrows pointing in their direction of transcription, normalised so the
query always points right. Arrows sharing a colour and number are one family.
The query itself is outlined in black at the centre of each row.

Grey means the gene had no family (a singleton). RNA genes keep a green outline,
pseudogenes navy.

### Reading `results_operon.tsv`

The `query` column is the query protein and `assembly` is the genome it was found
in. They are separate columns because one protein appearing in several genomes
produces one row per genome, so the pair — not the protein alone — identifies a
row. The `accession` column is the flanking gene itself. `offset` is the position
relative to the query: `0` is the query, negative upstream, positive downstream.
`contig` is the sequence the gene lies on, which is what Sismis hits are matched
against.

With `--tree_order`, rows appear in tree leaf order rather than input order.

---

## When something goes wrong

**"No flanking neighbourhoods could be extracted for any query."**
Nothing matched. Check `results_accessionIssues.txt` — it separates *no genome
resolved* from *genome found but the protein was not in it*. The second usually
means the accession does not appear in that genome's annotation, which for
MGnify genomes normally means the locus tag is wrong.

**Queries silently missing from the figure.**
`results_QueryStatus.txt` lists every query and whether it produced a row.

**Downloads look stuck.**
Run with `-vb`; downloads report as they complete. Requests are rate-limited
(5/s, or 10/s with `-api`) so a large list takes a while by design — this keeps
NCBI and EBI from rejecting the run.

**The diagram is unexpectedly wide.**
A neighbour lying very far from the query stretches the canvas. This normally
means a fragmented assembly where the query sits near a contig edge.

**A tool was skipped.**
Warnings name the missing dependency and the install command. The run continues
without that feature.