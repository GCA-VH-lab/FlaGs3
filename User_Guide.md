# FlaGs3 — User Guide

Predicting protein functional association by analysis of conservation of genomic context (Flanking Genes).

---

## Installing

With conda:

```bash
bash build.sh
conda activate FlaGs3
```

`build.sh` creates the environment, verifies it, then offers each optional tool
in turn: the Pfam-A database, the DefenseFinder HMM profiles, MMseqs2, geNomad
and its database, DefenseFinder, PadLoc, and the two licensed tools. It is safe
to re-run — it detects what is already installed and skips it.

To answer up front instead of being asked:

```bash
bash build.sh --all                       # everything that needs no licence
bash build.sh --none                      # environment only
bash build.sh --with genomad,padloc       # just these
```

SignalP and DeepTMHMM need a package you obtained yourself, so `--all` reports
them as needing a path rather than stopping to ask. Give them one with
`--signalp PACKAGE.tar.gz` or `--deeptmhmm PACKAGE`, or run their installers
later.

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
the rest of the run. That applies to `-lth` and `-lsp` too: if either tool is
unusable where you pointed it, that feature is skipped and the reason is recorded
in `_runinfo.txt`, so a run never silently falls back to the cloud.

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
| `-rm`, `--remap` | off | Look a protein up again through IPG when the assembly it was paired with in the input produced nothing. Costs extra requests and downloads, so it is opt-in. |
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

`-g` and `-r` are two ways to say how big a neighbourhood is, and `-r` wins when
both are given. `-g 4` takes four genes either side whatever the distance; `-r
50000` takes everything within 50 kb whatever the gene count, which on a typical
bacterial genome is around 85 genes.

That difference is not free. Clustering compares every flanking protein against
every other, so its cost grows with the *square* of how many genes you take:
about 80x going from `-g 4` to `-r 50000`. On a large input that is the
difference between an hour and a fortnight.

`-sr` is the way out. The scanning tools take a genomic interval and do not care
how many genes are in it, so the span they see can be set independently of the
neighbourhood:

```
-g 5 -sr 50000
```

clusters and draws 11 genes per row while still giving the tools 100 kb of
context around each query. Use `-r` when you want the whole island coloured and
the input is small enough; use `-g` with `-sr` when it is not.

| Option | Default | Description |
|---|---|---|
| `-g`, `--gene N` | `4` | Flanking genes to take each side of the query. |
| `-r`, `--range BP` | off | Take every gene within this many bases of the query gene instead of a fixed count. Measured outwards from the query gene's own start and end; a gene straddling the edge is included, so the distance reached usually overshoots by part of one gene. Overrides `-g`. |
| `-sr`, `--scan_range BP` | neighbourhood span | Genomic span around the query handed to the scanning tools, independent of how many genes are clustered and drawn. |
| `-sm`, `--scan_margin BP` | `10000` | Extra sequence given to the scanning tools beyond the analysis range, so a system straddling the edge is called whole rather than cut in half. Hits reaching into the margin are kept and marked `partial`. `0` scans exactly the analysis range. |
| `-df`, `--defensefinder` | off | Call anti-phage defence systems with DefenseFinder, drawn as bands labelled with the system name. Needs `defensefinder_installer.sh`. |
| `-pl`, `--padloc` | off | The same with PadLoc. Both can run together. Needs `padloc_installer.sh`. |
| `-gn`, `--genomad` | off | Find proviruses and plasmids with geNomad, drawn as bands like Sismis'. Needs `genomad_installer.sh` then `genomad_loader.sh`. |
| `-gdb`, `--genomad_db DIR` | from tools table | geNomad database directory. |
| `-e`, `--ethreshold X` | `1e-3` | Inclusion E-value for clustering. Lower is stricter, giving more and smaller families. |
| `-n`, `--number N` | `3` | Jackhmmer iterations. More iterations find remoter homology but blur family boundaries. |
| `-cc`, `--cluster_collapse [ID[,COV]]` | off | Collapse near-identical flanking proteins with MMseqs2, cluster only the representatives, and give each member its representative's family. Bare `-cc` means `0.9,0.8`. Needs `mmseqs`; run `mmseqs_installer.sh`. |
| `-cr`, `--cluster_rna` | off | Also cluster flanking RNA genes into families. Uses nhmmer on RNA sequences where available, otherwise groups by product name. |

### Output and progress

| Option | Default | Description |
|---|---|---|
| `-o`, `--output DIR` | `output` | Result directory. A `_YYYYMMDD_HHMMSS` stamp of the run start is appended so repeated runs do not overwrite each other, and the stamped name is also the prefix on every file inside, so `-o myrun` produces `myrun_20260810_093134/myrun_20260810_093134_neighbors.svg`. |
| `-nt`, `--no_timestamp` | off | Use `-o` verbatim, without the stamp. Repeated runs then overwrite each other; use it when a pipeline needs a fixed path. |
| `-vb`, `--verbose` | off | Per-stage progress and a timing breakdown. Worth using on any long run. |
| — | always on | Every line the run prints is copied to `<prefix>_console.log` in the output directory. Nothing turns this off; it costs one open file. |
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
| `-lth`, `--local_tmhmm` | off | Predict transmembrane regions with a local DeepTMHMM rather than the BioLib cloud. Implies `--tmhmm`, so it is used on its own. |
| `-lsp`, `--local_signalp` | off | Predict signal peptides with a local SignalP rather than the BioLib cloud. Implies `--signalp`, so it is used on its own. |
| `--tools TSV` | `tools_table.tsv` | Table of external tool commands. |
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

### Reading the defence figure

Defence systems are drawn as bands behind the genes they span, labelled `D1`,
`D2` and so on rather than by name -- a row often carries several and the names
do not fit. The legend gives the full name for each number, and is split by the
tool that made the call, so a system under both DefenseFinder and PadLoc is one
both agreed on. It keeps the same number and colour in both panels.

Bands that overlap are stacked in separate lanes within the row. Three thin bars
means three systems sharing that stretch, not one wide one.

### Which flag controls which tool

Tools split by what they read, and that decides which flag sizes them:

| reads | tools | sized by |
|---|---|---|
| genes and proteins | domain scan, DeepTMHMM, SignalP, DefenseFinder, PadLoc | `-g` or `-r` |
| a stretch of DNA | Sismis, geNomad | `-sr`, whole genome without it |

So `-g 5 -sr 50000` clusters and draws 11 genes while Sismis and geNomad each see
100 kb around the query. Without `-sr` those two scan whole genomes, as they did
before.

If one tool needs a different span, put a number in the `scan_range` column of
`tools_table.tsv` for its row, or `genome` to give that one tool whole genomes.
That column overrides `-sr` for that tool only.

### When to collapse before clustering

Clustering compares every flanking protein against every other, so its cost grows
with the square of how many there are. `--cluster_collapse` reduces them to
MMseqs2 representatives first, clusters those, and hands each member its
representative's family.

Two honest caveats before you reach for it.

**The gain is usually smaller than the curve suggests.** FlaGs3 keys its sequence
table by accession, so proteins identical across genomes already collapse for
free — MMseqs2 only recovers what sits between your cutoff and 100% identity. On
real flanking proteins from unrelated genera that was 1.00x at 90% and 1.17x at
50%. It pays off when the input holds many strains of the same species; it does
very little when the genomes are all different species. Choosing `-g` over `-r`
is the far bigger lever.

**Collapsing hard can split families, not just merge them.** A sequence that
would have bridged two groups never gets searched once it is folded into a
representative, so the bridge disappears. Measured on a 950-protein run:

| setting | families | identical to an uncollapsed run |
|---|---|---|
| no collapse | 528 | — |
| `-cc` (0.9, 0.8) | 528 | every one |
| `-cc 0.5,0.8` | 533 | 521, with 71 proteins moved |

At the default the result was indistinguishable from not collapsing at all. At
50% it was not. Stay at `0.9,0.8` unless a run is otherwise impossible, and check
`_collapse.tsv` when a family looks wrong.

It is never switched on automatically, whatever the run size. If `mmseqs` is
missing the run stops rather than quietly clustering everything, since that turns
an hour into days with nothing in the output to say why.

### When a paired assembly holds nothing

Giving a protein with an assembly means FlaGs3 uses that assembly and never asks
NCBI where the protein lives. That is what makes a paired input fast. It also
means that if the assembly has since been withdrawn, fails to download, or simply
does not contain that accession, the protein is dropped with no second attempt —
and a verbose run says so:

```
>> 214 proteins whose paired assembly gave nothing; --remap looks them up again
   through IPG
```

`--remap` does exactly that, once, after extraction. Both failure modes look the
same at that point, so one pass covers a bad assembly and a good assembly missing
the protein alike. It is off by default because it costs extra NCBI requests and
possibly extra downloads; turn it on when the assembly column was compiled a
while ago.

Resolution through IPG goes out in chunks of 200 accessions. A single request for
thousands comes back slowly, truncated, or not at all, and a truncated report is
not an error — it looks exactly like those proteins having no assembly. Chunked,
a failure costs only its own 200 and the count is reported.

### Two tool tables

`tools_table.tsv` holds the shipped defaults and is what the repository tracks.
The installers never touch it. They write `tools_table.local.tsv` instead, which
is git-ignored, and FlaGs3 reads the default table first and lets the local one
override it row by row.

That split exists so a path like
`/home/you/miniconda3/envs/flags3-genomad/bin/genomad` stays on your machine
instead of turning up in every diff. A row in the local table only needs the
columns it changes; anything left blank falls back to the default.

To see what is actually in effect, read both files, or run with `-dbg`, which
records each external command as it runs.

## External tools

Every program FlaGs3 shells out to is defined in `tools_table.tsv` beside the code:
mafft, trimal, VeryFastTree, IQ-TREE, blastp, sismis, DeepTMHMM and SignalP. Edit a
row to change how one is invoked, or point `--tools` at your own copy. Rows you
leave out keep their defaults.

```
#name	command	directory
deeptmhmm	/opt/dtm-venv/bin/python3 predict.py --fasta {fasta} --output-dir {out}	/opt/DeepTMHMM
signalp	signalp6 --fastafile {fasta} --output_dir {out} --organism other --format txt --mode fast	/opt/signalp/bin
```

Two of these have installers, because they are licensed downloads with awkward
dependencies. Request the packages first — SignalP 6 from DTU, DeepTMHMM by
emailing `licensing@biolib.com` — then:

```bash
bash signalp_installer.sh   /path/to/signalp-6-package.tar.gz
bash deeptmhmm_installer.sh /path/to/deeptmhmm-package.tar.gz
```

Each builds a conda environment (`flags3-signalp`, `flags3-deeptmhmm`) with the
Python and PyTorch version that tool needs, and fills in its `tools_table.tsv` row
with the resulting paths, so `--local_signalp` and `--local_tmhmm` work without
further configuration.

The DeepTMHMM installer finishes by running `predict.py` on the bundled sample,
with its output shown rather than hidden. That both proves the install works and
pulls any model weights it needs, so a later `--local_tmhmm` run does not stop to
download anything. Set `DEEPTMHMM_TIMEOUT` (seconds, default 3600) if that run
needs longer.

SignalP 6 needs PyTorch below 2.0, and torch 1.x is built against NumPy 1.x, so
the installer pins `numpy<2` as well and re-checks after installing the package —
the package itself can pull NumPy 2 back in, which makes `signalp6` fail at import
with a NumPy 1.x/2.x mismatch.

DeepTMHMM's own `requirements.txt` pins `torch==1.5.0+cu92`, a CUDA 9.2 build from
2018 that exists only on PyTorch's index. The installer uses the CPU build instead:
cu92 predates current GPUs, and on one it fails inside cuBLAS as soon as it runs.
Set `DEEPTMHMM_GPU=1` to use the CUDA build anyway. Either way the torch pin is
stripped from the requirements before the rest are installed, so the two cannot
conflict.

`directory` is the working directory for that tool, and a relative program name is
looked up inside it. That is what makes tools with their own Python environments
workable: DeepTMHMM needs Python 3.8 and SignalP 6 needs PyTorch below 2.0, so
neither can share the FlaGs3 environment — naming their interpreter in `command`
avoids the clash.

Placeholders substituted per tool: `{in}`, `{out}`, `{fasta}`, `{threads}`,
`{mode}`, `{model}`, `{prefix}`, `{db}`, `{evalue}`, `{hits}`.

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

A figure taller than `-fh`/`--figure_height` (default 16383 px, the canvas limit
of Illustrator and librsvg) is written as `<name>_part1.svg`, `<name>_part2.svg`
and so on, splitting the rows in order. Raise the limit if you only ever open
figures in a browser, which has no such ceiling. Figures with a tree panel are
never split, because the tree spans every row; those are left whole with a
warning.

## Output files

`_runinfo.txt` and `_input.txt` are written before any work starts, so a run that
fails partway still records what it was asked to do. Your `--api_key` value is
masked in `_runinfo.txt`.

`_console.log` is a transcript of the run. Everything printed to the terminal
goes into it, and so does the output of the external tools FlaGs3 captures rather
than shows — mafft, trimal, blastp, sismis, the local feature tools and the
figure step. Lines that went to stderr, which is where `--debug` writes, carry a
`[stderr] ` prefix, so `grep -v '^\[stderr\]'` gives back the plain terminal
transcript and `grep '^\[stderr\]'` gives just the diagnostics. It is written
as the run goes, and survives a crash or Ctrl-C.

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
| `results_defence.tsv` | each defence system, the genes it spans, which tool called it, and which neighbourhoods it overlaps (`--defensefinder`, `--padloc`) |
| `results_defence_diagnostics.txt` | per-tool status and how many neighbourhoods were scanned |
| `results_genomad.tsv` | geNomad's proviruses and plasmids, each labelled with what it was called, and which neighbourhoods it overlaps (`--genomad`) |
| `results_genomad_diagnostics.txt` | per-genome geNomad status, windows scanned, and how much of the genome that came to (`--genomad`) |
| `results_sismis_diagnostics.txt` | per-genome Sismis status, how many windows were scanned, and how much of the genome that came to (`--sismis`) |
| `results_runinfo.txt` | how the run was invoked: version, host, command line, and every option split into those you set and those left at default |
| `results_collapse.tsv` | each MMseqs2 representative, its family, and its members, so a propagated family assignment can be traced back (`--cluster_collapse`) |
| `results_rangeReport.tsv` | per row: contig length, how much sequence was available up and downstream, how much the window actually reached, which sides were truncated, gene counts, and the span handed to the scanning tools |
| `results_console.log` | everything the run printed, plus the output of external tools that never reached the terminal |
| `results_input.txt` | a copy of the input list, so the results stay self-contained |
| `results_blast_hits.tsv` | BlastP hits used as queries: accession, E-value, bitscore, description |
| `results_blast_accessions.txt` | the same hits as a plain accession list — pass it to `-i` to repeat the run without searching again |
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

**geNomad says `Missing argument 'DATABASE'` or cannot read its database.**
The path in the third column of the `genomad` row of `tools_table.tsv` must be
the directory holding `genomad_marker_metadata.tsv`, not the directory holding
*that*. `genomad download-database DEST` creates `DEST/genomad_db`, so the path
you want is usually one level further in than the one you gave it. FlaGs3
descends that level for you and says so, but the table is clearer if it points
straight at the database.

**A tool reports its own dependencies as missing.**
Tools installed into a conda environment are run by absolute path, which does not
activate that environment. FlaGs3 puts the tool's own directory on its `PATH` so
sibling programs are found, but a dependency installed somewhere else entirely
still will not be. Check it is in the same environment: `conda list -n
flags3-genomad mmseqs2`.

**A tool produced nothing and the warning does not say why.**
Read `results_console.log`. External tools are run with their output captured,
so their own error messages never reach the terminal — the log has them in full,
under a `--- toolname (exit N) ---` header.

**A tool was skipped.**
Warnings name the missing dependency and the install command. The run continues
without that feature.