# FlaGs3 User Guide

This guide is written to be followed top to bottom the first time. Every
command in it is meant to be pasted as is.

## 1. Installing

FlaGs3 needs Python 3.11 or later. The recommended way is pipx, which
installs the `flags3` command into its own environment and puts it on your
PATH once and for all:

```
git clone https://github.com/GCA-VH-lab/FlaGs3.git
pipx install ./FlaGs3
flags3 --version
```

pipx itself comes from your package manager (`pacman -S python-pipx`,
`apt install pipx`, `brew install pipx`; then `pipx ensurepath` once).
To upgrade later: `pipx upgrade flags3`, or `pipx install --force` from
the updated checkout.

Without pipx, a virtual environment does the same job, except that it
must be activated (`source .venv/bin/activate`) in every new terminal
before `flags3` is found:

```
cd FlaGs3
python -m venv .venv
source .venv/bin/activate
pip install .
```

Developers working on the code use `pip install -e ".[test]"` in
that environment so edits take effect without reinstalling.

`--pdf` needs the cairo library, which Linux has and macOS gets from
`brew install cairo`; without it SVGs are still written and `--pdf` prints
a warning.

If you live in conda, the third route is one environment with FlaGs3 and
the core tools together:

```
cd FlaGs3
conda env create -f environment.yml
conda activate flags3
flags3 --version
```

That brings mafft, trimal, VeryFastTree, IQ-TREE, BLAST+, MMseqs2 and
cairo with it, so `flags3 install core` and `flags3 install mmseqs` find
them on PATH and build nothing. The environment must be activated in every
terminal, as any conda environment.

Whichever route, everything else — external tools, databases, downloaded
genomes — lives under `~/.flags3/`.

### External tools: `flags3 install`

Nothing beyond Python is needed for the basic run: genome download,
neighbourhood extraction, clustering (jackhmmer runs inside Python) and the
figures. Everything else is installed on demand:

```
flags3 install                 # list components and whether they are installed
flags3 install core            # mafft, trimal, VeryFastTree, IQ-TREE, BLAST+   (for --tree)
flags3 install pfam            # Pfam-A HMMs and clans, ~1.5 GB               (for --domains)
flags3 install mmseqs          # MMseqs2                                       (for -cm mmseqs_cluster)
flags3 install sismis          # Sismis                                        (for --sismis)
flags3 install genomad         # geNomad and its 1.6 GB database              (for --genomad)
flags3 install defensefinder   # DefenseFinder and its models                  (for --defensefinder)
flags3 install padloc          # PadLoc and its database                       (for --padloc)
flags3 install defence-hmm     # DefenseFinder HMM profiles alone              (for --domains -db)
flags3 install deeptmhmm       # DeepTMHMM2, local and licence-free            (for --tmhmm -lth)
flags3 install --all           # all of the above that are not yet installed
flags3 install --update        # refresh the databases of what is installed
flags3 install --force pfam    # reinstall one from scratch
```

One component is licensed and needs a package you download yourself from
DTU, then hand to the installer:

```
flags3 install signalp ~/Downloads/signalp-6.0h.fast.tar.gz
```

And one is a file that is not downloadable at all — the InterPro metadata
table used to annotate domains. If you have it:

```
flags3 install interpro /path/to/interpro_metadata_processed.tsv
```

The installer first looks for each tool on your PATH and uses it if it is
there — from conda, Homebrew, a cluster module, wherever. What is missing
it builds in its own environment under `~/.flags3/envs/` with micromamba,
using the one on your PATH if you have it and otherwise downloading a
single binary into `~/.flags3/tools/`. It never uses or changes a conda
installation of yours. MMseqs2 comes as the static release binary. After
each install the tool's location is written into
`~/.flags3/tools_table.tsv`, which every run reads. To remove a component,
delete its directory under `~/.flags3/envs/`, `~/.flags3/tools/` or
`~/.flags3/db/`.

Without `-lth`/`-lsp`, `--tmhmm` and `--signalp` run on the BioLib cloud:
one job per tool (up to 2,000 proteins for DeepTMHMM and 1,000 for SignalP
per job), both submitted at once, with a link printed for each so you can
watch it on biolib.com. A run's wait is one queue wait, not one per batch.

## 2. The genome directory

Genomes are downloaded once into `~/.flags3/genomes/` and reused by every
later run. Nothing is ever deleted from it. Three flags change that:

```
-gd DIR       use this directory instead: looked in first, downloaded into
--offline     use the directory only; never contact NCBI or MGnify
-lf           look bare queries up in the directory before asking NCBI
--no_cache    leave the cache alone: download into <run>/genomes/ instead
```

`-gd` is how you point a run at genomes you already have — a colleague's
cache when reproducing their run, a folder on shared storage, or your own
assemblies that are in no database. Files there are recognised by name:
`X_genomic.gff` (or `.gff3`, gzipped or not) and `X_protein.faa` sharing
the basename `X`, plus optionally `X_genomic.fna` for the sequence
scanners and `X_rna_from_genomic.fna` for RNA clustering. A query paired
with `X` in the input list uses those files without any network access.
With `--offline`, bare queries are looked up in the directory's protein
files instead of at NCBI, and nothing is downloaded. With `-lf`
(`--local_first`) they are looked up there first and only the ones not
found go to IPG, so one list can mix your own assemblies with NCBI
queries in an online run.

`--no_cache` is for a run whose genomes you don't want kept: a download
test, or a one-off on hundreds of assemblies. They still land in
`<run>/genomes/`, so every stage can be rerun; deleting the run directory
deletes them.

## 3. The input list

A text file, one query per line:

```
WP_005328829.1
NP_416433.1	GCF_000005845.2
MGYG000454827_00001	MGYG000454827
# comments and blank lines are ignored
```

- A bare protein accession is resolved through NCBI's IPG to the assembly
  that carries it (RefSeq first). `-m 3` takes up to three assemblies,
  `-nc` keeps RefSeq proteins in RefSeq assemblies only.
- A tab and an assembly accession pin the genome. Add `-rm` if a paired
  query might be annotated under a different accession in that assembly;
  IPG is then asked, and the query is matched under whatever name it has
  there.
- MGnify genome accessions (`MGYG…`) are recognised and fetched from the
  MGnify catalogue; their proteins are named by locus tag.
- Several lists can be given: `-i a.txt -i b.txt`.

### Starting from one protein: BlastP

Instead of, or in addition to, a list, give one starting point and let
BlastP find its homologues, which then become the queries:

```
flags3 run -bi start.txt -u you@example.org -o myrun
```

`start.txt` holds either one accession or a protein sequence (FASTA or
bare residues). The same can be done inside the list: a line with a tab
and `BLAST` after an accession or a sequence is expanded in place.

```
WP_005328829.1	BLAST
MKVLLAKQRTV...	BLAST
```

```
-bm remote      NCBI QBLAST, nothing to install, minutes per search (default)
-bm local       blastp from flags3 install core; -bd is then a local database path
-bd refseq_select   or refseq_protein, genbank, swissprot
-be 1e-5        E-value cutoff
-bh 50          hits carried forward as queries; each is a genome and a row
-bw 60          minutes to wait on NCBI's queue
```

The hits go to `blast/hits.tsv` with E-values and descriptions, and
`blast/accessions.txt` is a ready input list for repeating the run without
searching again.

## 4. Running

```
flags3 run -i list.txt -u you@example.org -o myrun
```

`-u` is your email address, which NCBI requires. `-o` names the run
directory; a timestamp is appended unless you add `-nt`. The run prints one
line as each stage starts and finishes, more with `-vb`.

### What to analyse

```
-g 4            flanking genes on each side of the query (default 4)
-r 5000         instead of -g: every gene within 5000 bp of the query
-sr 20000       span each side of the query handed to the sequence scanners
                (Sismis, geNomad); without it they see the whole contig
-sm 10000       extra sequence beyond -sr so a system at the edge is called whole
```

### Clustering

```
-cm jackhmmer   the method: a row of the tools table (default jackhmmer;
                mmseqs_cluster or mmseqs_cluster_exhaustive after
                flags3 install mmseqs)
-cr             cluster flanking RNA genes too (nhmmer; RNAs without a
                sequence are grouped by product name)
```

Thresholds, iterations and sensitivity live in the tools table, not on the
command line — see section 9.

### Tree

```
-t              mafft → trimal → VeryFastTree on the query proteins
-iq             IQ-TREE with ModelFinder and 1000 ultrafast bootstraps instead
-to             order rows in figures and tables by the tree
-tm gt -tv 0.1  trimal mode and value (also cons, st, gappyout, strict, …)
```

### Domains and protein features

```
-d              scan flanking proteins against Pfam-A (from flags3 install pfam)
-db NAME=PATH   another HMM source: a .hmm file, or a directory of profiles;
                repeat for several, e.g. -db defensefinder=~/.flags3/db/defensefinder/profiles
-hc NAME=Q,H    minimum fraction of the protein (Q) and of the model (H) a hit
                must cover; without NAME= it applies to every source
-e 1e-3         E-value for models without a gathering threshold
-cl FILE        Pfam-A.clans.tsv(.gz): colour domains by clan instead of family
-ip FILE        an InterPro metadata table, if you have one
-th             transmembrane regions with DeepTMHMM2, on the BioLib cloud or,
                with -lth, locally from flags3 install deeptmhmm (no queue, no
                licence, but slower on a CPU); --tmhmm_app DTU/DeepTMHMM picks
                the 1.0 model on the cloud
-sp             signal peptides (SignalP 6 on BioLib); -lsp runs the licensed
                local SignalP from flags3 install signalp
```

### Sequence and system scanners

```
-ss             secretion systems (Sismis)
-gn             proviruses and plasmids (geNomad); -gdb DIR for another database
-df             anti-phage defence systems (DefenseFinder)
-pl             anti-phage defence systems (PadLoc); with -df, agreed calls are drawn once
```

### Figures and output

```
-f FILE         your own figure table (section 8)
-nf             no figures; draw later with flags3 figures
-fh 16383       split a figure into parts above this height in pixels
-no             classic style: leave the number out of a gene too small for it
-pdf            a PDF beside every SVG
-c 8            worker threads (default: all cores)
--tools FILE    a tool table other than ~/.flags3/tools_table.tsv
```

### A full example

```
flags3 run -i list.txt -u you@example.org -g 4 -sr 20000 -cr -t -to -d -ss -gn -df -pdf -o myrun -nt
```

## 5. The run directory

```
myrun/
    run.tsv           version, start time, command line
    config.tsv        every option, resolved; a rerun of a stage reads it
    console.log       everything printed, plus every external command and its output
    input/            copies of the input lists
    blast/            hits.tsv, accessions.txt (only with -bi or a BLAST line)
    fetch/            genomes.tsv (which files were used), queries.tsv, failures.tsv
    extract/          genes.tsv, windows.tsv, rows.tsv, unmatched.tsv, range_report.tsv,
                      proteins.faa, queries.faa, rna.fna
    cluster/          families.tsv, hits.tsv, annotations.tsv
    cluster_rna/      the same for RNA genes
    tree/             tree.nwk, alignment.aln, trimmed.aln, leaves.tsv, commands.txt
    domains/          domains.tsv, annotations.tsv
    features/         features.tsv, annotations.tsv
    sismis/           secretion.tsv, diagnostics.tsv, annotations.tsv, raw/
    genomad/          mobile_elements.tsv, diagnostics.tsv, annotations.tsv, raw/
    defence/          defence.tsv, diagnostics.tsv, annotations.tsv, raw/
    report/           the tables for reading (section 6)
    figures/          the SVGs and PDFs (section 7)
```

Every stage directory has a `status.tsv` saying `ok` or `failed`, with the
time it took and the error if any. An optional stage that fails (a tool not
installed, no queries for a tree) is recorded there and the run continues;
`report/run_summary.txt` lists them all.

### Rerunning a stage

Every stage is a subcommand that takes a run directory and re-does only its
own part, reading the options the run was made with, plus any you give:

```
flags3 figures myrun -f figures.tsv        # redraw
flags3 cluster myrun -cm mmseqs_cluster    # recluster, then: flags3 report myrun; flags3 figures myrun
flags3 domains myrun -db defensefinder=~/.flags3/db/defensefinder/profiles -hc defensefinder=0.7,0.5
flags3 extract myrun -g 8                  # then cluster, tree, domains …, report, figures
flags3 tree myrun -iq
```

A flag the run had on can be switched off for the rerun with `--no-<option>`,
e.g. `flags3 features myrun -lth --no-signalp` reruns only the local
transmembrane prediction.

A stage deletes and rewrites its own directory only. Stages that read it
(`report`, `figures`, and after `extract` all analysis stages) are not
rerun automatically; run them yourself in that order.

## 6. Reading the results

`report/neighbourhoods.tsv` is the main table: one line per gene per row.

| column | meaning |
|---|---|
| row_id | `query|assembly` |
| offset | position relative to the query in the query's reading direction; negative is upstream |
| accession | protein accession, or `pseudogene*` |
| family | the family label, as on the figure: `Q1` holds a query, `R1` is an RNA family, plain numbers are shared flanking families, `-` is a singleton |
| strand | the gene's strand on the contig (the figure flips rows so the query points right) |
| domains | domain names on the protein, in order, when `-d` was on |

`report/queries.tsv` says what happened to every input query: the
assemblies it was resolved to, how many rows it produced, and why none if
none. `report/families.tsv` lists every family with a representative
product. `flanking.faa`, `queries.faa` and `all.faa` carry the sequences;
every FASTA id is the identifier used in the tables and the product or
species follows after a space.

`report/legacy/` holds the FlaGs2-era files with the same names and formats
as before (`<run>_operon.tsv`, `_outdesc.txt`, `_speciesInfo.txt`,
`_QueryStatus.txt`, `_flankgene_Report.log`, `_clusters.tsv`,
`_rangeReport.tsv`, `_accessionIssues.txt` and the three FASTAs) for
scripts that depend on them.

The stage tables are one level down. `cluster/hits.tsv` says which
sequences each protein hit; `domains/domains.tsv` has every hit with
E-value, clan and InterPro fields; `sismis/secretion.tsv`,
`genomad/mobile_elements.tsv` and `defence/defence.tsv` give every called
system with absolute contig coordinates, the rows it overlaps, and the
tool's own columns.

## 7. Reading the figures

Each row is one query in one genome. The query is drawn pointing right,
with a black outline; the row is mirrored when the query is on the minus
strand. Genes are coloured by family and labelled with the family label
above them; grey genes belong to no shared family. RNA genes have a green
outline, pseudogenes a blue one.

| figure | what it shows |
|---|---|
| `neighbors` | the neighbourhoods |
| `tree` | the same, aligned to the tree of the queries, with bootstrap values in red |
| `domains` | domain wedges inside the genes, numbered; TM regions as red hatching, signal peptides as a black triangle; family numbers prefixed `G` so they don't clash with domain numbers |
| `sismis`, `genomad`, `defence` | genes in grey, the called systems as coloured bands under the row with codes `S1`, `M1`, `D1` |
| `all-in-one` | everything, one band lane per tool |
| `classic`, `tree_classic` | the original FlaGs look |

A figure is drawn only when a stage it shows produced something. Every
figure has a legend panel per layer; a domain panel can be long.

## 8. Your own figures

Figures are rows of a table. Copy the shipped one and edit it:

```
python -c "import flags3.data, importlib.resources as r; print(r.files(flags3.data) / 'visualisation_table.tsv')"
cp <that path> figures.tsv
flags3 figures myrun -f figures.tsv
```

Columns:

| column | values |
|---|---|
| name | file name of the figure |
| layers | stages to draw, comma-separated: `cluster`, `cluster_rna`, `domains`, `features`, `sismis`, `genomad`, `defence` |
| mode | `versatile` (default look) or `classic` (FlaGs) |
| tree | `1` to add the tree panel and order rows by it |
| numbers | `true` to label genes with their family |
| palette | `bright`, `pastel`, `classic`, `colourblind`, `monochrome` |
| monochrome | `true` to grey the gene fills so bands stand out |
| font_size, row_height, gene_height, bases_per_pixel, pad, domain_height, label_step, arrow_head, min_gene_width, band_opacity, tree_width | geometry; empty means the mode's default |

To pin single colours, put a `colours.tsv` in the run directory:

```
#stage	category	colour
cluster	family:3	#d40000
sismis	T3SS	#0066cc
```

Categories are the values in each stage's `annotations.tsv`.

## 9. The tools table

External programs and clustering settings are rows of a tab-separated
table. The shipped one is the default; `flags3 install` writes your
machine's paths into `~/.flags3/tools_table.tsv`, which is read
automatically; `--tools FILE` points at another.

```
#name	command	directory	scan_range	engine	options
jackhmmer				jackhmmer	iterations=3;incE=1e-3;chunk=100;chunks_per_worker=16
mmseqs_cluster	/home/me/.flags3/envs/mmseqs/bin/mmseqs easy-search {in} {in} {out} {tmp} ...		mmseqs	sensitivity=7.5;coverage=0.5
sismis	/home/me/.flags3/envs/sismis/bin/sismis run -g {in} -o {out}
```

To change a clustering threshold, edit the `options` of its row; to use a
tool installed elsewhere, put its absolute path in `command`. Rows can be
edited but not invented: the set of tool names is fixed by the package.

## 10. When something goes wrong

- **"no query could be resolved"**: check `fetch/queries.tsv` for the
  reason per query and `fetch/failures.tsv` for network errors. IPG needs
  a valid email in `-u`.
- **A stage says failed**: `<stage>/status.tsv` has the error and
  `<stage>/traceback.txt` the full trace. Most often the tool is not
  installed: `flags3 install <name>`.
- **A genome file is broken** (an interrupted download): it is detected
  and downloaded again on the next run. To force it, delete the file from
  `~/.flags3/genomes/`.
- **Figures look wrong after rerunning a stage**: `report` and `figures`
  read the other stages' tables and need rerunning too.
- **Everything printed** is in `console.log`, including every external
  command and what it wrote to stderr.
