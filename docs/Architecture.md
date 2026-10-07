# FlaGs3 — Architecture

Developer-facing. Records the decisions the rewrite is built on; the sections
grow as modules are ported. For usage see `User_Guide.md`.

## Why a rewrite

FlaGs3 2.3.0 grew by attaching modules to a two-module pipeline. Its data bus
was the live `NeighborhoodExtractor` object, read by five modules; its
pipeline was a 530-line `main()`; the three band-producing tools were one
module written three times; the renderer drew every layer in one method; and
figures were produced by shelling out to a second script. None of it could be
rerun on an existing output directory.

The rewrite keeps what worked (network access in one module, families as
connected components with pluggable edge-finding, the two configuration
tables) and replaces the glue with three contracts.

## The three contracts

### Identity

- A **row** is one query in one genome: `row_id = query|assembly`.
- A **gene** is identified by its accession (or locus tag for local and
  MGnify genomes). One protein in five rows is scanned once and drawn five
  times.
- A **bp subject** is `assembly|contig`.

### Coordinate spaces

Exactly two: `bp`, an absolute position on a contig, and `aa`, a residue on a
protein. Nothing stored is row-relative, so a change of window size after the
fact means re-extract and re-render, never rescan. Row-relative offsets exist
only inside the renderer's layout.

### Annotation

Every analytical stage writes one `annotations.tsv` in a single schema:

| column | meaning |
|---|---|
| subject | accession, or `assembly\|contig` for bp |
| space | `aa`, `bp`, or `-` for a whole-gene fact |
| start, end | 1-based inclusive; `-` when space is `-` |
| kind | glyph class: `fill`, `outline`, `band`, `wedge`, `segment`, `triangle` |
| category | what the visualisation table maps to a colour |
| label | text shown in the figure |
| tool | which program made the call |
| score | e-value, probability, or `-` |

`kind` is a closed vocabulary owned by the renderer. Stages emit facts and a
glyph class; colour, opacity, lane and font are rendering decisions made from
`category` at draw time. Clustering fits the same schema (`fill` with
`category = family:N`), so the gene layer is not special.

The extractor's `genes.tsv` and `windows.tsv` are the only other schemas.
Every table is TSV with a plain header line, `-` for missing, `true`/`false`
for booleans. `flags3.schema` holds the dataclasses and does all reading and
writing; no module parses TSV itself.

## Run directory

```
<run>/
    run.tsv        version, start time, command line
    config.tsv     every resolved option as key<TAB>value
    console.log
    input/         copies of the input lists and blast input
    extract/       genes.tsv  windows.tsv  proteins.faa  rna.fna  range_report.tsv
    <stage>/       annotations.tsv  status.tsv  raw/
    figures/
    report/
```

A stage owns its directory and nothing else. Running a stage deletes and
recreates its directory, so a rerun is an overwrite by construction.
`status.tsv` records `status` (`running`, `ok`, `failed`), start time,
seconds, and the error for a failure; `traceback.txt` beside it holds the
full trace. `report/` and `figures/` are derived from the stage directories
and are rebuilt by their own stages.

`config.tsv` is written once by `flags3 run` from the command line. A stage
run later on the directory (`flags3 domains <run>`) reads its options from
there; flags given on that command line override the file for that run only.
Flags remain the interface; the file exists so a rerun needs none.

## Stages

```python
class Stage:
	name = ""
	requires = ()
	optional = False
	def wanted(self, config) -> bool
	def run(self, run, config, out) -> None
```

`requires` names the stages whose `status.tsv` must say `ok`; the runner
refuses to start otherwise rather than guessing. `wanted` is the opt-in
switch, read from config. `optional` decides whether a failure stops the run
or is recorded and stepped over. The runner (`flags3.stage.Runner`) does the
timing, the status file, the directory reset and the console lines, so a
stage body contains only its work.

Stages import `schema`, `run`, `tools` and `log`; never each other. The
pipeline order is one list in the CLI's `run` subcommand. Stages run
in-process; only external tools are subprocesses.

## Package

```
flags3/
    __init__.py     VERSION
    cli.py          subcommands: run, one per stage, redraw, install
    run.py          RunDir, Config, Status
    schema.py       Gene, Window, Annotation
    stage.py        Stage, Runner
    tools.py        external command table
    log.py          console transcript
    stages/         one file per stage
    render/         layout, one function per glyph kind, figure assembly, pdf
    data/           default tools_table.tsv, visualisation_table.tsv, font metrics
```

Python 3.11 or later. Tabs. No docstrings; this file and `User_Guide.md`
carry the explanations.

## Decided, not yet built

- Rendering: a `Layout` is computed once from `genes.tsv` (row → y, bp → x,
  aa → x inside an arrow, clip paths). Each glyph kind is a pure function of
  the layout and its annotations, returning one `<g id="layer-<stage>">`.
- Configuration tables: the package ships defaults; a user copy is merged at
  load time by adding missing rows and never touching existing ones. That
  merge, with `pip install --upgrade`, is the update mechanism.
- Installers: one `flags3 install <tool>` driven by a registry, replacing the
  eight shell scripts.
- Report: new human tables designed around the schema, plus a `legacy/` set
  reproducing the FlaGs2-era files for users who rely on them.
- Output-file compatibility with 2.3.0 is not preserved.

## Fetch and extract

`fetch` turns the input lists into two tables: `genomes.tsv` (assembly,
source, and the four file paths) and `queries.tsv` (query, assembly, the
accessions IPG says the query goes by in that assembly — identical
sequence by IPG's definition, never aligned locally — a status, and `via`:
how the query was placed), plus
`failures.tsv` for everything that went wrong on the way (unreachable IPG
chunks, download errors, remaps, cross-database exclusions). Resolution and
download are one stage because both answer "where are this query's genome
files".

### The genome directory is a cache

2.3.0 had `-tmp` (download here, delete after), `-k` (keep) and `-ul` (look
here, never download here). Rerunnable stages need the genomes to outlive
the run, so there is one directory, `-gd/--genomes`, default
`~/.flags3/genomes`. Fetch looks there first; what is missing is downloaded
into it; nothing is ever deleted. Files are downloaded to a `.part` name and
renamed into place, so two runs can share the directory without a lock: a
half-written file is never visible under its final name. Cached files this
run needs are checked before use (gzip streams are read to the end); a
truncated file is deleted and downloaded again. No checksums are kept.

`--no_cache` sets the run's genome directory to `<run>/genomes/` before
the config is written; fetch and every later stage read that path from
`config.tsv`, so nothing else knows the difference.

### Resolution order

- A bare protein goes to IPG, which chooses the assembly (`-m` many,
  RefSeq first, `-nc` to stay in the protein's own database); then the cache
  is checked for that assembly. 2.3.0's `-ul` did the reverse, taking
  whichever local genome contained the protein, which would let a growing
  cache change an assembly choice between runs. That behaviour is opt-in:
  under `--offline` for everything, or under `-lf/--local_first` for the
  bare queries the directory's protein files contain, the rest still going
  to IPG in the same run.
- A paired query uses its assembly as given. Under `-rm/--remap` IPG is asked
  about it too: if the protein is in that assembly under another accession,
  the alias lands in `queries.tsv` and extract matches on it; if IPG says the
  protein is not in that assembly at all, IPG's first choice replaces it and
  `failures.tsv` records the remap. That is 2.3.0's fetch-extract-fetch
  retry loop done in one pass.
- A version-less assembly (`GCF_000422825`) matches a versioned cache entry
  by prefix and the row is named after the cache entry.

### What gets downloaded

Each stage declares `needs` (`"rna"`, `"genome"`); `run` writes the union
over wanted stages into `config.tsv` as `fetch.slots`. GFF and protein FASTA
are always fetched. Fetch never learns a stage's name, and a stage never
learns how a file arrived. MGnify offers no RNA FASTA; extract cuts RNA
sequences from the genome FASTA when it has one.

`net.py` holds the rate limiter (5/s, 10/s with an API key), the retrying
session and the atomic download. `ncbi.py` holds IPG, the `XP_` BioProject
path and the partitioned FTP layout; `mgnify.py` the API v2 lookup. Tests
replace the two network seams (`IpgMapper._fetch_ipg`, `NcbiGenomes.listing`
and `.stream`) with recorded responses and run offline.

`extract` reads those two tables and nothing else. Per assembly it parses the
GFF once (`genome.GeneTable`: genes sorted by contig and start, a running
maximum of gene ends per contig so a bp range can walk left without scanning
the whole contig, and pseudogene/RNA handling inherited from 2.3.0), then cuts
one window per query and writes:

| file | content |
|---|---|
| `genes.tsv` | the neighbourhood, one `Gene` per row and gene |
| `windows.tsv` | per row: neighbourhood span and scan span, both absolute |
| `rows.tsv` | per row: the accession that matched, contig, strand, species |
| `unmatched.tsv` | queries with no row and why |
| `range_report.tsv` | how much contig was available and reached, up and down |
| `proteins.faa`, `queries.faa`, `rna.fna` | sequences; queries keyed by row id |
| `scan_genes.tsv`, `scan_proteins.faa` | genes inside the scan range, when `-sr` is set |

Three changes from 2.3.0, all deliberate:

- `strand` is the gene's true strand. 2.3.0 stored it flipped relative to the
  query so the renderer could ignore strands; that is a layout decision and
  now happens in the layout. `offset` is still query-relative (negative is
  upstream in the query's reading direction), because it is a fact about the
  neighbourhood, not about the drawing.
- The scan span is one per row, from `-sr` plus `-sm`, or the whole contig
  without `-sr`. The per-tool `scan_range` column of the tools table no longer
  changes it; a tool that must see the whole genome reads the genome file.
- RNA sequences are written whenever a source is available. Downloading the
  RNA FASTA is still opt-in at fetch time; extract writes what it finds.

Merging overlapping scan spans across rows of one assembly is the scanning
stage's business, not the extractor's.

## Clustering

A family is a connected component of a graph whose edges mean "these two
sequences are homologous"; only the edge-finding differs between methods,
which is why `cluster.py` has one `connected_components` and three small
`Clusterer` subclasses (jackhmmer and nhmmer through pyhmmer, mmseqs through
the tools table). Components are taken on the symmetrised graph, so a hit in
either direction joins two sequences and a sequence never lands in two
families. The method is a row of the tools table (`-cm`), and every setting
(iterations, inclusion E-value, chunking, mmseqs sensitivity) lives in that
row's `options` column; there is no `-e` or `-n` any more.

The jackhmmer row sets `incE` and `incdomE` to 1e-10, FlaGs2's values. FlaGs3
1.0 to 2.3.0 ran at `incE` 1e-3 with HMMER's default `incdomE` of 0.01, which
admits weak partial hits as family edges and lets marginal domains into the
iterated profile; on the fusion case that raised the question the two
settings gave identical families (every edge below 1e-42), but on promiscuous
domains the looser setting joins families FlaGs2 keeps apart. `-ce` sets
both thresholds together.

jackhmmer calls carry many queries at once (`chunk`, `chunks_per_worker`),
which is the fix for the per-call memory growth that 2.3.0 hit on dense
subfamilies. mmseqs writes its input and hit table under `cluster/raw/`.

Two stages share the code. `cluster` runs on `extract/proteins.faa`;
`cluster_rna` runs on `extract/rna.fna` under `--cluster_rna`, declares
`needs = ("rna",)`, and groups RNAs that have no sequence by product name.
Each writes:

| file | content |
|---|---|
| `families.tsv` | family number, label, size, occurrences across rows, members |
| `hits.tsv` | per accession: its family and the accessions it hit |
| `annotations.tsv` | one `fill` per member of a labelled family, `category = family:N` |

Labels follow 2.3.0: a family is labelled only if its members occur more
than once across rows; labelled families are ordered by occurrence count;
a family holding a query is `Q<n>`, an RNA family `R<n>`, the rest `1, 2, …`.
Within a query family the proteins that were in the input list carry `Q`
and every other member — a homologue met as someone else's flanking gene —
carries `q`; same colour and number, the case says which arrows are queries.
Unlabelled families get no annotation, so the renderer draws them in its
default style without knowing why. Family numbers in `families.tsv` are
ordered by size and are stable identifiers; labels are what the figure shows.

Verified against the 2.3.0 reference run: the same 130 families, the same
label on every gene, and the same hit set for every accession.

## Tree

`tree` aligns `extract/queries.faa` with mafft, trims with trimal (or an
internal gap-threshold trim when trimal is absent and the mode is `gt`),
and builds the tree with VeryFastTree or, under `--iqtree`, IQ-TREE with
ModelFinder and 1000 ultrafast bootstraps. It is wanted under `--tree`,
`--iqtree` or `--tree_order` and is optional: fewer than three queries or a
missing binary is recorded in `tree/status.tsv` and the run continues.

Outputs: `tree.nwk`, `alignment.aln`, `trimmed.aln`, `commands.txt`, and
`leaves.tsv`, the midpoint-rooted, ladderized leaf order as `position`,
`row_id`. Anything that wants rows in tree order reads `leaves.tsv`; nothing
re-parses the Newick to decide an order. Intermediate files live in
`tree/raw/`.

## The home directory

Everything FlaGs3 keeps between runs lives under `~/.flags3/` (`home.py`):

```
~/.flags3/
    genomes/            the genome cache (-gd overrides)
    db/pfam/            Pfam-A.hmm (pressed), Pfam-A.clans.tsv.gz
    db/interpro/        interpro_metadata_processed.tsv(.gz), optional
    tools/              binaries the installer downloads (mmseqs, ...)
    envs/               conda environments the installer creates
    tools_table.tsv     the user's overrides, written by the installer;
                        read by default when it exists (--tools overrides)
```

Defaults that used to be "next to FlaGs3.py" or "./pfam_db" now point here,
so a run behaves the same from any working directory and the package
directory is never written to.

## Domains

`domains` (`-d`) scans `extract/proteins.faa` with pyhmmer against one or
more HMM sources (`-db`, repeatable; a `.hmm` file, pressed or not, or a
directory of `.hmm` files such as DefenseFinder's profiles). Models that
carry a gathering threshold are scored by it; the rest use `-e`.
`-hc` drops hits that cover too little of the protein or of the model, per
source or for all. Profile names with a `__` separator (`Cas__Cas9`) get a
`group` (`Cas`).

Outputs: `domains.tsv` (protein, database, domain, group, Pfam accession,
clan, span, E-value, and the InterPro columns when the metadata table is
present) and `annotations.tsv` with one `wedge` in `aa` space per hit.
`category` is the clan when `-cl` is given, otherwise the group, so the
renderer colours by clan or by domain without knowing which was chosen;
`label` is the profile name; `tool` is the source's name.

The family column 2.3.0 put in `_domains.tsv` is not here: joining families
onto domains is the report stage's job, and this stage does not depend on
clustering.

## Installing tools: `flags3 install`

One command replaces `build.sh` and the eight installer scripts. It manages
what FlaGs3 owns under `~/.flags3/` and nothing else: the Python
environment FlaGs3 itself runs in is the user's (`pip install flags3`).

`flags3 install` lists the components and their state; `flags3 install
core pfam`, `flags3 install --all` (everything downloadable), `flags3 install
signalp ~/signalp-6.tar.gz` (a component that takes a file), `--update`
refreshes installed ones.

A tool already on PATH is used as it is: `core` locates each binary
through its tools-table row and then `shutil.which` before building
anything, which is how a conda environment made from the shipped
`environment.yml`, a Homebrew install or a cluster module are picked up.
What is missing comes from bioconda through micromamba — the one on PATH
if present, otherwise one 18 MB static binary the installer downloads into
`~/.flags3/tools/`; an existing conda is neither used nor touched. Each
built component is one environment under `~/.flags3/envs/`, so components
update independently and a broken one is deleted by removing a directory.
MMseqs2 is the static release binary from GitHub, as in 2.3.0, not a conda
package. On Darwin `DYLD_LIBRARY_PATH` is unset for the installer's
subprocesses, as `build.sh` did; environments are built for the native
platform, and only if that solve fails on an arm64 Mac is the same
command retried with `--platform osx-64`. Databases go under `~/.flags3/db/`. Every component finishes by
writing its rows into `~/.flags3/tools_table.tsv` with absolute binary paths
through `Tools`, so nothing depends on an environment being activated. The
run reads that table by default.

| component | what it does |
|---|---|
| `core` | one env with mafft, trimal, VeryFastTree, IQ-TREE, BLAST+; rewrites their rows |
| `pfam` | Pfam-A.hmm and clans from EBI, pressed with pyhmmer; or from a file you give |
| `interpro` | copies the lab's metadata table into place (not downloadable) |
| `defence-hmm` | DefenseFinder profiles from the latest GitHub release, for `-db` |
| `mmseqs` | MMseqs2 release binary under `~/.flags3/tools/`; rewrites the three mmseqs rows |
| `genomad`, `padloc` | bioconda env plus each tool's own database download |
| `defensefinder` | env with Python 3.10 and hmmer, then `pip install mdmparis-defense-finder` as 2.3.0 did (the bioconda package resolves onto a Python without a working MacSyFinder), then `defense-finder update` |
| `sismis` | env with hmmer, `pip install sismis` |
| `deeptmhmm` | DeepTMHMM2 (MIT, `fteufel/DeepTMHMM2`) in its own env via pip from GitHub, weights fetched by a warm-up prediction; replaces the licensed DeepTMHMM 1.0 package |
| `signalp` | licensed: env built from the package you downloaded |

A component is a class with `installed()`, `install(source)` and
`update()`; the conda-backed ones differ only in package list, binary name
and command template (`CondaTool`). Components that take a file accept it
in place of a download, which is also how the tests install a one-model
Pfam and a two-file DefenseFinder tarball without network.

The geNomad database path goes into the `options` column of its row as
`db=…`, the same place clustering keeps its settings, rather than the
`directory` column 2.3.0 used for it.

Not verified in the sandbox: bioconda package resolution (the sandbox cannot
reach conda channels), the Pfam and DefenseFinder downloads, and whether
Sismis' PyPI package needs more than hmmer beside it. The micromamba
bootstrap itself was run for real.

## Scanning windows: sismis, genomad

Tools that read sequence rather than genes get one window per row, cut from
the genome FASTA. `windows.tsv` carries three spans per row: the
neighbourhood (`lo..hi`), the analysed span (`scan_lo..scan_hi`, the query
plus `-sr`), and the cut span (`cut_lo..cut_hi`, the analysed span plus
`-sm`, so a system straddling the edge is called whole). Overlapping cut
spans on one contig are merged (`windows.merge`), and `windows.Batches`
writes them as FASTA batch files of at most 200 Mb each with records named
`w0, w1, …` and a `cuts.tsv` beside them that maps each record back to
assembly, contig and offsets. A hit in record coordinates is placed with
`Cut.place`, which also decides `full` (inside the analysed span) or
`partial` (reaching into the margin).

Each scanner can have its own span — `-ts tool=SPAN` per run or `span=`
on its row — in bp each side of the query or `contig` for the whole
contig. The stage recomputes that scanner's windows from the query
coordinates in `genes.tsv` (`windows.respan`, the same arithmetic as the
extractor's `scan_span`, with `-sm` as margin) before the ordinary
merge-and-batch path; a tool without a span follows `-sr`. Other contigs
of the assembly are never scanned: geNomad classifies each sequence on
its own, so they would add nothing to the query's contig. Defence has no
span of its own: it follows `-sr` or the neighbourhood, since systems more
than a few genes from the query are not what the figure is about.

`scan.WindowScan` is the stage base: it resolves genome files, cuts and
batches, runs the tool per batch, places the hits, and writes
`annotations.tsv` (`band` in `bp` space, `category` = the called type),
a report table with the tool's own columns passed through under a prefix,
and `diagnostics.tsv` per assembly. A subclass supplies the tool name and
`scan_batch`, which runs the command and parses its output into `Hit`s in
record coordinates. `Sismis` is one such subclass; `Genomad` will be the
second. Each stage cuts its own windows under its own `raw/`; 2.3.0 shared
the cut files between tools, which saved a few seconds and broke the rule
that a stage owns only its directory.

Verified: the 27 cut windows are identical in coordinates and sequence to
the 2.3.0 reference batch, and placing the reference run's two Sismis
clusters through `Cut.place` reproduces their absolute coordinates and
coverage.

## geNomad, defence systems, TM and signal peptides

`genomad` (`-gn`) is a second `WindowScan`: the same cut windows, a
different parser. The band's `category` — hence its code, colour, legend
entry and summary line — is what geNomad called: the last named viral rank
(`Caudoviricetes`) or the plasmid type (`conjugative plasmid (AMR)`), as
in 2.3.0; the coarse kind (`virus`/`plasmid`) is a column of
`mobile_elements.tsv`. Its row carries `batch_mb` (default 50), the size of
each FASTA batch it is called on, and `splits`, passed as geNomad's own
`--splits` to divide the marker search. Batching caps the query-side
memory (translated ORFs, per-thread k-mer tables, prefilter result
buffers); `splits` caps the DB-index side. A 766-window, 89 Mb batch at
22 threads peaked around 20 GB and was killed by the OOM killer on a
laptop with a browser open; 50 Mb is chosen from that measurement, and a
signal-9 death is reported with a pointer to the two options. The
database path is the `db=` option of the genomad row, which
`flags3 install genomad` writes, or `-gdb`.

`defence` (`-df`, `-pl`) runs on genes rather than windows. With `-sr`
set it takes the protein-coding genes of the scan span (`scan_genes.tsv`,
`scan_proteins.faa`), so a system that straddles the edge of the drawn
neighbourhood is called whole, as it is for the sequence scanners; without
`-sr` it uses the neighbourhood genes, as 2.3.0 did. Each row's set
becomes a synthetic replicon `r<n>` in one GFF
and one FASTA under `defence/raw/`, DefenseFinder and PadLoc run on those,
and their calls are mapped back through the `r<n>_<order>` tags to real
accessions and coordinates. A system both tools call on the same genes with
the same name is one row with `called_by = defensefinder,padloc`. Outputs:
`defence.tsv`, `diagnostics.tsv` (one line per tool), and `band`
annotations whose `tool` column carries the caller list.

`features` (`-th`, `-sp`) runs DeepTMHMM and SignalP 6 on the neighbourhood
proteins, on BioLib by default or locally with `-lth`/`-lsp` and the tools
table. On BioLib each tool is one job for up to 2,000 (DeepTMHMM) or 1,000
(SignalP) proteins, submitted with `blocking=False` so both tools queue at
once; 2.3.0 sent batches of 25, each paying a queue wait. Scanners have a
`start`/`finish` pair so a local run and a cloud job share one shape. Both
the cloud app (`DTU/DeepTMHMM2`, arguments `<fasta> results --simplify-io`
as its `biolib/run.sh` expects; `--tmhmm_app DTU/DeepTMHMM` selects the
1.0 app with its `--fasta` syntax) and the local predictor are DeepTMHMM2
with `--simplify-io`, whose 3-line output
uses `S`, `M`, `B` as DeepTMHMM 1.0 did, plus `R` (reentrant), `F`
(interfacial) and `>` (transit peptide), which are not drawn. DeepTMHMM's
own signal-peptide calls are used only when SignalP is not requested. Outputs: `features.tsv` and annotations in `aa` space:
`segment` for a transmembrane region (`category = tm`), `triangle` for a
signal peptide (`category = signal`).

Every optional stage fails softly: a missing binary, database or package is
recorded in the stage's `status.tsv` and the run continues. With all ten
stages registered the pipeline is `fetch, extract, cluster, cluster_rna,
tree, domains, features, sismis, genomad, defence`; `fetch.slots` comes out
as `genome,rna` when the window scanners and RNA clustering are on.

## Figures

`figures` is the last stage and the one `flags3 figures <run>` reruns. It
reads only the run directory: `extract/genes.tsv`, `windows.tsv` and
`rows.tsv` for geometry and labels, `tree/leaves.tsv` and `tree.nwk` for
order and the tree panel, and every `<stage>/annotations.tsv` whose stage
is `ok`. It imports no stage module; a stage that isn't there is a layer
that isn't drawn.

### Layout, once

`render/layout.py` computes, for a figure's rows, everything the layers
share: y per row, x per contig base (`x(row, bp)`), the arrow polygon per
gene, the residue-to-pixel mapping inside an arrow (`x_in_gene`), and one
clip path per arrow. The strand flip lives here: a row whose query is on
`-` is mirrored so the query points right, and `drawn_strand` says what a
gene looks like after the mirror. Rows are in tree order when the figure
carries a tree, else input order (extract keeps input order).

### Layers, one function per glyph kind

`render/layers.py`: `draw_genes` (the arrows, filled from the cluster
layer's `fill` annotations or grey), `draw_bands` (`band`: a hatched
rectangle across the called span drawn behind the genes, 5 px strokes with
5 px gaps in the category colour, written as plain `<line>` elements
inside a clip path rather than as an SVG `<pattern>`, because cairosvg
rasterises patterns into a washed-out tile in PDFs while every renderer
draws clipped lines the same way; each band stage hatches at its own angle
— sismis 45°, genomad −45°, defence horizontal — so overlapping calls
cross rather than hide each other, and the row height does not grow with
the number of tools; the band codes `S1`, `M1`, `D1` of a row are written
after its last gene in their colours),
`draw_wedges` (`wedge`, clipped to the arrow), `draw_features` (`segment`
hatch and `triangle`), `draw_outlines` (re-stroke arrows that carry
overlays), the label placers, and the row labels. Each returns one
`<g id="layer-<stage>">`, so a layer can be found, hidden or replaced in the
SVG. Z-order is fixed: band hatches, genes, domain wedges, outlines, protein
features, text. Transmembrane hatching is clipped to its arrow; the
signal-peptide triangle is not clipped and is drawn above the outline at
the true N-terminal position, so it is never cut by the arrow's edge. A bp annotation is clipped to its row's neighbourhood span.

Labels above a row (family numbers and domain codes) are laid out in one
line by measured text width: a label that would overlap the previous one
is pushed right until it clears it. The gene area starts where
the longest *label-plus-its-own-genes* row needs it, not where the longest
label ends: a row with a long species name and few upstream genes lets its
label run under the space other rows use for genes, which keeps the figure
narrow.

### Style

`render/style.py` reads `visualisation_table.tsv` (shipped in
`flags3/data`, `-f` overrides): one figure per row, with `layers` (stages
to draw), `mode` (versatile or classic), `tree`, `numbers`, `palette`
(bright, pastel, classic, colourblind, monochrome), `monochrome` (grey gene
fills only, so bands stand out), and geometry columns that fall back to the
mode's defaults. A figure is drawn when any of its non-cluster layers has
data; a legend panel lists only the categories that figure (or part)
actually drew. `Colours` assigns each stage's categories — and so the codes `1, 2, …`,
`S1, S2, …` — in order of occurrence count by default (`numbering =
occurrence`; `appearance` restores first-seen order); the `bright` and `pastel` palettes step the hue by the
golden angle rather than sweeping it, so categories numbered next to each
other — which on a domain figure are usually neighbours on the same gene —
get hues far apart, and after seven entries the value drops so the second
lap is told from the first; `colours.tsv` in the run directory (`stage, category,
colour`) overrides single entries, opt-in.

Family numbers come from the cluster layer's labels and are placed above
the arrows in versatile mode, inside them in classic mode; on a figure that
also draws domains, plain numbers get a `G` prefix so they don't collide
with the domain codes. Domains are numbered per category and listed in the
Domains legend panel; bands carry `S1`, `M1`, `D1` codes with their panels.

Beside the figures, `figures/` gets companion files written by the same
stage from the same data: `legend.tsv` (stage, code, name, occurrences —
every domain and band code used across the run's figures),
`protein_clusters_legend.txt` and `rna_clusters_legend.txt` (the `_outdesc`
layout: label(occurrences), accession, product per cluster, members ordered
by subcluster letter, then bridges, then `?`, and by occurrence within a
letter), `systems.tsv` (one line
per drawn band: stage, code, system name, assembly, contig, span, the rows
it overlaps, the genes inside it as `accession(label)`, and their
products), and `systems_summary.txt` (per band stage: how many of the
run's neighbourhoods carry a system, then each type with the number of
neighbourhoods and of calls). Codes there are assigned over the whole run with the first
figure's numbering, so they match what the figures show.

Files in `figures/` and `report/` are named with the run's start time
(`RunDir.stamp`, `YYMMDD_HHMMSS`) unless the run was made with `-nt`, in
which case names are plain and stable for scripts. Stage tables are never
stamped; stages address each other by fixed names.

Figures taller than `-fh` are split by rows into SVG parts; a figure with a
tree panel is never split. The PDF of a split figure is one file rendered
from the whole figure, since PDF has no height limit — the limit is the
browser's.

A band is drawn behind its genes with a margin (`band_margin`, 6 px) above
and below the arrow. That margin is all that is visible of a system that
spans exactly its genes, so it is deliberately wider than the gene outline;
a legend entry with no visible band on the page means the band lies
entirely under a gene, not that it failed to render. `--pdf` converts with cairosvg (an optional
extra); the three other converters 2.3.0 tried are gone.

## Report

`report` joins the stage tables into what a person reads, and is the only
stage that reads several other stages' outputs. It runs before `figures`
and is rerun with `flags3 report <run>`.

`report/`:

| file | content |
|---|---|
| `neighbourhoods.tsv` | one line per row and gene, with family label, true strand, product and the domain names on that protein |
| `queries.tsv` | one line per input query: assemblies, rows found, status or reason |
| `families.tsv` | protein and RNA families with a representative product |
| `flanking.faa`, `queries.faa`, `all.faa` | sequences under the one header convention: id, space, description |
| `run_summary.txt` | version, command, counts, and each stage's status and seconds |

The FASTA convention: the id is the identifier the tables use (accession,
or `query|assembly` for a row); the product or species follows after a
space in the description field. No `|product` suffixes.

`report/legacy/<run>_*` reproduces the 2.3.0 files for people and scripts
that depend on them: `_operon.tsv` (with the query-relative strand of old),
`_clusters.tsv`, `_outdesc.txt`, `_speciesInfo.txt`, `_QueryStatus.txt`,
`_flankgene_Report.log`, `_accessionIssues.txt`, `_rangeReport.tsv`, and
the three FASTAs with unwrapped sequences and `accession|product` headers.
Checked against the reference run: byte-identical, except that queries are
listed in input order where 2.3.0 listed paired lines first.

Found while porting: 2.3.0's `_clusters.tsv` listed the RNA families twice
(once from nhmmer, once again from the product-name fallback). The legacy
file lists each family once.

## BLAST expansion

`blast` is the only stage before `fetch`. It is wanted when `-bi` names a
starting point or an input list has a line marked `BLAST`; `run` records
the latter in `config.tsv` as `blast_inline` so the stage can decide
without re-reading the lists. It fetches the sequence of an accession
through Entrez, runs QBLAST (remote) or the tools table's `blastp` row
(local), dedupes hits by version-less accession, and writes `hits.tsv` and
`accessions.txt`. `fetch` reads `accessions.txt` when the stage is `ok` and
appends its accessions as bare queries after the list's own, skipping any
the list already has. Nothing else knows BLAST happened; a rerun of `fetch`
on the directory picks the accessions up again.

## Background stages

A stage that mostly waits on the network can declare `background = True`.
The runner starts it in a thread once its requirements are met and goes on
with the next stages; a stage that declares `barrier = True` (report and
figures) makes the runner join every background thread first. Only
`features` is a background stage: DeepTMHMM and SignalP jobs queue on
BioLib while clustering, the tree, domains and the window scanners use the
cores. CPU-bound stages are not run concurrently with each other on
purpose: each already takes all cores, and running two at once would only
add memory pressure. Rerunning a single stage from the command line is
always synchronous.

## Subfamilies

A family is a connected component of "these two align", which is the
sensitive answer and also the one a fusion protein can lie to: one protein
carrying a ThiF domain and a UvrA/ABC pair joins two unrelated families
into one. Rather than cut the graph by shape (articulation points) or by a
community model (Louvain), the cluster stage uses a fact every hit already
carries — how much of each sequence the alignment covers — and asks the
same question twice at two strictnesses.

- The **family** is unchanged: any hit is an edge.
- A **core group** is a connected component of the full-length graph inside
  the family: an edge counts only when the alignment covers at least
  `subcov` (tools-table option, default 0.6; `-sc` overrides, `0` switches
  the whole thing off) of *both* sequences. Groups of one are not cores.
- A core group whose every member is **contained** in at least half of
  another group's members (the alignment covers ≥ `subcov` of the smaller
  protein) is folded into that group: fragments and truncations join what
  they are fragments of instead of forming a group of their own.
- A family with two or more cores after folding is split. Cores are
  lettered by size (`a`, `b`, …); a member's label is its family label plus
  its letters: its own core, plus every other core it *contains* a member
  of — unless 95% or more of its core-mates hold that core too, in which
  case the held group is a fragment of what the whole core shares and the
  letter would say nothing — plus, for a protein in no core, every core
  that contains it. A held core that partitions the holder's core is kept
  on purpose: in the run that raised it, `Q1af` (136) against `Q1a` (18)
  turned out to be ThiF proteins with and without a ~110-residue
  N-terminal extension, which the Pfam scan alone does not show. A
  ThiF–UvrA fusion is therefore `Q1ab`: full-length with the UvrA core and
  holding whole ThiF proteins. A protein with no letters at all is `?` in
  `families.tsv` and keeps the plain family label.

A protein at the end of a contig — every one of its occurrences in the run
is the first or last gene of its contig (`contig_edge` in `genes.tsv`) —
may be incomplete, so in a split family it is labelled `?` regardless of
its alignments; it still helps form the cores. A protein seen mid-contig
anywhere is complete and is labelled normally.

`families.tsv` gains `subfamilies` (`a:29;b:16`) and `bridges`
(`WP_260604591.1:ab`); `hits.tsv` gains `full_hits`. Annotations carry
`category = family:1/a` for core members and `family:1` for bridges, and
the renderer draws subfamilies as luminance steps of the family colour
(`a` the colour, `b` lighter, `c` lighter still), so a figure still reads
"all G1" from a distance and "two kinds" up close.

Verified on the collaborators' run that raised it: family Q1 (232 members,
one fusion) comes out as `a` = 158 ThiF proteins plus 3 fragments, `b` =
the 50 UvrA/ABC proteins plus the 20 unannotated ones, and `WP_260604591.1`
as `Q1ab`; the 46-protein subset gives `a:29`, `b:16`, `Q1ab`.

### The user tools table holds locations, not settings

Until 3.0.0rc1 `flags3 install` copied every row of the shipped table into
`~/.flags3/tools_table.tsv`, search options included, so a shipped default
changed later never reached a machine that had installed before — the
jackhmmer thresholds stayed at 1e-3 for one user after they had become
1e-10 in the package, and a 2,571-member family chained through 17 edges in
exactly that band was the symptom. The user table is now a table of
overrides: `write` emits only columns that differ from the shipped table, and
when the installer writes it, only tool commands, directories and the
options the installer itself sets (`db=`); any other option found there is
dropped with a note. Search settings therefore come from the shipped table,
from `-n`/`-ce`/`-sc`, or from a hand edit made after the last install.


## Home directory

`home.py` resolves the FlaGs3 home once at import: `$FLAGS3_HOME` if set,
otherwise `~/.flags3`. Everything below derives from it, so a server whose
`$HOME` cannot hold a 100 GB cache sets one variable and installs normally.


## Batch

`flags3 batch DIR -o OUT -- options` is deliberately the smallest thing
that removes the by-hand reload: a sorted loop over the files in `DIR`, one
`flags3 run -i file -o OUT/stem -nt options` per file as a subprocess, so a
crash in one run cannot take the loop down. "Done" is `report/status.tsv`
saying `ok`; done runs are skipped, which is what makes a killed batch
restartable with the same command; an incomplete run directory is reported
and redone only with `--retry`, so a batch never silently deletes work. No
parallelism — one run already uses every core — and no daemon: `nohup` or
`tmux` is the daemon.
