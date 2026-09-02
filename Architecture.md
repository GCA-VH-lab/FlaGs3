# FlaGs3 — Architecture Overview

A developer-facing map of how the pipeline fits together. For usage, see
`User_Guide.md`.

---

## Shape of the program

A linear pipeline. Each stage takes the previous stage's output and adds to it;
nothing loops back.

```
input list
    ↓  AccessionListReader
protein accessions (+ optional genome)
    ↓  LocalGenomeResolver → ProteinAssemblyMapper
protein → [genome ids]
    ↓  AssemblyDownloader / MgnifyGenomeDownloader
genome id → GenomeFiles(gff, faa, rna, genome)
    ↓  NeighborhoodExtractor
[FlankingGene] + protein/RNA sequence tables
    ↓  NeighborhoodClusterer / RnaClusterer
families: [[accession, ...], ...]
    ↓  OperonView / NeighborhoodVisualizer / ReportWriter
SVG figures + TSV tables + FASTA sequences
```

Optional stages (`--tmhmm`, `--signalp`, `--sismis`) hang off the side after
extraction and feed extra layers into the figures.

---

## Modules

| Module | Responsibility | Imported |
|---|---|---|
| `FlaGs3.py` | data pipeline and CLI | always |
| `flags_pdf.py` | SVG to PDF, one backend of several | on `--pdf` |
| `flags_log.py` | debug output, shared by every module | always |
| `flags_view.py` | shared styling, `OperonView` renderer | always |
| `flags_redraw.py` | figure table, reload from an output dir, render dispatch, CLI | subprocess |
| `flags_tree.py` | MAFFT + trimming + VeryFastTree or IQ-TREE, tree figure | on `--tree`/`--iqtree` |
| `flags_domains.py` | pyhmmer domain scan | on `--domains` |
| `flags_features.py` | DeepTMHMM / SignalP via BioLib | on `--tmhmm`/`--signalp` |
| `flags_secretion.py` | Sismis secretion-system scan | on `--sismis` |

Optional modules are imported lazily inside the branch that needs them, so a
default run never touches `mafft`, `pybiolib` or `sismis` — and a missing
dependency degrades to a warning rather than an import error at startup.

---

## Core types

Three small carriers move data between stages:

```python
GenomeFiles(gff, faa, rna, genome)   # paths; any may be None
FlankingGene(accession, strand, start, end, product, offset, query, is_rna, contig)
families: List[List[str]]            # each inner list is one family's accessions
```

`GenomeFiles` is what every genome source returns, whether downloaded from NCBI,
downloaded from MGnify, or found on disk. Downstream code never learns where a
genome came from.

`FlankingGene.offset` is signed distance from the query (`0` = query), already
flipped when the query is on the minus strand, so renderers can lay rows out
without knowing about strands. `FlankingGene.query` holds a **row id**
(`protein|genome`), not a bare protein — one protein in several genomes makes
several independent rows.

---

## Resolution and download

Genomes are resolved in priority order: local directory, then NCBI. Anything
matching `MGYG\d+` is routed to MGnify instead.

**`ProteinAssemblyMapper`** turns bare protein accessions into genome ids via
NCBI's IPG database in one batched request. `XP_` proteins are not in IPG, so
they take a separate BioProject → assembly path. Under `--no_cross_db` it drops
assemblies from the other database (RefSeq protein → `GCF_`, INSDC protein →
`GCA_`) before truncating to `-m`. Filtering happens after the fetch rather than
in the query, which leaves the IPG round trip unchanged and lets
`dropped_cross_db` record what was excluded for the `-vb` summary. Assemblies
supplied in the input file bypass the mapper entirely and are never filtered.

**`_GenomeDownloader`** holds everything the two remote sources share: HTTP
session, retry policy, worker pool, rate limiter, per-file streaming. Subclasses
supply only their URL scheme:

- `AssemblyDownloader` — NCBI's partitioned FTP layout, so it must first list a
  directory to discover the versioned name, then fetch gzipped files from it.
- `MgnifyGenomeDownloader` — one API call returns direct URLs; files are plain
  text. MGnify has no RNA-only FASTA, so the RNA slot stays empty and the
  extractor cuts RNA sequences out of the genome FASTA instead.

Each assembly's files are fetched concurrently, and assemblies are fetched
concurrently with each other.

**`RateLimiter`** caps requests per second across all threads by handing out
timed slots. This is deliberately separate from the worker count: more threads
or a faster machine would otherwise mean a faster request rate and rejected
downloads. The cap is 5/s, or 10/s when an NCBI API key is supplied.

---

## Extraction

`NeighborhoodExtractor` parses GFF and FASTA into gene records, then takes a
window of `±flank` genes around the query.

The GFF parser handles two annotation styles without being told which it has:

- **NCBI PGAP** — a `gene` feature followed by a `CDS` carrying the accession in
  `protein_id=`. The CDS attaches itself to the preceding gene record.
- **Prokka / Prodigal** (MGnify) — bare `CDS` features with no parent gene, where
  the locus tag *is* the accession. With no gene record to attach to, the CDS
  becomes its own record.

`_cds_accession` resolves the accession by trying `protein_id=`, then `ID=`
(stripping NCBI's `cds-` prefix), then the locus tag, and only then `Name=`.
`Name=` is last because Prodigal often sets it to a gene symbol such as `hisZ_2`,
which matches no FASTA record.

Genes are sorted by `(contig, start)` after parsing. Neighbours are selected by
list index, so file order must equal genomic order — true for NCBI, not for
MGnify, which writes all CDS features before all ncRNA features.

Parsed GFF and FASTA tables are cached per assembly, since one genome is usually
queried by several proteins. Genome FASTAs are the exception: they are loaded
only when an RNA actually needs slicing, and only one is held at a time, because
they are orders of magnitude larger than the other tables.

---

## Clustering

`NeighborhoodClusterer` runs jackhmmer with every flanking protein as a query
against all of them, building an adjacency map from the included hits, then takes
connected components as families. Clustering is symmetric by construction: if A
hits B, they end up in the same component regardless of direction. The adjacency
map is kept on the clusterer after the run and handed to `ReportWriter`, which
writes it as `_jackhits.tsv` — the only record of *why* two proteins share a
family.

`_outdesc.txt` orders its family blocks by total occurrences descending, which is
how main's file read. That ordering is local to the file: family *labels* come
from `family_numbers()` and stay consistent with the figures and `_clusters.tsv`,
so the labels in `_outdesc.txt` are not ascending. Main had both because it
derived its numbering from the same sort.

`RnaClusterer` mirrors this with nhmmer for RNA genes, falling back to grouping
by normalised product name when no nucleotide sequence is available.

---

## Tree building

`TreeBuilder` runs three stages: MAFFT `--auto`, a gap-threshold column trim, then
an inference engine. The trim keeps columns where at least `gap_threshold` of the
sequences carry a residue — the same rule as `trimal -gt`, which is what ete3's
`trimal01` ran in the old pipeline. It is implemented in `_trim` rather than
shelled out because `-gt` is a deterministic column filter and adding a binary
dependency for ten lines is not worth it; the heuristic modes (`-automated1`,
`-gappyout`) would need the real trimal. If trimming would empty the alignment,
the untrimmed one is used instead.

The engine is `veryfasttree` by default and `iqtree` under `--iqtree`
(ModelFinder, plus 1000 ultrafast bootstrap replicates when there are at least
four taxa — below that IQ-TREE refuses to bootstrap). The binary is looked up as
`iqtree3`, `iqtree2`, then `iqtree`, since distributions disagree on the name.

The trimmed alignment stays on the builder as `.alignment` so `main()` can write
it as `_tree.aln` without re-running anything.

`write_run_info()` runs immediately after the output path is resolved, before any
network or file work, so a crashed run still leaves `_runinfo.txt` and a copy of
the input list. It reads defaults off the parser rather than hardcoding them, so
the "set explicitly" split stays correct as options are added. Anything named in
`SECRET_ARGS`/`SECRET_FLAGS` is masked in both the option dump and the recorded
command line — `_redact_command_line` handles `--flag value` and `--flag=value`,
since only covering the first form leaks the key in the second.

## Concurrency, locking and diagnostics

`InstanceLock` guards `-tmp`, not the output directory: output paths are already
unique per run because of the timestamp, while the temporary directory is the
resource two runs actually fight over. The lock is created with `O_CREAT|O_EXCL`
and records the owning pid and host.

It is deliberately self-healing. A lock whose pid no longer exists on this host
is reclaimed automatically, and an unwritable lock path is a debug note rather
than a fatal error. The failure mode being avoided is pybiolib's, documented
below: a lock that outlives its owner and blocks every later run forever. Locks
from another host are never reclaimed, since liveness cannot be checked there.
Release is via `atexit`, so it covers exceptions and Ctrl-C; `SIGKILL` leaves the
file behind and the staleness check handles that on the next run.

`--debug` sets a module-level flag read by `debug()`, which writes to stderr and
optionally appends a traceback. `flags_tree.py` imports it lazily inside a
helper so the module keeps working standalone. Every external command goes
through `TreeBuilder._run`, which logs the command line and exit code — the
subprocess calls are the parts that fail silently.

### Repeating a BlastP run

`_blast_accessions.txt` holds the hits as a bare accession list that `-i` accepts,
so a run can be repeated without going back to NCBI. It is deduplicated on the
version-stripped accession and its two header lines start with `#`, which the
input reader skips.

### HMM sources

`HmmSource` describes one database: a `.hmm` file or a directory of them, a name
for the outputs, optional coverage cutoffs and the separator that splits a profile
name into a group. Nothing about it is DefenseFinder-specific -- that set is just a
directory of profiles named `System__Profile`, and any other collection works the
same way.

Models are scored by their own gathering threshold when they carry one and by
`-e` when they do not, so a database is split into two searches rather than forced
into one rule. Coverage cutoffs are per database because they mean opposite things
for the two kinds of model: DefenseFinder profiles are full-length proteins where a
partial match is noise, Pfam entries are domains where partial coverage is the
point.

Grouping keys the figure and legend on the part of the name before `__`, falling
back to the whole name, while `_domains.tsv` keeps database, profile and group so
the exact model that matched is recoverable.

### InterPro annotation of domains

`InterProAnnotator` joins Pfam hits to InterPro metadata on the **Pfam accession**,
which is why `DomainHit` carries `accession` alongside `name`: the source table has
no Pfam names, only `PF` ids in its `pfam_members` column, so a name-keyed join is
impossible. `pyhmmer` exposes the profile's `ACC` line as `query.accession`; the
version suffix is stripped for the lookup.

Only six short columns are kept. The table is ~85 MB with description fields up to
9 kB per entry, so pulling the description through into a TSV would make the output
unusable; the categorical summaries carry the useful signal. One `pfam_members`
cell can list several ids, each of which becomes its own key; if two InterPro
entries claim the same Pfam id the first wins and the count is reported under
`-vb`, so the join stays deterministic rather than depending on row order.

`resolve_interpro` treats a missing default and a missing explicit path
differently: the default is optional enrichment, so its absence is a debug note
rather than a warning on every run, while a path the user typed is a typo and
exits. It also looks next to `FlaGs3.py` and for a `.gz`, so the table works
whether it sits with the code or in the working directory.

### PDF output

`flags_pdf.py` tries `cairosvg`, then `svglib`+`reportlab`, then `rsvg-convert`,
then `inkscape`, and reports which it used under `--debug`. Availability is
checked with `importlib.util.find_spec` rather than a trial import, so probing
has no side effects. `cairosvg` is first because it embeds font subsets, so the
PDF renders identically on a machine without Liberation Sans; `svglib` is a pure
Python fallback needing no system cairo. Both were checked against real figures
for clip paths and fill opacity, which the domain and secretion overlays depend
on. A missing backend warns and leaves the SVGs, like every other optional tool.

### Fill and outline carry different facts

A gene can be both clustered and something else. The fill always carries the
family; the outline carries whatever else the gene is -- RNA, pseudogene or query
-- and falls back to the fill's own colour, or mid grey when there is no family.
`_gene_style` used to return early for RNA and pseudogenes, which discarded the
family colour entirely. Stroke width doubles for RNA genes, pseudogenes and
queries, but it only emphasises what the outline colour already says -- it is not
the sole carrier of anything.

### What counts as a family

`family_shared` decides whether a family is worth colouring, and it counts
**occurrences across neighbourhoods**, not distinct accessions. RefSeq gives
identical proteins one `WP_` accession, so a gene conserved in every genome
collapses to a one-member family; counting accessions called that a singleton and
greyed out precisely the genes a conservation tool exists to show. Families are
then numbered in descending order of occurrences, so family 1 is the most
widespread and `_outdesc.txt` reads in order. The renderer colours whatever the
numbering labelled rather than applying its own size rule, so colour and number
can no longer disagree -- including the trailing grey-out pass, which must key on
"was it labelled", not on family size, or it silently undoes the colour a
conserved one-accession family was just given.

### SVG 1.1, not CSS

Figures are opened in Illustrator and Inkscape as often as in a browser, so the
output stays inside SVG 1.1. `rgba()` is CSS colour syntax: browsers accept it,
Illustrator does not and falls back to opaque black, which turned every gene
outline into a filled black arrow and hid the whole domain figure. Transparent
fills are `fill="none"`, and translucency uses the `fill-opacity` attribute, which
is part of SVG 1.1. Keep new colour values as hex or a colour keyword.

### Where gene styling lives

`_gene_style`, `_accent`, `_stroke_width`, the pastel helpers and the `classic`
property all sit on `_FlaGsBase`, not on `OperonView`. Both renderers draw genes,
so a helper added to one of them leaves the other raising `AttributeError` at
render time -- which happened four times, each caught only when a figure silently
failed to draw. Anything about how a gene looks belongs on the base class.

### Pseudogenes and annotation styles

A pseudogene is flagged by `gene_biotype=pseudogene` on the gene row, but whether
it also has a CDS child varies by annotation: Salmonella's have none, so the gene
kept a null accession and picked up the `pseudogene*` marker, while Bradyrhizobium's
do, and the CDS overwrote the accession with its own id. The marker was then absent
and the renderer had nothing to key the navy border on, so 381 pseudogenes in one
genome drew as ordinary genes. The biotype now survives the CDS, and `pseudo=true`
on a CDS is honoured too.

Classic mode fills unclustered genes white rather than grey, matching the original
FlaGs. The translation happens in `_gene_style` after the family lookup, because
`_family_colors` assigns grey to every unlabelled accession, so a default only
reached genes absent from that map.

### Overlays compose

`OperonView` draws family colour, family numbers, domain wedges, protein features
and secretion bands as independent layers rather than choosing one. It used to
pick a single exclusive `mode`, so `sismis` silently suppressed domains and
colouring, and `domains` suppressed family numbers. The gene fill is an opaque pastel tint of the family colour when wedges sit on
top: opaque so a secretion band underneath is covered rather than showing
through, and light so it does not compete with the saturated wedges. The outline
takes the same tint, except on a gene whose outline carries an accent -- an RNA,
pseudogene or query marker lightened to match would stop marking anything. Domain
wedges are drawn N-to-C along the gene as displayed; a wedge is thin at its start
and tall at its end, so each new start hides only a sliver of the previous end. Legends stack: one
panel per overlay actually present.

The standard figure set lives in `visualisation_table.tsv` beside the code, not in
a string inside it, so a site can define its own standard by editing that file.

### Figures are table-driven, and separate

Every figure comes from a row of `visualisation_table.tsv`: `FigureSpec` describes
one figure, `render_figure` draws it, and both the main run and `flags_redraw.py`
go through `render_all`, so a redraw reproduces what the run produced rather than
approximating it.

Nothing outside `flags_view.py` and `flags_redraw.py` draws anything. The pipeline
writes tables and then runs `flags_redraw.py` as a subprocess, so there is exactly
one rendering path and a run's own figures are produced by the same code a later
redraw uses. `flags_tree.py` builds trees and returns Newick; the drawing that used
to live there is now `NeighborhoodVisualizer` in `flags_view.py`. `family_numbers`
moved into `FlaGs3.py` because family labels are a property of the tables; the
renderers read them back from the `family` column rather than numbering a second
time, which is what used to let the two drift apart.

`_draw_tree` sits on `_FlaGsBase`, so any mode can allocate a tree gutter and draw
into it -- previously only `triangles` could, which is why `tree_width` did nothing
elsewhere.

Only two files are involved, and the split between them is a real seam rather than
a layer for its own sake: `flags_redraw.py` knows about files, tables and the CLI
and nothing about SVG; `flags_view.py` takes objects and returns SVG text and knows
nothing about where anything is stored. The table spec and the CLI were separate
modules until the pipeline started shelling out, at which point the CLI was the
only caller and the split bought nothing.

`FlaGs3.py` imports nothing from either. It does not even name the table: if no
table is given, `flags_redraw.py` writes the default one itself. The pipeline's
entire contact with visualisation is a subprocess call and the list of `.svg`
lines it prints.

Three outputs exist purely so a redraw is possible: `_features.tsv` (tmhmm/signalp
regions used to be rendered and discarded), the `is_rna` column in `_operon.tsv`,
and the normalised `contig/start/end/type` columns leading `_secretion.tsv` --
Sismis' own column names vary, so the coordinates were otherwise unfindable.

`features_allowed` carries a `domains` token that was not in the original sketch of
the table. Without it `neighbors` and `domains` differ only by `family_numbers`,
and nothing says which figure gets the domain overlay.

### BlastP entry point

BlastP queries arrive from two places -- `--blast_input` and lines in the main
list whose second column is `BLAST` -- so `resolve_blast` takes a list of queries
and one `BlastSearcher` serves all of them. `parse_query` works on lines rather
than a path, which is what lets the same accession-or-sequence rules apply to a
file and to an inline entry. A query with no hits warns and the run continues;
only every query failing is fatal.

`flags_blast.py` turns one accession or sequence into a query list, mirroring
webFlaGs' two input boxes. `read_query` decides accession-vs-sequence by pattern
rather than by asking the user which they supplied; anything that is neither
raises rather than being silently BLASTed as nonsense.

Both modes share one path because an accession is resolved to residues through
Entrez first. QBLAST would accept an accession directly, but resolving first means
`--blast_mode local` works with accessions too, and the same bytes are sent either
way.

Database names differ between the web service (`refseq_select_prot`) and the
downloadable sets (`refseq_select_protein`), so `DB_ALIASES` maps stable FlaGs3
names per mode and passes anything unrecognised through unchanged — which is also
how a local database path is given.

Hits are appended to the `-i` queries, deduped on the version-stripped accession so
`NP_414542.2` from BlastP does not duplicate a `NP_414542.1` the user paired with a
specific assembly. Order is best-hit-first and preserved.

### Identifying to NCBI

`flags3` is the registered tool name and is sent on every NCBI request:
`Entrez.tool` for the E-utilities in `FlaGs3.py` and `flags_blast.py`, the `tool`
parameter on QBLAST, and a `User-Agent` on the download session. NCBI uses this
plus the address from `-u` to contact a tool's author before blocking it, so a
new code path that reaches NCBI should set it too.

### MGnify and API v2

`MgnifyGenomeDownloader` uses MGnify's APIv2 (`/metagenomics/api/v2/genomes/{acc}`).
v1 was deprecated in June 2026, switched off from September 2026, and served frozen
data before that, so it is not a fallback worth keeping.

The response parsing is deliberately loose. `_url_map` keys files on the **basename
of the download URL** rather than on an `alias`/`id` metadata field, and
`_download_url` accepts a URL under any of several key names. This survives v2
renaming its metadata fields and works unchanged whether a download URL points at
the API or at the FTP server. If the detail endpoint carries no downloads, the
`/downloads` sub-resource is tried once.

Compressed files are preferred: for each slot the `.gz` name is tried before the
plain one, and the local filename keeps whatever extension was served. That needs
no parser change because `NeighborhoodExtractor._open` already dispatches on the
`.gz` extension. It matters because the v1 API served these files uncompressed
(its own metadata said `"compression": false`) which is why MGnify downloads ran
an order of magnitude slower than NCBI's gzipped equivalents.

One diagnostic bug is worth remembering: `_stream` used to `return False` on any
non-200 without recording anything, so HTTP 429 throttling was invisible and
looked like a slow network. Failures now always land in `self.failures` with the
status code and any `Retry-After` value.

---

## Rendering and reporting

`OperonView` draws all three figures — neighbourhoods, domains, secretion — from
one code path with a `mode` switch, so layout stays identical between them.
`flags_tree.NeighborhoodVisualizer` is separate because it must align rows to
tree leaves.

Text width is measured from a per-character Arial metrics table rather than
estimated, and the font stack is pinned to metrically identical faces (Arial,
Liberation Sans, Helvetica). A generic `sans-serif` would resolve differently per
platform and render text at a width the layout never reserved.

`ReportWriter` owns every TSV, text and FASTA output. It takes the leaf order under
`--tree_order` and re-keys its per-row map, so the tables and the figures agree
on row order rather than only the figures being sorted. The FASTA outputs come
from the extractor's sequence tables; queries are excluded from
`_flankgene.fasta` because `extractor.sequences` holds every gene in the window,
query included.

---

## Concurrency

| Where | Model | Bound by |
|---|---|---|
| Downloads | thread pool, ≤10 assemblies, files within each in parallel | network, rate limiter |
| Clustering | thread pool over queries, `cpus=1` each | CPU |
| Domain scan | pyhmmer internal threads | CPU |
| tmhmm / signalp / sismis | 3-thread pool, run together | cloud / subprocess |

The last group runs concurrently because each spends its time waiting on a remote
service or a subprocess rather than on local CPU. Each measures its own elapsed
time, so the `-vb` timing table reports true per-tool cost even though the three
overlap in wall time. `TOTAL` is measured wall clock, not the sum of stages —
the stages overlap, so summing them would over-count.

### External commands live in a table

The scripts split by what they actually do: `pfamA_loader.sh` and
`defenceFinder_loader.sh` fetch data, while `signalp_installer.sh` and
`deeptmhmm_installer.sh` build software environments.

`signalp_installer.sh` and `deeptmhmm_installer.sh` build the fixed environments
`flags3-signalp` and `flags3-deeptmhmm` and then rewrite their own row of
`tools_table.tsv` with the interpreter and directory they produced. Neither can
download its tool: both are licensed, so the scripts take a package the user has
already obtained. Fixing the environment names is the point -- the table row then
means the same thing on every machine.

Both installers pin more than the tool's own instructions ask for, because the
instructions predate their dependencies moving on: SignalP needs `numpy<2` to go
with `torch<2.0`, since torch 1.x is compiled against the NumPy 1.x ABI and fails
at import against NumPy 2.

DeepTMHMM's `requirements.txt` pins a `+cu92` torch build that is not on PyPI, so
the installer resolves torch itself and strips the torch lines before installing
the rest. It installs the CPU build rather than the pinned one: CUDA 9.2 does not
support current GPU architectures, and torch then dies in cuBLAS on the first
matmul. The test run also clears `CUDA_VISIBLE_DEVICES`, so a CUDA build that is
already installed still verifies on CPU. Doing it the other way round -- letting pip read the file as written --
fails on any machine without that CUDA index.

`flags_tools.py` holds one row per external program, defaulting to the built-in
invocation and overridable in `tools_table.tsv`. Each row carries a command
template and a working directory; a relative program name is resolved inside that
directory. The point is not cosmetic: DeepTMHMM pins Python 3.8 and SignalP 6 pins
PyTorch below 2.0, so neither can run in the FlaGs3 environment, and naming the
interpreter in the command is the only way to reach them. `-lth` and `-lsp`
therefore take no arguments -- where a tool lives does not change between runs --
and each implies its own feature, so `-lsp` alone is enough and `-sp -lsp` is not
a thing anyone has to type.

### Running the feature tools locally

`_LocalScanner` passes an output path that does not yet exist and lets the tool
create it: DeepTMHMM's `predict.py` refuses to run if its `--output-dir` is already
there, so pre-creating it broke every local run. Tools that expect to create their
own output directory are the norm; one that requires an existing directory would
need the scanner to make it.

`_LocalScanner` shells out and hands the output to the same parsers the cloud path
uses, so a local and a cloud run produce identical `_features.tsv`. The two tools
are switched independently because they install differently: SignalP 6 has a real
`signalp6` CLI, whereas DeepTMHMM has no pip entry point and local use means a
licensed checkout with `predict.py`, run from its own directory. pybiolib once
offered `machine="local"`; current versions raise on it, so that route is closed.

`available()` is checked before running, and an unusable tool skips that feature
rather than falling back to the cloud -- a run asked to stay local should not
start uploading sequences because of a typo. The reason is appended to
`_runinfo.txt`.

Batching applies to the cloud path only -- it exists for the submission limit, not
for the tools -- so a local run sends everything in one call.

### BioLib submission constraints

Three things in `flags_features.py` look removable and are not:

- **`_SUBMIT_LOCK` around `load()` + `cli()`.** pybiolib's `attempt_sign_in()` is
  a check-then-act on a shared singleton, and its `UserState` lock file is
  created with `fail_fast_on_lock_acquire=True`, so a second concurrent caller
  gets no retry. Worse, the loser sets `_is_in_memory_only` on the shared object,
  which makes the winner's `__exit__` skip releasing the lock file — leaving a
  stale lock that degrades *every later run on that machine*. `warm_up()` burns
  the one-shot sign-in on the main thread before the pool starts; the lock covers
  the rest. `job.wait()` stays outside it so the cloud jobs still overlap.
- **The `chdir` into the scratch directory.** The fasta argument must reach
  `app.cli()` as a bare relative name. Relative names are mounted at
  `/query.fasta` and passed through unchanged; an absolute path is mounted under
  a hashed directory and rewritten *without* its leading slash, so the app
  receives `hash/query.fasta` and never finds its input. `chdir` is process-wide,
  which is safe only because it is held inside `_SUBMIT_LOCK` and because
  `main()` absolutises `args.output` and `args.temporary` before any thread
  starts — the concurrent sismis task holds paths under both.
- **SignalP's `--output_dir output`.** The app's own `generate_output.py` reads
  `output/output.json`. Any other directory name makes the remote step fail.

---

## Extending it

**Adding a genome source.** Subclass `_GenomeDownloader`, implement `_fetch_one`
returning `GenomeFiles`, define `SUFFIX`, and add a routing test alongside
`is_mgnify_accession` in `main()`. Nothing downstream needs to change. If the
source uses an unusual GFF dialect, extend `_cds_accession` rather than branching
on source anywhere else.

**Adding a per-protein annotation layer.** Follow `flags_features.py`: return
`{accession: [(kind, start, end)]}`, add a lazy import and a `_run_*` function in
`main()`, and register it in the parallel task group. If it goes through BioLib,
subclass `_BioLibScanner` rather than calling `biolib` directly, so it inherits
`_SUBMIT_LOCK` and the scratch-directory handling described above.