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
window around the query: `±flank` genes by default, or every gene within
`±range_bp` bases under `-r`.

### Two ways to size a window

`_window` returns a pair of indices either way, so everything downstream is
unchanged -- `offset` stays **ordinal** (the n-th gene from the query, sign
flipped on the minus strand), not a distance. That is what lets `-r` cost the
renderers nothing: `classic` and `triangles` lay genes out on an offset grid and
simply get wider, while `versatile` was already drawing to genomic scale.

Range mode walks outwards from the query index rather than bisecting, because the
window is contiguous around a known point. Walking left cannot stop at the first
gene whose end falls short of the boundary: a gene nested inside a longer
neighbour would cut the window early. `_reach` is a per-contig prefix maximum of
`end`, so the walk stops only when *no* earlier gene on that contig can reach the
boundary. Walking right needs no such array, since genes are sorted by `start`.
Both arrays are built once per assembly in `_genes` and cached with it; range
extraction measured the same as gene extraction on real genomes.

A gene straddling the boundary is included, so `up_reached` routinely exceeds the
requested distance by part of one gene. That is the point: a gene half inside the
window is context, not noise.

### The range report

Contig lengths come from the `##sequence-region` pragma, read in the same pass
that parses the features -- the parser used to drop every `#` line. `region`
features are a second source, and the longest gene end on the contig is the
fallback for annotations carrying neither. The fallback under-reports, so a
Prokka genome with no pragma can show truncation that is not there; it reports
what the file actually says rather than guessing.

`up`/`down` in `_rangeReport.tsv` are in **query orientation**, flipped with the
query strand exactly as `offset` is. Genomic left/right would make "20000 up and
5000 down" mean nothing biologically, which is the whole reason the report exists.

Every row is written, not only the truncated ones, because the same table answers
"how many genes did this row actually get" -- the question that comes up as soon
as two rows in a figure are different widths.

### Scanning span is separate from neighbourhood size

`-sr` exists because clustering cost is quadratic in the number of genes drawn
while the scanning tools take a genomic interval and do not care how many genes
are in it. On real data, ±50 kb holds about 85 genes against 9 at `-g 4`, which
is roughly 80x the clustering time; `-g 5 -sr 50000` clusters 11 genes and still
hands the tools 100 kb. Without the split, wide tool context and affordable
clustering are mutually exclusive.

`RangeInfo.scan_start`/`scan_end` are computed at extraction time and clipped to
the contig, so the tools receive an interval that is already known to exist. When
`-sr` is absent the scan window is the neighbourhood's own genomic span.

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

### Clustering memory scales with homology, not with size

A subfamily drawn from few genera has flanking proteins that nearly all hit each
other, so the adjacency approaches a complete graph: 8397 sequences at 90%
density is 63 million edges. Two things then made that far more expensive than it
had to be.

`_name` decoded `hit.name` for every hit, so each edge held its own string object
rather than a pointer to the one canonical name. Measured at that scale, 7.9 GB
instead of 4.4. `search_one` now maps the hit's raw bytes back to the key string
the block was built from.

`_connected_components` built a symmetrised copy holding every edge in both
directions, while the original adjacency was still live -- about 4 GB more at the
same scale. Union-find gives the same components from one parent entry per node,
which measured 3800x smaller on a dense graph and produced identical families.

Together these took a dense 8397-sequence run from roughly 17 GB to roughly 4 GB.
The symptom was the OOM killer on a subfamily spanning 11 genera, while a
subfamily six times larger but spanning 503 genera had run comfortably.

### Components join whichever direction the hit came from

jackhmmer adjacency is genuinely one-directional in a few percent of pairs: A's
profile finds B while B's does not find A, which is normal when one of them sits
in a larger or more diverse family. On a real 950-protein run, 235 hit pairs were
one-directional.

`_connected_components` used to traverse that raw map. Skipping a node only as a
*starting* point is not enough -- a later node reaches back into a component that
has already closed and pulls its members in a second time. 26 of those 950
proteins ended up in two families, so `_clusters.tsv` double-counted them and
`family_numbers` picked one assignment arbitrarily, letting a figure and a table
disagree about the same gene.

The map is now symmetrised before traversal, which is what the "symmetric by
construction" claim always assumed. Components become a genuine partition and the
result no longer depends on iteration order.

### Collapsing before clustering

Clustering is quadratic: each of N queries is searched against a block of all N,
so cost fits `N * (0.137 + 4.4e-5 * N)` seconds per core, measured. `flags_collapse`
runs MMseqs2 first, clusters only the representatives, and gives every member its
representative's family.

It is **opt-in and never automatic**. A size threshold that silently switched
algorithms would make two runs disagree with nothing in the output explaining why,
and the whole point of `_collapse.tsv` is that every propagated family assignment
is traceable: member to representative there, representative to evidence in
`_jackhits.tsv`.

Two things about the ratio are worth knowing before reaching for it. First,
`extractor.sequences` is keyed by accession, so RefSeq's own non-redundancy
already collapses proteins identical across genomes -- MMseqs2 only recovers the
band between the identity cutoff and 100%. Second, the gain is smaller than the
quadratic curve suggests: measured collapse on real flanking proteins was 1.00x
at 90% identity and 1.17x at 50%, because those genomes were different species.
The neighbourhood size (`-g` versus `-r`) remains by far the larger lever.

Collapsing can **split** a family as well as merge one. A sequence that would have
bridged two groups is never searched once it is folded into a representative, so
the bridge is lost. At 90%/80% on real data the families came out identical to an
uncollapsed run; at 50% they did not -- 533 families against 528, with 71 of 950
proteins landing somewhere different. That is the reason the default is 0.9,0.8
and the reason the flag takes the numbers explicitly.

A missing MMseqs2 is a hard error rather than the usual degrade-and-continue,
which is a deliberate exception. Every other optional tool adds a layer; this one
changes what the core stage costs, and quietly falling back to an uncollapsed run
means a job the user expected to take an hour runs for days instead. The error
prints the estimate so the choice is informed.

---

## Tree building

`TreeBuilder` runs three stages: MAFFT `--auto`, a gap-threshold column trim, then
an inference engine. The trim keeps columns where at least `gap_threshold` of the
sequences carry a residue — the same rule as `trimal -gt`, which is what ete3's
`trimal01` ran in the old pipeline. It shells out to trimal when that is on
PATH, since the heuristic modes (`-automated1`, `-gappyout`) have no short
equivalent. When it is not, `_trim` stands in for the default `-gt`, which is a
deterministic column filter worth ten lines rather than a hard dependency; the
heuristic modes fall back to the untrimmed alignment and say so. If trimming
would empty the alignment, the untrimmed one is used instead.

`_trim` was described here as the implementation while nothing called it, and it
would have raised `NameError` if anything had: its first parameter was named
`cls` with no `@classmethod` decorator, and the body referenced `cls.GAP_CHARS`.
Wiring it exposed both.

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

### The console transcript

`flags_log` tees `sys.stdout` and `sys.stderr` from the first line of `main()`,
before `parse_args()`, and buffers what it captures in memory until `attach()` is
given a path. The buffer exists because the output directory's name is not known
until the timestamp is applied and the path absolutised, by which point several
validation messages have already been printed. Capturing from process start and
flushing later is the only way those end up in the file; installing the tee after
the path is resolved would silently lose them.

Lines are split on newlines and tagged per stream, with an unterminated tail held
back until the rest of its line arrives -- `print` writes the text and the newline
as two calls, so tagging each `write()` would break every line in half. Stderr
lines get a `[stderr] ` prefix rather than a second file: one file keeps the
interleaving, which is what makes a `--debug` run readable, and the prefix is
enough to recover either stream with `grep`.

The tee is installed by an explicit call rather than on import, since
`flags_redraw.py` imports `flags_log` too and a module that replaces `sys.stdout`
merely by being imported is a trap.

`record_command` exists because every external tool is run with `capture_output=True`,
so a tool's own diagnostics never reach the terminal and were previously discarded
entirely -- the run reported "produced nothing" with no way to find out why. Tools
whose stdout *is* the payload (mafft's alignment, VeryFastTree's newick, blastp's
hit table) record only stderr. `TreeBuilder._run` records on `CalledProcessError`
and `OSError` as well as on success, because a command that failed or was missing
is the one worth having in the log.

`atexit` closes the transcript, and it is registered first so it runs last: the
lock release and any shutdown message are captured before the file closes. An
uncaught `SystemExit` prints its message through the still-installed tee, so the
`sys.exit("Error: ...")` paths and a Ctrl-C both leave a complete log.

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

### Defence systems run on a synthetic GFF

DefenseFinder and PadLoc call *systems*, which means they need gene order and
gene adjacency, not just sequences. They are gene-based tools, so they take the
neighbourhood from `-g`/`-r` and never look at `-sr`.

Each neighbourhood is written as its own synthetic replicon, `r0`, `r1` and so
on. That boundary is the point: two queries on the same contig but 400 kb apart
would otherwise sit next to each other in the file and invite a system call
across the gap that does not exist. One replicon per neighbourhood makes gene
adjacency mean what it should.

Genes keep their **real coordinates**, and `##sequence-region` runs from the
first gene's start rather than 1, so gene spacing and contig-edge tests stay
truthful. The consequence is that nothing has to be mapped back: a system's band
is the span of the genes it named, looked up by tag. Unlike the Sismis and
geNomad paths there is no offset table here at all, because there is no slicing.

Both tools produce a system-to-protein assignment and differ only in how they
spell it, so one reader each is enough. A system both tools call on the same span
is emitted once, with `called_by` naming both -- two identical bands stacked on
one another would just look like a rendering fault.

### Two tool tables, one of them ignored

`tools_table.tsv` is the shipped default and the repository tracks it. Every
installer writes `tools_table.local.tsv`, which `.gitignore` covers, and `load()`
reads the default first and layers the local file over it row by row. A local row
only carries what it changes; blank cells fall back to the default rather than
blanking it.

The point is that an installed path is a property of a machine, not of the
project. Writing them into the tracked table meant every install produced a diff
full of somebody's home directory, and those diffs got committed.

### One jackhmmer call carries many queries

Each `jackhmmer` call leaves roughly 2 MB behind whatever is done with the
result. At one call per query that is invisible on a small run and fatal on a
large one: a subfamily of 8397 flanking proteins makes 8397 calls, which is about
19 GB, and the OOM killer took three attempts at it.

Queries now go in chunks of `CHUNK`, so the same subfamily makes 84 calls rather
than 8397. Measured on a dense block: 2.28 MB per call one at a time against
0.16 MB per query batched, with the wall time unchanged. The thread pool is
untouched, since that is where the parallel speedup comes from and it was never
the problem.

Chunk size is set so the pool has several chunks per worker rather than one
each. One chunk per worker leaves nothing to balance with, and jackhmmer cost per
query varies widely -- a query hitting hundreds of relatives builds a large MSA
while one hitting nothing converges immediately -- so the slowest chunk sets the
pace for the whole stage. Simulated on 300 queries across 26 workers, one chunk
each wastes about 160% against the ideal split; four each wastes about 40%.

Below a few hundred queries this falls to one query per chunk, which is exactly
what small runs did before batching and is the fastest thing for them. The call
count only needs reducing at thousands of queries, which is also the only place
memory is a constraint.

### The window's genes are output too

The scanning tools work on the whole scan window while only the `-g` neighbourhood
is clustered and drawn, so a defence system can name genes that appear nowhere
else in the output. On a test run 91 of 91 genes named as system members were
absent from every other file: the hit said which accessions it spanned and
nothing said what they were.

`_window.tsv` lists every gene inside the window -- coordinates, strand, offset,
product, and whether it was in the drawn neighbourhood -- and `_window.fasta`
carries their proteins. Both come from `window_genes` and `window_sequences`,
which extraction already built, so writing them costs a pass over data already in
memory.

`_window.fasta` headers carry the same facts as key=value pairs -- query,
assembly, contig, protein length, position and strand -- so the file describes
itself without the table beside it. Position is relative to the query gene in
query orientation, negative upstream, matching `offset` and the range report's
up and down columns rather than genomic left and right. A protein sitting in two
windows gets one record per window, since its position differs in each.

`in_neighbourhood` is the column that matters when reading a figure beside the
tables: a `no` gene is real and was scanned, it simply fell outside what was
drawn.

### A paired assembly is not the last word

A protein given with an assembly in the input skipped IPG entirely, so a
withdrawn assembly, a failed download, or an accession that simply is not in that
assembly left the protein unresolved with no second attempt. On a dataset whose
assembly column was compiled at some earlier date that is a steady, silent loss.

`--remap` retries those proteins through IPG once, which is the only route that
can find where the protein actually lives, and downloads and extracts whatever
that finds. The retry runs after extraction rather than after download because
both failure modes -- no usable genome, and a usable genome without the protein --
look the same from there, so one code path covers both.

It is opt-in rather than automatic because it spends NCBI requests and possibly a
second round of downloads on a failure that may be expected. A verbose run always
reports how many proteins it would have applied to, so the cost can be judged
before it is paid.

### IPG resolution is chunked

`_map_ipg` put every unresolved accession into one `efetch`. That is fine for a
demo and wrong at scale: NCBI recommends 200 ids per request, and a report for
thousands comes back slowly, truncated, or not at all. A truncated IPG report is
the dangerous case -- it is not an error, it just looks like those proteins have
no assembly, so the run continues and quietly drops them.

Requests now go in chunks of 200 and a failed chunk costs only its own 200
accessions; the rest still resolve, and the count that failed is printed and
recorded. On this dataset the 1336 rows with no assembly written become 7
requests instead of one.

### Figures too tall to open are written as parts

A subfamily of 8000 rows is a 190,000 px SVG. Browsers open it; Illustrator and
librsvg refuse, so it cannot be edited or converted to PDF.

`render_all` renders once, measures, and splits only if the result is over
`-fh`/`--figure_height`, which defaults to 16383 px. Rather than guessing a row count, `_rows_per_part` renders one probe and
solves the layout: height is a fixed part -- title, legend, axis -- plus a cost
per row, and two measurements give both terms. A guessed 10% margin left parts 1%
over the limit, because the fixed part is repeated by every part and does not
shrink with the row count.

Figures carrying a tree panel are never split, since the tree spans every row and
cutting the rows without pruning the tree would produce a figure that lies. Those
are left whole with a warning saying why.

### Parsed genomes are freed as they are finished

`_gff_cache` and `_faa_cache` are keyed by assembly and were never evicted, so a
run held every genome it had ever touched. Measured at 4.6 MB per assembly, which
is nothing for a demo and 32 GB for a subfamily of 6913.

Extraction is planned in input order but performed grouped by assembly, so each
GFF and protein table is still parsed exactly once and `forget()` can drop it the
moment that assembly is done. Retained memory falls 226x, to 0.14 GB for the same
subfamily. `_contig_len` survives eviction because the scanning tools and the
range report read it after extraction is over, and it is small.

Grouping changes the order sequences arrive in, which used to change the output:
family numbers came out of a dict whose insertion order followed extraction
order. Components are now ordered by `(-size, first member)` and clustering
iterates sorted names, so the same input gives the same family numbers however it
was ordered. FASTA output is sorted for the same reason. The row-ordered tables
still follow the input, which is what they are for.

### Attributes are parsed once per line

`_attr` built and searched a regex for every attribute lookup, and every GFF line
needs three to six -- 609,534 regex searches for 26 genomes. Splitting the
attribute column into a dict once per line is 10x faster on that path and halves
the time in `_genes`. First occurrence wins, as the regex did, since GFF allows
repeated keys.

### A tool's cost can be per invocation, not per base

Measured on two genomes: geNomad took 1109 s on whole genomes and 918 s on 50 kb
windows. Cutting 74% of the sequence bought 17% of the time, because the run is
dominated by loading a 1.6 GB database and a neural model. At 459 s per
invocation, one call per assembly is 37 days for a subfamily of 6913.

So the batching, not the windowing, is what makes it possible. `write_batches`
puts every assembly's windows into as few FASTA files as it can, each record
carrying an entry in an offset table that names the assembly it came from, so one
invocation covers thousands of genomes and hits still map back to the right
genome and coordinates. `MAX_BATCH_BASES` caps a file at 200 Mb so a subfamily
becomes a handful of calls rather than one unbounded one.

Windowing still matters -- it cuts what the tool reads and it is what makes a
batch file a sensible size -- but on its own it was never going to be enough.

Sismis gets the same treatment for the same reason.

### No second way to scan

When Sismis and geNomad moved to one invocation over many genomes,
`SismisScanner.scan_windows` and `scan_assembly` -- the per-assembly path they
replaced -- stayed in the file, along with `write_windows`, `place`,
`scanned_bases` and `_decompressed`, which nothing else called. Two ways to do
the same thing, one of them unreachable and none of it exercised.

They are gone. A reachability pass over the whole tree now finds no unreferenced
definition that was not already in 1.4.0.

### The cut windows describe themselves

The records in `windows/<span>/batch*.fna` are what Sismis and geNomad actually
read, and they were named `w0`, `w1` and nothing else. The id has to stay a short
token, because the tools echo it back and the offset table is keyed on it, so
everything a reader needs goes in the description after it: the assembly, the
contig, the query the window was cut for, the length, the genomic slice, the
analysed span without the margin, and the position relative to the query gene.

`place_batched` falls back to the first whitespace-separated token, so a tool that
reports the whole header rather than the id still resolves.

`scan_windows` used to recompute the span rather than taking the one the range
report recorded. The report clips to the contig; the recomputation did not, so at
a short contig a window claimed an analysed range running past the end of the
sequence, and a hit near that end was called `full` when the report said the range
had been truncated. It now uses the recorded span.

### One copy of the windows

Sismis and geNomad both take nucleotide FASTA, and cutting the same windows twice
costs hundreds of gigabytes across a full batch. `shared_batches` keys the cut
files by the span that produced them, so tools asking for the same span reuse one
set. The key is the span rather than the tool, because a `scan_range` override in
`tools_table.tsv` genuinely does mean different windows and those must not be
shared. A lock guards the store, since the scanning tools run on separate threads.

The files are kept rather than cleaned up: they are the exact input the tools saw,
which is what makes a call reproducible after the fact.

### Defence tools read the window, not the drawing

`-g` decides what is clustered and drawn. It should not decide what a system
caller sees: a defence system spanning genes just outside an 11-gene neighbourhood
would be truncated into a fragment or missed. `window_genes` holds every gene
within `-sr` of the query, collected in the same pass as the drawn neighbourhood
and costing nothing extra since the GFF and FASTA are already parsed and cached.
DefenseFinder and PadLoc get those; clustering and the figures still get the
drawn set.

On real data this is 101 genes against 11. Systems found outside the drawn
neighbourhood appear in `_defence.tsv` and are clipped from the figure, which is
the right way round -- the table is the result, the figure is a view of part of it.

### Overlapping bands need lanes

Every band on a row was drawn at the same `y` with the same height. Two tools
calling the same locus produce two hits that `merge_calls` will not merge when
their spans differ by a gene, and the narrower one was then drawn inside the
wider one and invisible. On real data that reads as "DefenseFinder is in the
legend but has no band", which is exactly what it looked like.

`_band_lanes` does greedy interval partitioning per row: bands are placed in the
first lane whose last band ended before this one starts, so two bands share a
lane only when they do not overlap. The row's band height is divided by the lane
count, so a row with three overlapping systems shows three thin bars rather than
one bar and two lies. Rows are independent, and a row with no overlap looks
exactly as it did.

### Bands are numbered when names will not fit

A row can carry several defence systems and `CBASS_Type_II` beside `RM_Type_I`
beside `Gabija` does not fit next to a row. `BAND_CODES` maps a hit's source to a
prefix -- `defence` to `D` -- and types from those sources are numbered `D1`,
`D2` in the legend order, stamped on the band itself when it is wide enough for
the text, and used for the labels beside the row. Sources without a prefix keep
their names, so Sismis' `T6SS` and geNomad's `Caudoviricetes` are unchanged;
those rarely carry more than one or two per row. Adding a source to `BAND_CODES`
is the whole change if that stops being true.

Numbering is global within a source rather than per tool, so a system both tools
called keeps one number and one colour in both legend panels.

### The legend is split by what made the call

Band hits carry the tool that called them, from the `called_by` column, and
`_band_panels` groups the legend by it. A system both tools agree on appears
under both, with the same number, which is what makes agreement and disagreement
visible at a glance. Sources with no calling tool -- Sismis, geNomad -- fall into
one shared panel, since splitting a single-tool legend by tool would be noise.

### A database is identified by its contents

`genomad download-database DEST` creates `DEST/genomad_db`. The loader defaulted
`DEST` to `./genomad_db`, so the database landed in `genomad_db/genomad_db`, and
then looked for it with `find "${DEST}" -maxdepth 2 -type d -name genomad_db`.
`find` matches the directory it starts from, and depth 0 comes first, so `head -1`
returned the empty wrapper. `tools_table.tsv` was given a path containing nothing
but another directory.

Both halves are fixed, but only the second one matters: the loader now finds the
database by looking for `genomad_marker_metadata.tsv` or `version.txt` and taking
the directory holding it. A name can collide with its own parent; contents cannot.
The default destination no longer ends in `genomad_db` either, so nothing nests
in the first place.

`resolve_database` applies the same test at run time and descends one level into
`genomad_db` when the configured path is a wrapper, printing a note saying so.
Existing installs with the nested layout keep working without anyone editing the
table, and the note means the behaviour is visible rather than magical. A
directory that is neither says what it lacked and where the real one usually sits.

### Finding a binary is not the same as being ready to run

`available()` was written three times, once per tool module, and each copy ended
its absolute-path branch with `return True`. Finding the executable was treated as
the whole answer. For geNomad it is not: the command also needs a database, and
because the installer writes an absolute path, that branch returned early and the
database check below it never ran. geNomad was launched with an empty `{db}` and
answered `Missing argument 'DATABASE'`.

`flags_tools.locate` now owns finding the binary, and each module checks what else
that tool needs afterwards, unconditionally. Splitting the two questions is the
point -- one copy of the resolution logic cannot drift, and a readiness check
cannot be short-circuited by how the tool happened to be installed.

`missing_values` is the backstop. It compares a command's placeholders against the
values being substituted and refuses to run when one the command actually uses is
empty, because `{db}` expanding to nothing does not fail loudly -- it silently
produces a command with one fewer argument, which the tool then misreads.

### Installers must not ask Python where an environment is

The installers used to find their environment with
`conda run -n <env> python -c "import sys; print(sys.prefix)"`. That works only
when the environment contains Python. PadLoc is an R package, so `python`
resolved to whatever was next on `PATH` -- the caller's own interpreter -- and
answered `/usr`. The installer then looked for `/usr/bin/padloc`, did not find
it, and reported the install as broken when it had succeeded.

`env_prefix` asks `conda run -n <env> printenv CONDA_PREFIX` instead. `printenv`
is coreutils, so it is present whatever the environment holds, and `CONDA_PREFIX`
is set by the activation `conda run` performs. `conda env list` parsed by name is
the fallback, and a prefix that does not resolve to a directory is a hard error
rather than a path that silently points somewhere wrong.

Each installer verifies its binary exists **before** using it, and lists the
environment's `bin` when it does not, since "missing after the install" is not
much help on its own. Fetching models or databases happens after that check, so a
mislocated environment fails at the cheap step rather than partway through a
download.

`deeptmhmm_installer.sh` and `signalp_installer.sh` still use the Python form.
Both wrap Python tools, so their environments contain Python by construction and
the assumption holds there.

### A failed scan is not an empty scan

Per-genome statuses used to be counted only by length, so "scanned 2 assemblies,
found 0 mobile elements" was printed whether the tool had run and found nothing
or crashed on both. Those are opposite conclusions and the second one looked like
the first.

`scan_outcome` splits statuses into scanned, failed and skipped. Failure on
everything is a warning naming the first error and where to read the rest;
partial failure says how many. The verbose line only claims what was actually
scanned.

Tool output goes through `flags_tools.brief`, which keeps the **last** few
non-empty lines rather than a character slice. Tools announce themselves before
they fail, so the head of their output is a banner and the tail is the reason --
slicing from the front reported "Executing geNomad annotate (v1.12.0)" as the
error. Statuses also have their whitespace collapsed before being written, since
a status holding a newline splits its own row in a TSV.

### Absolute tool paths need their own bin on PATH

`tools_table.tsv` rows written by the installers point at an absolute path inside
a conda environment, like `.../envs/flags3-genomad/bin/genomad`. Calling that path
runs the right binary but does **not** activate the environment, so the
subprocess inherits the parent's `PATH`.

Self-contained tools do not care. Tools that shell out to siblings do: geNomad
looks up `mmseqs` and `aragorn` by name, finds neither, and stops partway with
"These dependencies are missing" even though both are installed in the very
directory it was launched from.

`flags_tools.env_for` prepends the program's own directory to `PATH` whenever the
command is an absolute path, which is exactly what activation would have done for
sibling lookups. It returns `None` for relative commands and when the directory is
already on `PATH`, so `subprocess.run(env=None)` inherits normally and nothing
changes for tools that were already fine. Every call site that runs a
`tools_table` command passes it.

The alternative was `conda run -n <env>`, which activates properly but adds a
wrapper process and its own output buffering to every call.

### One rule for where a tool looks

Tools split by what they consume, and the split decides which flag controls them:

| consumes | tools | controlled by |
|---|---|---|
| individual genes or proteins | domain scan, DeepTMHMM, SignalP | `-g` / `-r`, i.e. the neighbourhood |
| a stretch of genomic sequence | Sismis, geNomad | `-sr`, whole genome when absent |

There is no third case and no per-tool special pleading in the code. A tool that
wants a different span than `-sr` says puts a number in the `scan_range` column of
`tools_table.tsv`, or `genome` to opt out of windowing entirely -- the same table
that already decides how each tool is invoked.

`flags_scan` holds the window machinery once: merging, slicing, and mapping
coordinates back. `flags_secretion` and `flags_genomad` differ only in what they
shell out to and how they name a hit.

### Scanning windows instead of whole genomes

`--sismis` used to hand each assembly's whole genomic FASTA to Sismis. It now
cuts the contig around each query instead, which on real data scans 2-6% of a
genome rather than all of it.

Windows are **merged before scanning**, not scanned per row. Two queries 30 kb
apart on one contig produce one interval, not two overlapping scans, and every
window for an assembly goes into a single FASTA as separate records, so an
assembly is one subprocess however many queries land on it.

`--scan_margin` exists because a system straddling the window edge would
otherwise be cut in half and called as a fragment or missed. The margin is scanned
but not treated as in-range: a hit lying wholly inside the analysis window is
`full`, one reaching into the margin is `partial`, and both are kept. Merging
happens on the padded intervals, so adjacent windows whose margins touch also
collapse into one.

Coordinates come back relative to the slice, so each record is written as `w{n}`
and mapped through an offset table on the way out. Nothing downstream sees slice
coordinates.

Without `-sr` the tools get whole genomes, which is what they got before windows
existed, so the flag is the whole of the switch.

### geNomad names what it found

geNomad calls proviruses and plasmids, and a band saying "mobile element" wastes
the call. The virus band is labelled with the lowest rank geNomad assigned --
`Caudoviricetes`, not `virus` -- and a plasmid band says `plasmid` or
`conjugative plasmid` depending on whether conjugation genes were found. The
label is the hit's `type`, which is what the legend and the band colouring key
on, so two different taxa get two different colours for free.

Proviruses carry a `coordinates` field relative to their record; a whole-record
call has none and spans the record. Both go through the same offset table as
Sismis' hits, so nothing downstream sees slice coordinates.

geNomad reuses `_secretion.tsv`'s renderer through a separate `_genomad.tsv` and
a `genomad` token in `features_allowed`. Bands are bands; the only reason they are
separate files is that mixing secretion systems and proviruses in one legend
would be unreadable. `all-in-one` takes both.

#### The passthrough columns are prefixed now

Sismis' own columns are appended to `_secretion.tsv` after the normalised ones,
and Sismis emits `start`, `end` and `type` -- the same names. `csv.DictReader`
keeps the **last** value for a duplicated field, so `row["start"]` in
`flags_redraw` was already reading Sismis' column rather than the normalised one.
That was harmless while both held the same number. Windowing broke it: the
normalised column holds a genome coordinate and Sismis' holds a slice-relative
one, so a band would have been drawn at position 201 instead of 1942310. The
passthrough columns are now prefixed `sismis_`.

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
`defensefinder_hmm_loader.sh` fetch data, while `signalp_installer.sh` and
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