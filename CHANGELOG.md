# Changelog

## 3.0.0

A rewrite. Same analysis, same figures, different construction.

### For users

- One command, `flags3`, with subcommands: `run`, one per stage
  (`fetch`, `extract`, `cluster`, `tree`, `domains`, `features`, `sismis`,
  `genomad`, `defence`, `report`, `figures`), `install`.
- Every stage can be rerun on a finished run directory with new options;
  `--no-<option>` switches an option off for the rerun.
- Genomes live in a cache (`~/.flags3/genomes/`) shared by all runs.
  `-gd` points elsewhere, `--offline` never touches the network, `-lf`
  looks bare queries up in the directory before NCBI, `--no_cache` keeps a
  run's genomes with the run. `-tmp`, `-k` and `-ul` are gone.
- `flags3 install` replaces `build.sh` and the installer scripts: every
  tool and database goes under `~/.flags3/`, binaries already on PATH are
  used as they are, micromamba is fetched only when something must be
  built. DeepTMHMM is now DeepTMHMM2, local and licence-free.
- DeepTMHMM and SignalP on BioLib run as one job each, submitted together
  and concurrently with the local stages.
- Clustering settings live in the tools table; `-n` and `-ce` override them
  per run. The jackhmmer inclusion thresholds are back at FlaGs2's values
  (`incE` and `incdomE` 1e-10); FlaGs3 up to 2.3.0 used 1e-3 and HMMER's
  default `incdomE`, which clustered more inclusively.
- Output layout: one directory per stage, `report/` for the human tables
  and `report/legacy/` for the 2.3.0 files (`_operon.tsv` and friends,
  byte-compatible), `figures/` for the SVGs.
- Figures: one layer per stage in the SVG, hatched bands drawn behind the
  genes with codes at the row's right, palettes `bright`, `pastel`,
  `classic`, `colourblind`, `monochrome`, `colours.tsv` for overrides.

### Since 3.0.0rc1

- GFF attribute values are URL-decoded (`%2C` → `,`), and a row's species
  falls back to the organism named by the rest of that assembly's proteins.
- `fetch/queries.tsv` has a `via` column (`ipg`, `paired`, `local`, `blast`);
  IPG aliases are identical sequences by NCBI's definition and are not
  aligned locally.
- `run_summary.txt` counts the exported protein sequences and says why that
  differs from the gene count (RNA genes, pseudogenes, repeated proteins).
- A gene running past its contig's declared length widens the contig with a
  warning instead of stopping the run.
- Older run directories stay readable: new columns with a default
  (`contig_edge`, `via`) are filled when absent.

- DefenseFinder and PadLoc are called on the `-sr` span's genes when `-sr`
  is given, not only on the drawn neighbourhood.

- Files in `figures/` and `report/` carry the run's start time in their
  names; `-nt` runs keep plain names.
- `families_legend.txt` is `protein_clusters_legend.txt` (plus
  `rna_clusters_legend.txt`), members ordered by subcluster.
- In a query cluster, `Q` marks proteins that were in the input list and
  `q` the other members.

### Migration from 2.3.0

| 2.3.0 | 3.0 |
|---|---|
| `python FlaGs3.py -i list -u mail …` | `flags3 run -i list -u mail …` |
| `-tmp DIR -k` | genomes are always kept, in `~/.flags3/genomes/` or `-gd DIR` |
| `-ul DIR` | `-gd DIR --offline` (or `-lf` for a mixed list) |
| `-e`, `-n` | `options` column of the clustering row in the tools table |
| `-db` per-tool `scan_range` column | `-sr`/`-sm` apply to every scanner |
| `build.sh`, `*_installer.sh`, `*_loader.sh` | `flags3 install <component>` |
| `<prefix>_operon.tsv` etc. in the output directory | `report/legacy/<run>_operon.tsv` etc. |
| `flags_redraw.py` | `flags3 figures <run>` |
| `-bl` BLAST expansion | `-bi FILE` or a `BLAST`-marked line, unchanged |

### Not carried over

- The per-tool `scan_range` column of the tools table.
- The `_runinfo.txt` timing table (now `report/run_summary.txt`).
- PDF conversion through svglib, rsvg-convert or Inkscape (cairosvg only).

## 2.3.0

Last script-based version. See the `v2.3.0` tag.
