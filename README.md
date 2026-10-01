# FlaGs3

Predicting protein functional association by analysis of conservation of
genomic context (Flanking Genes).

Give it a list of proteins. It finds each one in its genome, takes the genes
around it, groups the neighbours into families by sequence similarity, and
draws the neighbourhoods so conserved gene arrangements are visible at a
glance. On top of that it can build a tree of the queries, scan the
neighbours for domains, transmembrane regions and signal peptides, and scan
the surrounding sequence for secretion systems, proviruses, plasmids and
anti-phage defence systems, each drawn as its own layer on the same figure.

`User_Guide.md` says how to use it. `docs/Architecture.md` says how it is
built.

## Install

Python 3.11 or later.

```
git clone https://github.com/GCA-VH-lab/FlaGs3.git
pipx install ./FlaGs3
flags3 install core pfam
```

`pipx` (`pacman -S python-pipx`, `apt install pipx`, `brew install pipx`)
puts a `flags3` command on your PATH in its own environment, so nothing
needs activating afterwards. Without pipx, a virtual environment works the
same way except that it must be activated in every new shell:

```
cd FlaGs3 && python -m venv .venv && source .venv/bin/activate && pip install .
```

Conda users can instead do `conda env create -f environment.yml`, which
brings FlaGs3 and the core tools into one environment.

`flags3 install core` sets up mafft, trimal, VeryFastTree, IQ-TREE and
BLAST+; `pfam` downloads Pfam-A for `--domains`. Everything lands under
`~/.flags3/` and nothing else on the machine is touched. `flags3 install`
alone lists the other components (MMseqs2, geNomad, DefenseFinder, PadLoc,
Sismis, DeepTMHMM2, SignalP) and their state.

## Run

```
flags3 run -i proteins.txt -u you@example.org -o myrun
```

`proteins.txt` holds one protein accession per line, optionally with a
tab and an assembly accession to pin the genome. Genomes are downloaded once
into `~/.flags3/genomes/` and reused by every later run.

Common additions:

```
-g 6            six flanking genes each side (default 4)
-r 5000         or every gene within 5 kb of the query
-cr             cluster RNA genes too
-t              tree of the queries, -to to order the rows by it
-d              Pfam domains on the flanking proteins
-ss -gn -df     secretion systems, mobile elements, defence systems
-pdf            PDFs beside the SVGs
```

## What you get

```
myrun/
    report/         neighbourhoods.tsv, queries.tsv, families.tsv, FASTAs,
                    run_summary.txt, and legacy/ with the FlaGs2-era files
    figures/        neighbors.svg, tree.svg, domains.svg, sismis.svg, ...
    <stage>/        one directory per analysis step, with its tables
    config.tsv      every option the run used
    console.log     everything printed
```

Any step can be rerun on a finished directory without repeating the rest:

```
flags3 figures myrun -f my_figures.tsv     # redraw with another figure table
flags3 domains myrun -db defensefinder=~/.flags3/db/defensefinder/profiles
flags3 extract myrun -g 8                  # then flags3 cluster myrun, ...
```

## Coming from 2.3.0

Options and outputs moved; `CHANGELOG.md` has the table. The short
version: `flags3 run` replaces `python FlaGs3.py`, genomes are kept in a
cache instead of `-tmp`/`-k`/`-ul`, `flags3 install` replaces `build.sh`,
and the old tables are in `report/legacy/` under their old names.

## Citing

FlaGs3 is the successor of FlaGs (Saha CK, Sanches Pires R, Brolin H,
Delannoy M, Atkinson GC. *Bioinformatics* 2021, 37(9):1312–1314,
doi:10.1093/bioinformatics/btaa788). Please cite that paper until the FlaGs3
paper is out.

## Licence

GPL-3.0, as FlaGs.
