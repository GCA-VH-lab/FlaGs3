#!/usr/bin/env bash

set -e

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_TABLE="${HERE}/tools_table.local.tsv"

# Installers write here, never to the committed tools_table.tsv, so machine
# paths stay out of the repository. FlaGs3 reads the default table first and
# lets this one override it row by row.
if [[ ! -f "${TOOLS_TABLE}" ]]; then
    printf '#name\tcommand\tdirectory\tscan_range\n' > "${TOOLS_TABLE}"
fi
DEST="${1:-${HERE}}"

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<'USAGE'
Usage: genomad_loader.sh [destination]

Downloads the geNomad database, about 1.6 GB, and records where it went in the
db column of the genomad row of tools_table.tsv.

geNomad creates a genomad_db directory inside whatever destination it is given,
so the default destination is the directory holding this script and the database
lands in ./genomad_db. Passing ./genomad_db yourself would give you
genomad_db/genomad_db, which still works but reads badly.

Run genomad_installer.sh first; this needs the genomad command.
USAGE
    exit 0
fi

GENOMAD=$(python3 - "${TOOLS_TABLE}" <<'PY'
import sys, os, shlex
path = sys.argv[1]
if os.path.isfile(path):
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        if cols and cols[0].strip().lower() == "genomad" and len(cols) > 1:
            parts = shlex.split(cols[1])
            if parts:
                print(parts[0])
            break
PY
)
[[ -z "${GENOMAD}" ]] && GENOMAD="genomad"

if ! command -v "${GENOMAD}" &> /dev/null && [[ ! -x "${GENOMAD}" ]]; then
    echo "ERROR: ${GENOMAD} not found. Run genomad_installer.sh first."
    exit 1
fi

mkdir -p "${DEST}"
echo "Downloading the geNomad database into ${DEST}"
echo "  this is about 1.6 GB and takes a while"
"${GENOMAD}" download-database "${DEST}"

# Find the database by what is in it, not by what it is called. The destination
# may itself be named genomad_db, in which case matching on the name finds the
# wrapper rather than the database.
MARKER=$(find "${DEST}" -maxdepth 3 -type f \
         \( -name "genomad_marker_metadata.tsv" -o -name "version.txt" \) \
         2>/dev/null | head -1)
if [[ -z "${MARKER}" ]]; then
    echo "ERROR: no geNomad database under ${DEST} after the download."
    echo "       Looked for genomad_marker_metadata.tsv and version.txt. What is there:"
    find "${DEST}" -maxdepth 2 2>/dev/null | head -10 | sed "s/^/         /"
    exit 1
fi
DB=$(dirname "${MARKER}")
echo "  database at ${DB}"
echo "  contains $(ls "${DB}" | wc -l) files"

python3 - "${TOOLS_TABLE}" "${DB}" <<'PY'
import sys, os
path, db = sys.argv[1], sys.argv[2]
rows, seen = [], False
default = "genomad end-to-end --cleanup --threads {threads} {in} {out} {db}"
if os.path.isfile(path):
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        if cols and cols[0].strip().lower() == "genomad" and not line.startswith("#name"):
            cmd = cols[1] if len(cols) > 1 and cols[1].strip() else default
            rows.append("genomad\t{}\t{}\t\n".format(cmd, db)); seen = True
        else:
            rows.append(line)
if not seen:
    rows.append("genomad\t{}\t{}\t\n".format(default, db))
open(path, "w").writelines(rows)
print("  updated the genomad row of {}".format(path))
PY

echo "DONE"
echo
echo "Use with FlaGs3:  --genomad -sr 50000"
