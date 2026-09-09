#!/usr/bin/env bash

set -e

ENV_NAME="flags3-genomad"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_TABLE="${HERE}/tools_table.local.tsv"

# Installers write here, never to the committed tools_table.tsv, so machine
# paths stay out of the repository. FlaGs3 reads the default table first and
# lets this one override it row by row.
if [[ ! -f "${TOOLS_TABLE}" ]]; then
    printf '#name\tcommand\tdirectory\tscan_range\n' > "${TOOLS_TABLE}"
fi

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<'USAGE'
Usage: genomad_installer.sh

Builds the conda environment flags3-genomad and installs geNomad into it.
geNomad is open source, so nothing has to be obtained first.

This installs the software only. Fetch its database separately with
genomad_loader.sh, which is a 1.6 GB download.
USAGE
    exit 0
fi

if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install miniconda first."
    exit 1
fi

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Environment ${ENV_NAME} already exists, reusing it"
else
    echo "Creating environment ${ENV_NAME}"
    conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda genomad
fi

env_prefix() {
    local name="$1" prefix=""
    prefix=$(conda run -n "${name}" printenv CONDA_PREFIX 2>/dev/null | tr -d '\r' | tail -1) || true
    if [[ -z "${prefix}" || ! -d "${prefix}" ]]; then
        prefix=$(conda env list | awk -v e="${name}" '$1==e {print $NF}' | tail -1) || true
    fi
    if [[ -z "${prefix}" || ! -d "${prefix}" ]]; then
        echo "ERROR: could not locate the ${name} environment. Try: conda env list" >&2
        exit 1
    fi
    printf '%s' "${prefix}"
}

PREFIX=$(env_prefix "${ENV_NAME}")
echo "  environment at ${PREFIX}"
BIN="${PREFIX}/bin/genomad"
if [[ ! -x "${BIN}" ]]; then
    echo "ERROR: ${BIN} is missing after the install."
    echo "       What the environment does contain:"
    ls "${PREFIX}/bin" 2>/dev/null | head -20 | sed "s/^/         /"
    exit 1
fi

echo "Checking geNomad's own dependencies"
MISSING=()
for dep in mmseqs aragorn; do
    [[ -x "${PREFIX}/bin/${dep}" ]] || MISSING+=("${dep}")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "  installing into ${ENV_NAME}: ${MISSING[*]}"
    declare -A PKG=( [mmseqs]=mmseqs2 [aragorn]=aragorn )
    TO_INSTALL=()
    for dep in "${MISSING[@]}"; do TO_INSTALL+=("${PKG[$dep]}"); done
    conda install -y -n "${ENV_NAME}" -c conda-forge -c bioconda "${TO_INSTALL[@]}"
fi
for dep in mmseqs aragorn; do
    if [[ -x "${PREFIX}/bin/${dep}" ]]; then
        echo "  ${dep}: ${PREFIX}/bin/${dep}"
    else
        echo "  ERROR: ${dep} is still missing from ${ENV_NAME}."
        echo "         geNomad shells out to it and will stop partway without it."
        echo "         Try: conda install -n ${ENV_NAME} -c bioconda ${dep}"
        exit 1
    fi
done

echo "Verification"
if "${BIN}" --version > /dev/null 2>&1; then
    echo "  genomad runs: $("${BIN}" --version 2>&1 | head -1)"
else
    echo "  ERROR: ${BIN} did not run:"
    "${BIN}" --version 2>&1 | tail -5 | sed "s/^/    /"
    exit 1
fi

DB=$(python3 - "${TOOLS_TABLE}" <<'PY'
import sys, os
path = sys.argv[1]
if os.path.isfile(path):
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        if cols and cols[0].strip().lower() == "genomad" and len(cols) > 2:
            print(cols[2].strip())
            break
PY
)

CMD="${BIN} end-to-end --cleanup --threads {threads} {in} {out} {db}"
python3 - "${TOOLS_TABLE}" "${CMD}" "${DB}" <<'PY'
import sys, os
path, cmd, db = sys.argv[1], sys.argv[2], sys.argv[3]
rows, seen = [], False
if os.path.isfile(path):
    for line in open(path):
        if line.split("\t")[0].strip().lstrip("#").lower() == "genomad" and not line.startswith("#name"):
            rows.append("genomad\t{}\t{}\t\n".format(cmd, db)); seen = True
        else:
            rows.append(line)
if not seen:
    rows.append("genomad\t{}\t{}\t\n".format(cmd, db))
open(path, "w").writelines(rows)
print("  updated the genomad row of {}".format(path))
PY

echo "DONE"
echo
if [[ -z "${DB}" ]]; then
    echo "Next: run genomad_loader.sh to fetch the database (1.6 GB)."
else
    echo "Database already set to ${DB}"
fi
echo "Then use with FlaGs3:  --genomad -sr 50000"
