#!/usr/bin/env bash

set -e

TOOL="padloc"
ENV_NAME="flags3-${TOOL}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_TABLE="${HERE}/tools_table.local.tsv"

# Installers write here, never to the committed tools_table.tsv, so machine
# paths stay out of the repository. FlaGs3 reads the default table first and
# lets this one override it row by row.
if [[ ! -f "${TOOLS_TABLE}" ]]; then
    printf '#name\tcommand\tdirectory\tscan_range\n' > "${TOOLS_TABLE}"
fi

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: ${TOOL}_installer.sh"
    echo
    echo "Builds the conda environment ${ENV_NAME} and installs ${TOOL} into it,"
    echo "including its models. Both tools are open source, so nothing has to be"
    echo "obtained first."
    echo
    if [[ "${TOOL}" == "defensefinder" ]]; then
        echo "This is the whole DefenseFinder tool, which calls complete defence"
        echo "systems. It is separate from defensefinder_hmm_loader.sh, which only"
        echo "fetches the HMM profiles for the --hmm_db route."
    fi
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
    if [[ "${TOOL}" == "defensefinder" ]]; then
        conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda python=3.10 hmmer
        conda run -n "${ENV_NAME}" pip install mdmparis-defense-finder
    else
        conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda padloc
    fi
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

if [[ "${TOOL}" == "defensefinder" ]]; then
    BIN="${PREFIX}/bin/defense-finder"
else
    BIN="${PREFIX}/bin/padloc"
fi

echo "Verification"
if [[ -x "${BIN}" ]]; then
    echo "  ${TOOL} installed at ${BIN}"
else
    echo "  ERROR: ${BIN} is missing after the install."
    echo "         What the environment does contain:"
    ls "${PREFIX}/bin" 2>/dev/null | head -20 | sed "s/^/           /"
    exit 1
fi

if [[ "${TOOL}" == "defensefinder" ]]; then
    echo "Fetching DefenseFinder models"
    "${BIN}" update || echo "  WARNING: 'defense-finder update' failed; run it by hand"
    CMD="${BIN} run --db-type gembase -o {out} {faa}"
else
    echo "Fetching PadLoc database"
    "${BIN}" --db-update || echo "  WARNING: 'padloc --db-update' failed; run it by hand"
    CMD="${BIN} --faa {faa} --gff {gff} --outdir {out} --cpu {threads}"
fi

python3 - "${TOOLS_TABLE}" "${TOOL}" "${CMD}" <<'PY'
import sys, os
path, tool, cmd = sys.argv[1], sys.argv[2], sys.argv[3]
rows, seen = [], False
if os.path.isfile(path):
    for line in open(path):
        if line.split("\t")[0].strip().lstrip("#").lower() == tool and not line.startswith("#name"):
            rows.append("{}\t{}\t\t\n".format(tool, cmd)); seen = True
        else:
            rows.append(line)
if not seen:
    rows.append("{}\t{}\t\t\n".format(tool, cmd))
open(path, "w").writelines(rows)
print("  updated the {} row of {}".format(tool, path))
PY

echo "DONE"
echo
if [[ "${TOOL}" == "defensefinder" ]]; then
    echo "Use with FlaGs3:  --defensefinder"
else
    echo "Use with FlaGs3:  --padloc"
fi
echo "Both can run together; a system they agree on is drawn once."
