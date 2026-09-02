#!/usr/bin/env bash

set -e

ENV_NAME="flags3-signalp"
PY_VERSION="3.9"
PACKAGE="$1"
TOOLS_TABLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools_table.tsv"

if [[ -z "${PACKAGE}" ]]; then
    cat <<'USAGE'
Usage: signalp_installer.sh /path/to/signalp-6-package[.tar.gz]

SignalP 6 is licensed software and cannot be downloaded automatically.
Request it (free for academic use) at:
  https://services.healthtech.dtu.dk/services/SignalP-6.0/
then pass the downloaded archive or unpacked directory to this script.
USAGE
    exit 1
fi

unpack() {
    local archive="$1" dest="$2"
    python3 - "${archive}" "${dest}" <<'PY'
import shutil, sys
try:
    shutil.unpack_archive(sys.argv[1], sys.argv[2])
except Exception as exc:
    sys.exit("could not unpack {}: {}".format(sys.argv[1], exc))
PY
}

if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install miniconda first."
    exit 1
fi

if [[ ! -e "${PACKAGE}" ]]; then
    echo "ERROR: no such file or directory: ${PACKAGE}"
    exit 1
fi

WORK=""
if [[ -f "${PACKAGE}" ]]; then
    echo "Unpacking ${PACKAGE}"
    WORK=$(mktemp -d)
    unpack "${PACKAGE}" "${WORK}"
    SRC=$(find "${WORK}" -maxdepth 6 -type d -name "signalp*package*" | head -1)
    [[ -z "${SRC}" ]] && SRC=$(find "${WORK}" -maxdepth 6 -name setup.py -exec dirname {} \; | head -1)
else
    SRC="${PACKAGE}"
fi

if [[ -z "${SRC}" ]] || [[ ! -e "${SRC}/setup.py" && ! -e "${SRC}/pyproject.toml" ]]; then
    echo "ERROR: ${PACKAGE} does not look like the signalp-6-package."
    [[ -n "${WORK}" ]] && rm -rf "${WORK}"
    exit 1
fi
echo "  package: ${SRC}"

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Environment ${ENV_NAME} already exists, reusing it"
else
    echo "Creating environment ${ENV_NAME} (python ${PY_VERSION})"
    conda create -y -n "${ENV_NAME}" "python=${PY_VERSION}"
fi

echo "Installing SignalP 6"
conda run --no-capture-output -n "${ENV_NAME}" pip install "torch<2.0" "numpy<2"
conda run --no-capture-output -n "${ENV_NAME}" pip install "${SRC}/"

echo "Checking numpy"
NPV=$(conda run -n "${ENV_NAME}" python -c \
    "import numpy; print(numpy.__version__)" 2>/dev/null | tail -1) || true
case "${NPV}" in
    2.*)
        echo "  installing the package pulled numpy ${NPV}; torch 1.x needs numpy 1.x"
        conda run --no-capture-output -n "${ENV_NAME}" pip install "numpy<2"
        ;;
    *)
        echo "  numpy ${NPV}"
        ;;
esac

echo "Copying model weights"
TARGET=$(conda run -n "${ENV_NAME}" python -c \
    "import signalp, os; print(os.path.dirname(signalp.__file__))" | tail -1)
if [[ -d "${SRC}/models" ]]; then
    mkdir -p "${TARGET}/model_weights"
    cp -r "${SRC}"/models/* "${TARGET}/model_weights/"
elif [[ -d "${SRC}/signalp/model_weights" ]]; then
    mkdir -p "${TARGET}/model_weights"
    cp -r "${SRC}"/signalp/model_weights/* "${TARGET}/model_weights/"
else
    echo "  no models/ directory in the package; if signalp6 reports missing"
    echo "  weights, copy them into ${TARGET}/model_weights/ by hand"
fi

BIN=$(conda run -n "${ENV_NAME}" python -c \
    "import sys, os; print(os.path.join(sys.prefix, 'bin'))" | tail -1)
[[ -n "${WORK}" ]] && rm -rf "${WORK}"

echo "Verification"
if "${BIN}/signalp6" --version > /dev/null 2>&1; then
    echo "  signalp6 runs: $("${BIN}/signalp6" --version 2>&1 | head -1)"
else
    echo "  ERROR: ${BIN}/signalp6 did not run:"
    "${BIN}/signalp6" --version 2>&1 | tail -5 | sed "s/^/    /"
    exit 1
fi

CMD="${BIN}/signalp6 --fastafile {fasta} --output_dir {out} --organism other --format txt --mode fast"
if [[ -f "${TOOLS_TABLE}" ]]; then
    python3 - "${TOOLS_TABLE}" "${CMD}" <<'PY'
import sys
path, cmd = sys.argv[1], sys.argv[2]
rows, seen = [], False
for line in open(path):
    if line.split("\t")[0].strip().lstrip("#").lower() == "signalp" and not line.startswith("#name"):
        rows.append("signalp\t{}\t\n".format(cmd)); seen = True
    else:
        rows.append(line)
if not seen:
    rows.append("signalp\t{}\t\n".format(cmd))
open(path, "w").writelines(rows)
print("  updated the signalp row of {}".format(path))
PY
else
    echo "  add this row to tools_table.tsv:"
    printf '    signalp\t%s\t\n' "${CMD}"
fi

echo "DONE"
echo
echo "Use with FlaGs3:"
echo "  --signalp --local_signalp"
