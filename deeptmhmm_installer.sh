#!/usr/bin/env bash

set -e

ENV_NAME="flags3-deeptmhmm"
PY_VERSION="3.8"
PACKAGE="$1"
INSTALL_DIR="${2:-$(pwd)/deeptmhmm}"
TOOLS_TABLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools_table.tsv"

if [[ -z "${PACKAGE}" ]]; then
    cat <<'USAGE'
Usage: deeptmhmm_installer.sh /path/to/deeptmhmm-package[.tar.gz] [install_dir]

DeepTMHMM is licensed software and cannot be downloaded automatically.
Request the academic version by emailing licensing@biolib.com, then pass
the downloaded archive or unpacked directory to this script.

The package is installed into ./deeptmhmm unless install_dir is given.
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

if [[ -f "${PACKAGE}" ]]; then
    echo "Unpacking ${PACKAGE}"
    WORK=$(mktemp -d)
    unpack "${PACKAGE}" "${WORK}"
    SRC=$(find "${WORK}" -maxdepth 6 -name predict.py -exec dirname {} \; | head -1)
else
    SRC=$(find "${PACKAGE}" -maxdepth 6 -name predict.py -exec dirname {} \; | head -1)
fi

if [[ -z "${SRC}" ]]; then
    echo "ERROR: no predict.py found in ${PACKAGE}."
    [[ -n "${WORK:-}" ]] && rm -rf "${WORK}"
    exit 1
fi

echo "Installing to ${INSTALL_DIR}"
mkdir -p "${INSTALL_DIR}"
cp -r "${SRC}"/. "${INSTALL_DIR}/"
[[ -n "${WORK:-}" ]] && rm -rf "${WORK}"
INSTALL_DIR=$(cd "${INSTALL_DIR}" && pwd)

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Environment ${ENV_NAME} already exists, reusing it"
else
    echo "Creating environment ${ENV_NAME} (python ${PY_VERSION})"
    conda create -y -n "${ENV_NAME}" "python=${PY_VERSION}"
fi

echo "Installing build dependencies"
conda run --no-capture-output -n "${ENV_NAME}" python -m pip install \
    wheel "Cython==0.29.37" "pkgconfig==1.5.5"

echo "Installing PyTorch"
TORCH_INDEX="https://download.pytorch.org/whl/torch_stable.html"
LINUX64=""
[[ "$(uname -s)" == "Linux" && "$(uname -m)" == "x86_64" ]] && LINUX64="yes"

HAVE=$(conda run -n "${ENV_NAME}" python -c \
    "import torch; print(torch.__version__)" 2>/dev/null | tail -1) || true

if [[ "${DEEPTMHMM_GPU:-0}" == "1" ]]; then
    WANT="cu92"
else
    WANT="cpu"
fi

if [[ "${HAVE}" == *"${WANT}"* ]]; then
    echo "  torch ${HAVE} already installed"
    TORCH_OK="${WANT}"
else
    [[ -n "${HAVE}" ]] && echo "  replacing torch ${HAVE} with the ${WANT} build"
    TORCH_OK=""
    if [[ -n "${LINUX64}" ]]; then
        conda run --no-capture-output -n "${ENV_NAME}" python -m pip install \
            --force-reinstall "torch==1.5.0+${WANT}" -f "${TORCH_INDEX}" \
            && TORCH_OK="${WANT}" || true
    fi
    if [[ -z "${TORCH_OK}" ]]; then
        echo "  falling back to torch 1.5.0 from PyPI"
        conda run --no-capture-output -n "${ENV_NAME}" python -m pip install \
            "torch==1.5.0" && TORCH_OK="pypi" || true
    fi
    if [[ -z "${TORCH_OK}" ]]; then
        echo "ERROR: could not install torch 1.5.0 for this platform."
        echo "Install it by hand following the package README, then re-run."
        exit 1
    fi
    echo "  installed the ${TORCH_OK} build"
fi
if [[ "${WANT}" == "cpu" ]]; then
    echo "  CPU build by default: the CUDA 9.2 build named in the package README"
    echo "  predates current GPUs and fails with a cublas error on them."
    echo "  Set DEEPTMHMM_GPU=1 to use it anyway."
fi

if [[ -f "${INSTALL_DIR}/requirements.txt" ]]; then
    echo "Installing package requirements"
    REQS=$(mktemp)
    grep -viE "^[[:space:]]*torch([=<>!~[:space:]]|$)" \
        "${INSTALL_DIR}/requirements.txt" > "${REQS}" || true
    conda run --no-capture-output -n "${ENV_NAME}" python -m pip install -r "${REQS}"
    rm -f "${REQS}"
fi

PYBIN=$(conda run -n "${ENV_NAME}" python -c \
    "import sys, os; print(os.path.join(sys.prefix, 'bin', 'python3'))" | tail -1)

echo "Verification"
if (cd "${INSTALL_DIR}" && "${PYBIN}" -c "import torch, esm" 2>/dev/null); then
    echo "  torch and fair-esm import cleanly"
else
    echo "  ERROR: torch or fair-esm did not import in ${ENV_NAME}:"
    (cd "${INSTALL_DIR}" && "${PYBIN}" -c "import torch, esm" 2>&1 | tail -5 \
        | sed "s/^/    /") || true
    exit 1
fi

if [[ -f "${INSTALL_DIR}/sample.fasta" ]]; then
    TESTOUT="$(mktemp -d)/out"
    echo "Running predict.py on the bundled sample"
    echo "  this also downloads any model weights, so later runs need no network"
    echo "  output follows; it is not hidden"
    echo "  ----------------------------------------------------------------"
    GPU_ENV=""
    [[ "${DEEPTMHMM_GPU:-0}" != "1" ]] && GPU_ENV="env CUDA_VISIBLE_DEVICES="
    RUNNER=""
    command -v timeout &> /dev/null && RUNNER="timeout ${DEEPTMHMM_TIMEOUT:-3600}"
    set +e
    (
        set -o pipefail
        cd "${INSTALL_DIR}" && ${GPU_ENV} ${RUNNER} \
            "${PYBIN}" -u predict.py --fasta sample.fasta --output-dir "${TESTOUT}" \
            2>&1 | sed "s/^/  | /"
    )
    STATUS=$?
    set -e
    echo "  ----------------------------------------------------------------"

    RESULT=$(find "${TESTOUT}" -name "predicted_topologies.3line" | head -1)
    if [[ "${STATUS}" -eq 124 ]]; then
        echo "  ERROR: predict.py was still running after ${DEEPTMHMM_TIMEOUT:-3600}s."
        echo "  Raise it with DEEPTMHMM_TIMEOUT=7200 bash deeptmhmm_installer.sh ..."
        exit 1
    elif [[ "${STATUS}" -ne 0 ]]; then
        echo "  ERROR: predict.py exited ${STATUS}. The output above is the reason."
        exit 1
    elif [[ -z "${RESULT}" ]]; then
        echo "  ERROR: predict.py exited 0 but wrote no predicted_topologies.3line."
        echo "  It produced:"
        find "${TESTOUT}" -type f | sed "s/^/    /" | head -10
        echo "  FlaGs3 looks for predicted_topologies.3line; if this build names it"
        echo "  differently, that name has to change in flags_features.py."
        exit 1
    fi
    echo "  predict.py wrote $(basename "${RESULT}"), $(grep -c "^>" "${RESULT}") record(s)"
    rm -rf "${TESTOUT}"
else
    echo "  WARNING: no sample.fasta in the package, cannot test predict.py"
fi

CMD="${PYBIN} predict.py --fasta {fasta} --output-dir {out}"
if [[ -f "${TOOLS_TABLE}" ]]; then
    python3 - "${TOOLS_TABLE}" "${CMD}" "${INSTALL_DIR}" <<'PY'
import sys
path, cmd, wd = sys.argv[1], sys.argv[2], sys.argv[3]
rows, seen = [], False
for line in open(path):
    if line.split("\t")[0].strip().lstrip("#").lower() == "deeptmhmm" and not line.startswith("#name"):
        rows.append("deeptmhmm\t{}\t{}\n".format(cmd, wd)); seen = True
    else:
        rows.append(line)
if not seen:
    rows.append("deeptmhmm\t{}\t{}\n".format(cmd, wd))
open(path, "w").writelines(rows)
print("  updated the deeptmhmm row of {}".format(path))
PY
else
    echo "  add this row to tools_table.tsv:"
    printf '    deeptmhmm\t%s\t%s\n' "${CMD}" "${INSTALL_DIR}"
fi

echo "DONE"
echo
echo "Use with FlaGs3:"
echo "  --tmhmm --local_tmhmm"
