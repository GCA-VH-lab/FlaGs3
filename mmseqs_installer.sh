#!/usr/bin/env bash

set -e

ENV_NAME="flags3-mmseqs"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_TABLE="${HERE}/tools_table.tsv"
FALLBACK_DIR="${HERE}/mmseqs"

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<'USAGE'
Usage: mmseqs_installer.sh [--binary]

Installs MMseqs2, used by FlaGs3's --cluster_collapse to reduce flanking
proteins to representatives before clustering. MMseqs2 is open source, so
this script downloads it for you; nothing has to be obtained first.

By default it builds the conda environment flags3-mmseqs. With --binary, or
when conda is missing, it downloads the official static build into ./mmseqs
instead, which needs no conda and no compiler.

Either way the mmseqs row of tools_table.tsv is rewritten to point at what
was installed.
USAGE
    exit 0
fi

BIN=""

install_binary() {
    echo "Installing the static MMseqs2 build"
    local flavour=""
    if grep -qm1 avx2 /proc/cpuinfo 2>/dev/null; then
        flavour="avx2"
    elif grep -qm1 sse4_1 /proc/cpuinfo 2>/dev/null; then
        flavour="sse41"
    fi
    case "$(uname -s)" in
        Darwin) asset="mmseqs-osx-universal.tar.gz" ;;
        Linux)
            case "$(uname -m)" in
                aarch64|arm64) asset="mmseqs-linux-arm64.tar.gz" ;;
                *) asset="mmseqs-linux-${flavour:-sse2}.tar.gz" ;;
            esac
            ;;
        *)
            echo "ERROR: unsupported platform $(uname -s); install mmseqs yourself"
            echo "       and put its path in the mmseqs row of tools_table.tsv."
            exit 1
            ;;
    esac
    echo "  platform: $(uname -s) $(uname -m), using ${asset}"

    local url="https://github.com/soedinglab/MMseqs2/releases/latest/download/${asset}"
    local work
    work=$(mktemp -d)
    if command -v curl &> /dev/null; then
        curl -fsSL -o "${work}/mmseqs.tar.gz" "${url}"
    elif command -v wget &> /dev/null; then
        wget -qO "${work}/mmseqs.tar.gz" "${url}"
    else
        echo "ERROR: neither curl nor wget is available."
        rm -rf "${work}"
        exit 1
    fi
    rm -rf "${FALLBACK_DIR}"
    mkdir -p "$(dirname "${FALLBACK_DIR}")"
    tar xzf "${work}/mmseqs.tar.gz" -C "$(dirname "${FALLBACK_DIR}")"
    rm -rf "${work}"
    BIN="${FALLBACK_DIR}/bin/mmseqs"
    if [[ ! -x "${BIN}" ]]; then
        echo "ERROR: the archive did not contain bin/mmseqs"
        exit 1
    fi
}

install_conda() {
    if conda env list | grep -qE "^${ENV_NAME}\s"; then
        echo "Environment ${ENV_NAME} already exists, reusing it"
    else
        echo "Creating environment ${ENV_NAME}"
        conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda mmseqs2
    fi
    if ! conda env list | grep -qE "^${ENV_NAME}\s"; then
        return 1
    fi
    local prefix
    prefix=$(conda run -n "${ENV_NAME}" python -c \
        "import sys; print(sys.prefix)" 2>/dev/null | tail -1) || true
    if [[ -z "${prefix}" ]]; then
        prefix=$(conda env list | awk -v e="${ENV_NAME}" '$1==e {print $NF}')
    fi
    BIN="${prefix}/bin/mmseqs"
    [[ -x "${BIN}" ]]
}

if [[ "$1" == "--binary" ]]; then
    install_binary
elif command -v conda &> /dev/null; then
    if ! install_conda; then
        echo "  conda install did not produce a usable mmseqs, falling back"
        install_binary
    fi
else
    echo "conda not found, using the static build instead"
    install_binary
fi

echo "Verification"
if "${BIN}" version > /dev/null 2>&1; then
    echo "  mmseqs runs: $("${BIN}" version 2>&1 | head -1)"
else
    echo "  ERROR: ${BIN} did not run:"
    "${BIN}" version 2>&1 | tail -5 | sed "s/^/    /"
    exit 1
fi

echo "Self-test on a small set"
TEST=$(mktemp -d)
cat > "${TEST}/in.faa" <<'FASTA'
>a
MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRS
>b
MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRA
>c
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHL
FASTA
if "${BIN}" easy-linclust "${TEST}/in.faa" "${TEST}/clu" "${TEST}/tmp" \
        --min-seq-id 0.9 -c 0.8 --cov-mode 0 -v 1 > "${TEST}/log" 2>&1; then
    REPS=$(grep -c ">" "${TEST}/clu_rep_seq.fasta" || echo 0)
    echo "  3 sequences collapsed to ${REPS} representatives (2 expected)"
else
    echo "  ERROR: the self-test failed:"
    tail -10 "${TEST}/log" | sed "s/^/    /"
    rm -rf "${TEST}"
    exit 1
fi
rm -rf "${TEST}"

CMD="${BIN} easy-linclust {in} {out} {tmp} --min-seq-id {id} -c {cov} --cov-mode 0 --threads {threads} -v 1"
if [[ -f "${TOOLS_TABLE}" ]]; then
    python3 - "${TOOLS_TABLE}" "${CMD}" <<'PY'
import sys
path, cmd = sys.argv[1], sys.argv[2]
rows, seen = [], False
for line in open(path):
    if line.split("\t")[0].strip().lstrip("#").lower() == "mmseqs" and not line.startswith("#name"):
        rows.append("mmseqs\t{}\t\n".format(cmd)); seen = True
    else:
        rows.append(line)
if not seen:
    rows.append("mmseqs\t{}\t\n".format(cmd))
open(path, "w").writelines(rows)
print("  updated the mmseqs row of {}".format(path))
PY
else
    echo "  add this row to tools_table.tsv:"
    printf '    mmseqs\t%s\t\n' "${CMD}"
fi

echo "DONE"
echo
echo "Use with FlaGs3:"
echo "  --cluster_collapse            collapse at 90% identity / 80% coverage"
echo "  --cluster_collapse 0.5,0.8    collapse harder, for very large runs"
echo
echo "Swap easy-linclust for easy-cluster in tools_table.tsv if you want the"
echo "slower, more sensitive algorithm at low identity."
