#!/usr/bin/env bash

set -e

BASE="./defence_db"
REPO="https://github.com/mdmparis/defense-finder-models"
API="https://api.github.com/repos/mdmparis/defense-finder-models/releases/latest"

echo "Environment Check"

if command -v curl &> /dev/null; then
    DOWNLOAD_CMD="curl -fL -o"
    FETCH_STDOUT="curl -fsSL"
elif command -v wget &> /dev/null; then
    DOWNLOAD_CMD="wget -O"
    FETCH_STDOUT="wget -qO-"
else
    echo "ERROR: Neither 'curl' nor 'wget' was found. Please install one of them."
    exit 1
fi

echo "Finding the latest release"
RELEASE=$($FETCH_STDOUT "${API}" 2>/dev/null) || true
VERSION=$(printf '%s' "${RELEASE}" | grep -o '"tag_name": *"[^"]*"' | head -1 \
          | sed 's/.*"\([^"]*\)"$/\1/') || true
TARBALL=$(printf '%s' "${RELEASE}" \
          | grep -o '"browser_download_url": *"[^"]*\.tar\.gz"' | head -1 \
          | sed 's/.*"\(https[^"]*\)"/\1/') || true

if [[ -z "${TARBALL}" ]]; then
    echo "  release lookup failed (GitHub API rate limit is common); using the"
    echo "  master branch and taking the version from the archive"
    TARBALL="${REPO}/archive/refs/heads/master.tar.gz"
    VERSION="master"
fi

DF_DIR="${BASE}_${VERSION}"
echo "  version ${VERSION} -> ${DF_DIR}"

EXISTING=""
if [[ -d "${DF_DIR}" ]] && compgen -G "${DF_DIR}/*.hmm" > /dev/null; then
    EXISTING="${DF_DIR}"
elif [[ "${VERSION}" == "master" ]]; then
    EXISTING=$(compgen -G "${BASE}_*" 2>/dev/null | while read -r d; do
        compgen -G "${d}/*.hmm" > /dev/null && echo "${d}"; done | sort -V | tail -1)
    if [[ -n "${EXISTING}" ]]; then
        DF_DIR="${EXISTING}"
        echo "  version could not be confirmed; using what is already installed"
    fi
fi

if [[ -n "${EXISTING}" ]]; then
    echo "  already present with $(ls "${DF_DIR}"/*.hmm | wc -l) models, skipping download"
    echo "  delete ${DF_DIR} to force a fresh download"
else
    OTHER=$(compgen -G "${BASE}_*" 2>/dev/null | grep -v "^${DF_DIR}$" || true)
    if [[ -n "${OTHER}" ]]; then
        echo "  other versions present, left untouched:"
        printf '    %s\n' ${OTHER}
    fi

    echo "  downloading ${TARBALL}"
    $DOWNLOAD_CMD models.tar.gz "${TARBALL}"

    echo "Extracting profiles"
    rm -rf _unpack && mkdir _unpack
    tar -xzf models.tar.gz -C _unpack

    SRC=$(find _unpack -type d -name profiles | head -1)
    if [[ -z "${SRC}" ]]; then
        echo "ERROR: no profiles/ directory inside the archive."
        rm -rf _unpack models.tar.gz
        exit 1
    fi
    if [[ "${VERSION}" == "master" ]]; then
        META=$(find _unpack -maxdepth 2 -name metadata.yml | head -1)
        FROM_META=$(grep -oE "^vers(ion)?: *[0-9][^ ]*" "${META}" 2>/dev/null \
                    | sed "s/.*: *//") || true
        if [[ -n "${FROM_META}" ]]; then
            VERSION="${FROM_META}"
            DF_DIR="${BASE}_${VERSION}"
            echo "  version ${VERSION} read from the archive -> ${DF_DIR}"
        fi
    fi
    rm -rf "${DF_DIR}"
    mkdir -p "${DF_DIR}"
    find "${SRC}" -name '*.hmm' -exec mv {} "${DF_DIR}/" \;
    rm -rf _unpack models.tar.gz
fi

COUNT=$(ls "${DF_DIR}"/*.hmm 2>/dev/null | wc -l)
if [[ "${COUNT}" -eq 0 ]]; then
    echo "ERROR: no .hmm files were extracted."
    exit 1
fi

echo "Verification"
echo "  ${COUNT} HMM profiles in ${DF_DIR}"
du -sh "${DF_DIR}"

echo "DONE"
echo
echo "Use with FlaGs3:"
echo "  --domains --hmmdb defence=${DF_DIR} --hmm_coverage defence=0.3,0.3"
echo
echo "These are full-length protein models, so a coverage cutoff is worth setting;"
echo "without one you also get short partial matches. Raise it to tighten."
echo "Pfam can be given at the same time:"
echo "  --hmmdb pfam_db/Pfam-A.hmm --hmmdb defence=${DF_DIR}"
