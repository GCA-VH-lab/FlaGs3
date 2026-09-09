#!/usr/bin/env bash
set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info()  { printf "${GREEN}[INFO]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
error() { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; }
say()   { printf '%b\n' "$*"; }

ask() {
    printf "\n"
    for line in "$@"; do printf "  %s\n" "${line}"; done
    printf "  [y/N] "
    read -r REPLY
    case "${REPLY}" in
        [Yy]|[Yy][Ee][Ss]) return 0 ;;
        *) return 1 ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    error "Do not source this script. Run it directly: bash build.sh"
    return 1
fi

COMPONENTS="pfam defence-hmm mmseqs genomad defensefinder padloc signalp deeptmhmm"
SELECTED=""
ASSUME=""
SIGNALP_PKG=""
DEEPTMHMM_PKG=""

usage() {
    cat <<USAGE
Usage: bash build.sh [options]

Builds the FlaGs3 conda environment, then offers each optional tool in turn.
With no options every optional tool is offered interactively.

  --all                 install every optional tool without asking
  --none                install none of them, environment only
  --with LIST           install exactly these, comma separated, no prompts
  --signalp PATH        SignalP 6 package tarball, for --with signalp
  --deeptmhmm PATH      DeepTMHMM package tarball, for --with deeptmhmm
  -h, --help            this text

Components for --with:
  pfam            Pfam-A profiles, for --domains          (~1.5 GB)
  defence-hmm     DefenseFinder HMM profiles, for --domains --hmmdb defence=
  mmseqs          MMseqs2, for --cluster_collapse
  genomad         geNomad and its database, for --genomad (~1.6 GB)
  defensefinder   DefenseFinder itself, for --defensefinder
  padloc          PadLoc, for --padloc
  signalp         SignalP 6 locally, needs --signalp PATH  (licensed)
  deeptmhmm       DeepTMHMM locally, needs --deeptmhmm PATH (licensed)

Examples:
  bash build.sh --all
  bash build.sh --with mmseqs,genomad,defensefinder,padloc
  bash build.sh --none
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --all)  ASSUME="all"; shift ;;
        --none) ASSUME="none"; shift ;;
        --with) SELECTED="${2:-}"; shift 2 ;;
        --signalp)   SIGNALP_PKG="${2:-}"; shift 2 ;;
        --deeptmhmm) DEEPTMHMM_PKG="${2:-}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -n "${SELECTED}" ]]; then
    for item in ${SELECTED//,/ }; do
        if [[ " ${COMPONENTS} " != *" ${item} "* ]]; then
            error "Unknown component: ${item}"
            error "Known: ${COMPONENTS}"
            exit 1
        fi
    done
fi

# Decide whether a component is wanted: --all/--none win, then an explicit
# --with list, otherwise fall back to asking.
want() {
    local name="$1"; shift
    case "${ASSUME}" in
        all)  info "Selected by --all: ${name}"; return 0 ;;
        none) return 1 ;;
    esac
    if [[ -n "${SELECTED}" ]]; then
        [[ ",${SELECTED}," == *",${name},"* ]] || return 1
        info "Selected by --with: ${name}"
        return 0
    fi
    ask "$@"
}

# Run one installer inside the FlaGs3 environment and report clearly.
run_script() {
    local label="$1" script="$2"; shift 2
    if [[ ! -f "${script}" ]]; then
        error "${label}: ${script} not found next to build.sh"
        return 1
    fi
    info "Running $(basename "${script}") ..."
    if conda run --no-capture-output --name "${ENV_NAME}" bash "${script}" "$@"; then
        info "${label}: done"
        return 0
    fi
    warn "${label}: the installer failed. Run it by hand to see why:"
    warn "  conda run --name ${ENV_NAME} bash $(basename "${script}") $*"
    return 1
}

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${THIS_DIR}/environment.yml"
PFAM_SCRIPT="${THIS_DIR}/pfamA_loader.sh"
DF_SCRIPT="${THIS_DIR}/defensefinder_hmm_loader.sh"
MMSEQS_SCRIPT="${THIS_DIR}/mmseqs_installer.sh"
GENOMAD_SCRIPT="${THIS_DIR}/genomad_installer.sh"
GENOMAD_DB_SCRIPT="${THIS_DIR}/genomad_loader.sh"
DEFENSEFINDER_SCRIPT="${THIS_DIR}/defensefinder_installer.sh"
PADLOC_SCRIPT="${THIS_DIR}/padloc_installer.sh"
SIGNALP_SCRIPT="${THIS_DIR}/signalp_installer.sh"
DEEPTMHMM_SCRIPT="${THIS_DIR}/deeptmhmm_installer.sh"
ENV_NAME="FlaGs3"

if [[ "$(uname)" == "Darwin" ]]; then
    unset DYLD_LIBRARY_PATH || true
fi


info "Checking for Conda..."

if ! command -v conda &>/dev/null; then
    error "Conda executable not found in PATH."
    error "Please install Miniconda or Anaconda and re-run this script."
    error "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

info "Found: $(conda --version 2>&1)"


if [[ ! -f "${ENV_FILE}" ]]; then
    error "Environment file not found: ${ENV_FILE}"
    error "Expected it next to this script: ${THIS_DIR}/environment.yml"
    exit 1
fi

info "Using environment file: ${ENV_FILE}"
info "Detected platform: $(uname -s) / $(uname -m)"


if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    warn "Environment '${ENV_NAME}' already exists."
    if ask "Remove and recreate it?"; then
        info "Removing existing environment '${ENV_NAME}'..."
        conda env remove --name "${ENV_NAME}" --yes
    else
        info "Keeping existing environment. Skipping creation."
        info "Activate it with:  conda activate ${ENV_NAME}"
        exit 0
    fi
fi


info "Creating Conda environment '${ENV_NAME}' — this may take a few minutes..."

if ! conda env create --name "${ENV_NAME}" --file "${ENV_FILE}"; then
    error "Conda environment creation failed."
    error "Possible fixes:"
    error "  • Check your internet connection."
    error "  • Ensure the bioconda channel is reachable."
    error "  • Run:  conda clean --all  and retry."
    exit 1
fi

info "Environment '${ENV_NAME}' created successfully."


info "Verifying the installation..."

if ! conda run --name "${ENV_NAME}" python -c "import Bio, requests, pyhmmer" &>/dev/null; then
    error "Post-install check failed: core Python packages not importable in '${ENV_NAME}'."
    error "The environment may be incomplete. Remove and re-run:"
    error "  conda env remove --name ${ENV_NAME} && bash build.sh"
    exit 1
fi
info "Core packages OK (Bio, requests, pyhmmer) — default pipeline is ready."

for tool in mafft VeryFastTree; do
    if conda run --name "${ENV_NAME}" command -v "${tool}" &>/dev/null; then
        info "Found tree tool: ${tool}"
    else
        warn "Tree tool '${tool}' not found — --tree / --tree_order will be skipped."
    fi
done

if conda run --name "${ENV_NAME}" python -c "import biolib" &>/dev/null; then
    info "Found pybiolib — --tmhmm / --signalp available."
else
    warn "pybiolib not importable — --tmhmm / --signalp will be skipped."
    warn "  Install with:  conda run --name ${ENV_NAME} pip install pybiolib"
fi

if conda run --name "${ENV_NAME}" python -c "import sismis" &>/dev/null; then
    info "Found sismis — --sismis available."
elif ask "Install sismis now? It is only needed for secretion-system" \
          "detection (the --sismis option) and pulls in gecco," \
          "scikit-learn, scipy, polars and pyrodigal."; then
    info "Installing sismis..."
    if ! conda run --no-capture-output --name "${ENV_NAME}" pip install sismis; then
        error "sismis installation failed. Check the output above."
        error "You can retry at any time:"
        error "  conda run --name ${ENV_NAME} pip install sismis"
        exit 1
    fi
    if ! conda run --name "${ENV_NAME}" python -c "import sismis, pyhmmer" &>/dev/null; then
        error "sismis installed but the environment is now inconsistent"
        error "(sismis or pyhmmer no longer imports). Remove and re-run:"
        error "  conda env remove --name ${ENV_NAME} && bash build.sh"
        exit 1
    fi
    info "sismis installed successfully."
else
    info "Skipping sismis. Install it later with:"
    info "  conda run --name ${ENV_NAME} pip install sismis"
fi

# --- optional tools -----------------------------------------------------
# Each is offered in turn. want() honours --all / --none / --with and falls
# back to asking, so every component behaves the same way.

PFAM_HMM="${THIS_DIR}/pfam_db/Pfam-A.hmm"
if [[ -f "${PFAM_HMM}" && -f "${PFAM_HMM}.h3m" ]]; then
    info "Pfam-A already installed --domains is available."
elif want pfam "Install the Pfam-A database? Needed for --domains. ~1.5 GB."; then
    run_script "Pfam-A" "${PFAM_SCRIPT}" || true
else
    info "Skipping Pfam-A.  Later:  bash pfamA_loader.sh"
fi

DF_PROFILES="${THIS_DIR}/defence_db"
if [[ -d "${DF_PROFILES}" ]] && compgen -G "${DF_PROFILES}/*.hmm" > /dev/null; then
    info "DefenseFinder HMM profiles already installed."
elif want defence-hmm "Install the DefenseFinder HMM profiles? They let --domains" \
                      "annotate defence systems alongside Pfam. ~50 MB." ; then
    run_script "DefenseFinder profiles" "${DF_SCRIPT}" || true
else
    info "Skipping the profiles.  Later:  bash defensefinder_hmm_loader.sh"
fi

if conda run --name "${ENV_NAME}" command -v mmseqs &>/dev/null \
   || [[ -x "${THIS_DIR}/mmseqs/bin/mmseqs" ]]; then
    info "MMseqs2 already installed --cluster_collapse is available."
elif want mmseqs "Install MMseqs2? Needed for --cluster_collapse, which makes" \
                 "clustering affordable on very large inputs. ~20 MB."; then
    run_script "MMseqs2" "${MMSEQS_SCRIPT}" || true
else
    info "Skipping MMseqs2.  Later:  bash mmseqs_installer.sh"
fi

if want genomad "Install geNomad and its database? Needed for --genomad," \
                "which finds proviruses and plasmids. Database is ~1.6 GB."; then
    if run_script "geNomad" "${GENOMAD_SCRIPT}"; then
        run_script "geNomad database" "${GENOMAD_DB_SCRIPT}" || true
    fi
else
    info "Skipping geNomad.  Later:  bash genomad_installer.sh && bash genomad_loader.sh"
fi

if want defensefinder "Install DefenseFinder? Needed for --defensefinder," \
                      "which calls whole defence systems rather than single domains."; then
    run_script "DefenseFinder" "${DEFENSEFINDER_SCRIPT}" || true
else
    info "Skipping DefenseFinder.  Later:  bash defensefinder_installer.sh"
fi

if want padloc "Install PadLoc? Needed for --padloc, a second defence system" \
               "caller that can run alongside DefenseFinder."; then
    run_script "PadLoc" "${PADLOC_SCRIPT}" || true
else
    info "Skipping PadLoc.  Later:  bash padloc_installer.sh"
fi

# SignalP and DeepTMHMM are licensed, so they need a package the user already
# has. Non-interactive runs never stop to ask for a path: with --all or --with
# and no package given, they are skipped with a note.
licensed() {
    local name="$1" script="$2" package="$3" prompt="$4" hint="$5"
    if [[ -n "${package}" ]]; then
        run_script "${name}" "${script}" "${package}" || true
        return
    fi
    if [[ -n "${ASSUME}" || -n "${SELECTED}" ]]; then
        if want "$(basename "${script}" _installer.sh)" "${prompt}"; then
            warn "${name}: needs a package path, so it cannot be installed"
            warn "  non-interactively. Later:  bash $(basename "${script}") ${hint}"
        fi
        return
    fi
    if want "$(basename "${script}" _installer.sh)" "${prompt}"; then
        printf "    Path to the %s package (blank to skip): " "${name}"
        read -r reply
        if [[ -n "${reply}" ]]; then
            run_script "${name}" "${script}" "${reply}" || true
        else
            info "Skipping ${name}.  Later:  bash $(basename "${script}") ${hint}"
        fi
    else
        info "Skipping ${name}.  Later:  bash $(basename "${script}") ${hint}"
    fi
}

licensed "SignalP 6" "${SIGNALP_SCRIPT}" "${SIGNALP_PKG}" \
    "Install SignalP 6 locally? It is licensed, so you need the package tarball from DTU. Without it --signalp runs in the cloud." \
    "PACKAGE.tar.gz"

licensed "DeepTMHMM" "${DEEPTMHMM_SCRIPT}" "${DEEPTMHMM_PKG}" \
    "Install DeepTMHMM locally? It is licensed, so you need the package. Without it --tmhmm runs in the cloud." \
    "PACKAGE"

if [[ -f "${THIS_DIR}/tools_table.local.tsv" ]]; then
    info "Installed tool paths were written to tools_table.local.tsv."
    info "  That file is git-ignored; tools_table.tsv keeps the shipped defaults."
fi

printf "\n"
info "Installation complete."
printf "\n"
printf "  Activate the environment:\n"
say "    ${GREEN}conda activate ${ENV_NAME}${NC}"
printf "\n"
printf "  Run FlaGs3 (default — neighbours figure + data tables):\n"
say "    ${GREEN}python FlaGs3.py -i input.txt -u you@example.com -o myrun${NC}"
printf "\n"
printf "  Input list — one query per line, either form, freely mixed:\n"
printf "    WP_047256880.1                            protein only, genome found via NCBI\n"
say "    WP_047256880.1${YELLOW}<TAB>${NC}GCF_000001765.3      protein + NCBI assembly"
say "    MGYG000454827_00001${YELLOW}<TAB>${NC}MGYG000454827   protein + MGnify genome"
say "    ${YELLOW}MGnify genomes need the exact locus tag from that genome's annotation.${NC}"
printf "\n"
printf "  Figures and analysis:\n"
say "    Tree:             ${GREEN}--tree${NC} / ${GREEN}--tree_order${NC}   (mafft + VeryFastTree)"
say "    Domains:          ${GREEN}--domains --hmmdb pfam_db/Pfam-A.hmm${NC}"
say "    Clan colouring:   ${GREEN}--clans pfam_db/Pfam-A.clans.tsv.gz${NC}"
say "    SignalP / DeepTMHMM run in the cloud by default. To run them locally,"
say "    request the licensed packages and then:"
say "      ${GREEN}bash signalp_installer.sh   /path/to/signalp-6-package.tar.gz${NC}"
say "      ${GREEN}bash deeptmhmm_installer.sh /path/to/deeptmhmm-package.tar.gz${NC}"
say "    Each builds its own conda env and fills in tools_table.tsv."
say ""
say "    Defence systems:  ${GREEN}--domains --hmmdb defence=defence_db --hmm_coverage defence=0.3,0.3${NC}"
say "    Membrane/signal:  ${GREEN}--tmhmm${NC} / ${GREEN}--signalp${NC}       (pybiolib, uploads sequences)"
say "    Secretion:        ${GREEN}--sismis${NC}                  (sismis)"
say "    Cluster RNAs:     ${GREEN}--cluster_rna${NC}"
printf "\n"
printf "  Other useful options:\n"
say "    Local genomes:    ${GREEN}--use_local DIR${NC}   search a .gff/.faa directory first"
say "    NCBI API key:     ${GREEN}-api KEY${NC}          raises the download rate cap"
say "    Progress:         ${GREEN}-vb${NC}               per-stage progress and timings"
printf "\n"
say "  Deactivate:         ${GREEN}conda deactivate${NC}"
printf "\n"

exit 0