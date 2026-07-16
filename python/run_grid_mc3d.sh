#!/usr/bin/env bash
# run_grid_mc3d.sh -- Fill and compute the BE model grid from an explicit catalog.
#
# CHANGED vs the previous version:
#   The old script GENERATED each model by copying its neighbour and scaling
#   (bestem doubled, besden halved).  That only works for a regular geometric
#   lattice in (SD_emb, T_BE, rho_c).  In the pressure-truncated grid T_BE,
#   rho_c and xi_max are all DERIVED per node and take arbitrary values, so
#   copy-and-scale cannot reproduce them.  This version instead reads every
#   model's parameters from the grid catalog and writes them into the config.
#
#   Axes (directory tree, unchanged in shape):
#       i -> cSD_%02d   embedding column  SD_emb   (6 values)
#       j -> M_%02d     core mass         M_BE     (12 values)
#       k -> %02d       core size         FWHM     (10 values)
#   The tree is ragged: only 277 of 6x12x10 slots are physical; the catalog
#   lists exactly those, and only those are created.
#
# Usage:
#   ./run_grid_mc3d.sh [threads] [nphot] [force]
#     ./run_grid_mc3d.sh 2                 # 2 threads, config nphot, skip done
#     ./run_grid_mc3d.sh 2 1000000000      # override nphot
#     ./run_grid_mc3d.sh 2 1000000000 1    # force recompute of existing models
#
# Must be run from the top-level grid directory.  Needs:
#   - GRID_CATALOG  (default ./bes_grid_final.txt)
#   - TEMPLATE_DIR  (default ./+template) containing a reference +radmc3d.inp
#     and whatever static files a model needs (opacities, wavelengths, ISRF).

set -euo pipefail

THREADS=${1:-2}
NPHOT=${2:-}
FORCE=${3:-0}

GRID_CATALOG="${GRID_CATALOG:-./bes_grid_final.txt}"
TEMPLATE_DIR="${TEMPLATE_DIR:-./+template}"
OUTPUT_FILE="+mc3d.output"
CONFIG_FILE="+radmc3d.inp"
LOG_FILE="+log.MC3D"
CODE_CMD="mc3d all"
LOCK_SUFFIX=".lock"

TOPDIR="$(pwd)"
leaf_dir=""
trap 'echo; echo "Interrupted."; [[ -n "${leaf_dir:-}" ]] && release_lock "$leaf_dir"; exit 130' INT TERM

try_lock()     { mkdir "${1}${LOCK_SUFFIX}" 2>/dev/null; }
release_lock() { rmdir "${1}${LOCK_SUFFIX}" 2>/dev/null || true; }

# ---------------------------------------------------------------------------
# set_param <config> <key> <value>
#   Replaces the RHS of "key = <old>" with <value>, preserving the rest of the
#   line.  Works for both float and integer keys.  Adjust the key NAMES in
#   write_config() below to match your +radmc3d.inp.
# ---------------------------------------------------------------------------
set_param() {
    local cfg="$1" key="$2" val="$3"
    if ! grep -qE "^[[:space:]]*${key}[[:space:]]*=" "$cfg"; then
        echo "  ERROR: key '${key}' not found in ${cfg}" >&2
        return 1
    fi
    # replace only the VALUE, preserving any trailing "! comment"
    awk -v k="$key" -v v="$val" '
        $0 ~ "^[[:space:]]*"k"[[:space:]]*=" {
            comment = ""
            if (match($0, /![^\n]*$/)) comment = substr($0, RSTART)
            printf " %s = %s", k, v
            if (comment != "") printf "  %s", comment
            printf "\n"
            next
        }
        { print }
    ' "$cfg" > "${cfg}.tmp" && mv "${cfg}.tmp" "$cfg"
}

# ---------------------------------------------------------------------------
# write_config <leaf_dir> <rho_c> <T_BE> <xi_max> <SD_emb>
#
#   Only FOUR physical parameters are written; everything else the Fortran
#   derives for itself:
#       besden  <- rho_c    central BE density            [g/cm^3]
#       bestem  <- T_BE     construction temperature      [K]
#       ximax   <- xi_max   BE truncation (6.451=critical)
#       embden  <- SD_emb   embedding cloud SURFACE density [cm^-2]
#
#   NOT set here (deliberately):
#       Rbesph  is a FLAG (1.0 = BE sphere present), not a radius.  Leave at 1.0.
#       Router  = 0.0 -> chosen automatically by the code (now 30*R_BE).
#       Rinner  = 0.0 -> automatic.
#   R_out = ximax * r0 is derived inside the code from besden/bestem/ximax,
#   and rho_emb is derived from embden and the geometry.  The catalog columns
#   R_out_AU, R_cloud_AU and rho_emb are therefore cross-checks, not inputs:
#   compare them against the code's own values on the first few models.
# ---------------------------------------------------------------------------
write_config() {
    local dir="$1" rho_c="$2" tbe="$3" ximax="$4" sdemb="$5"
    local cfg="${dir}/${CONFIG_FILE}"

    set_param "$cfg" besden "$(printf '%.6e' "$rho_c")"
    set_param "$cfg" bestem "$(printf '%.6f' "$tbe")"
    set_param "$cfg" ximax  "$(printf '%.4f' "$ximax")"
    set_param "$cfg" embden "$(printf '%.6e' "$sdemb")"

    set_param "$cfg" setthreads "$THREADS" || true
    [[ -n "$NPHOT" ]] && set_param "$cfg" nphot "$NPHOT" || true
}

# ---------------------------------------------------------------------------
if [[ ! -f "$GRID_CATALOG" ]]; then
    echo " ERROR: grid catalog not found: $GRID_CATALOG" >&2; exit 1
fi
if [[ ! -d "$TEMPLATE_DIR" ]]; then
    echo " ERROR: template dir not found: $TEMPLATE_DIR" >&2; exit 1
fi

NMODELS=$(grep -vc '^[[:space:]]*#' "$GRID_CATALOG" | head -1)
echo
echo "   Grid catalog: ${GRID_CATALOG}"
echo "       Template: ${TEMPLATE_DIR}"
echo "        Threads: ${THREADS}"
[[ -n "$NPHOT"     ]] && echo "          NPHOT: ${NPHOT} (override)"
[[ "$FORCE" == "1" ]] && echo "          FORCE: recompute existing models"
echo

n_run=0; n_skip=0; n_fail=0

# Read the catalog: skip the header line (starts with ID or #)
while read -r ID i j k SD_emb T_BE M_BE FWHM_pc xi_max contr peak_exc contrast detect hole rho_c rho_edge rho_emb r0 R_out R_cloud N_tot a_BE stab; do
    [[ "$ID" =~ ^[0-9]+$ ]] || continue          # skip header / blanks

    # ---- skip undetectable nodes unless ALL_NODES=1 ----
    #   detect==1 means the model's predicted convolved contrast clears the floor.
    #   Undetectable models never enter a real extraction catalogue, so by default
    #   we do not spend RT time on them; set ALL_NODES=1 to compute the full grid.
    if [[ "${ALL_NODES:-0}" == "0" && "$detect" == "0" ]]; then
        n_skip=$((n_skip+1))
        continue
    fi

    csd=$(printf 'cSD_%02d' "$i")
    md=$(printf  'M_%02d'   "$j")
    kd=$(printf  '%02d'     "$k")
    leaf_dir="${TOPDIR}/${csd}/${md}/${kd}"

    # ---- create the leaf from the template if it does not exist ----
    created=0
    if [[ ! -d "$leaf_dir" ]]; then
        mkdir -p "$(dirname "$leaf_dir")"
        cp -a "$TEMPLATE_DIR" "$leaf_dir"
        rm -f "${leaf_dir}/${LOG_FILE}" "${leaf_dir}/${OUTPUT_FILE}"
        rm -rf "${leaf_dir}"*"${LOCK_SUFFIX}" 2>/dev/null || true
        created=1
        echo "  [mkdir] ${csd}/${md}/${kd}"
    fi

    # ---- always (re)write the parameters: they come from the catalog, not
    #      from a neighbouring model, so this is cheap and idempotent ----
    write_config "$leaf_dir" "$rho_c" "$T_BE" "$xi_max" "$SD_emb"

    # ---- skip if already done ----
    if [[ "$FORCE" == "0" ]]; then
        if [[ -f "${leaf_dir}/${OUTPUT_FILE}" ]] && \
           grep -q "MC3D: DONE" "${leaf_dir}/${LOG_FILE}" 2>/dev/null; then
            n_skip=$((n_skip+1))
            continue
        fi
    fi

    # ---- lock (allows several copies of this script to run in parallel) ----
    try_lock "$leaf_dir" || { echo "  [skip] ${csd}/${md}/${kd} — locked"; continue; }

    echo "  [run ] ${csd}/${md}/${kd}   T_BE=${T_BE}  xi=${xi_max}  M=${M_BE}"
    pushd "$leaf_dir" > /dev/null
    set +e
    $CODE_CMD
    set -e
    popd > /dev/null

    if grep -q "MC3D: DONE" "${leaf_dir}/${LOG_FILE}" 2>/dev/null; then
        echo "  [done] ${csd}/${md}/${kd}"
        n_run=$((n_run+1))
    else
        echo "  [FAIL] ${csd}/${md}/${kd} — mc3d did not complete" >&2
        n_fail=$((n_fail+1))
        [[ $created -eq 1 ]] && { echo "  [remove] ${csd}/${md}/${kd}"; rm -rf "$leaf_dir"; }
    fi

    release_lock "$leaf_dir"
    leaf_dir=""
done < "$GRID_CATALOG"

echo
echo "Done.  run=${n_run}  skipped=${n_skip}  failed=${n_fail}   (catalog: ${NMODELS} models)"
