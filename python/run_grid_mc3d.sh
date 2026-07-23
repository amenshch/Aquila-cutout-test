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
FAILED_LIST="${FAILED_LIST:-${TOPDIR}/failed_nodes.txt}"
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
[[ "${VERBOSE_SKIP:-0}" == "1" ]] && echo "   VERBOSE_SKIP: listing every skipped node"
echo

n_run=0; n_skip=0; n_fail=0
n_done=0; n_undet=0; announced=0

echo " Looking for nodes to compute (set VERBOSE_SKIP=1 to list skipped nodes) ..."


# Read the catalog: skip the header line (starts with ID or #)
while read -r ID i j k SD_emb T_BE M_BE FWHM_pc xi_max contr peak_exc contrast detect hole rho_c rho_edge rho_emb r0 R_out R_cloud N_tot a_BE stab; do
    [[ "$ID" =~ ^[0-9]+$ ]] || continue          # skip header / blanks

    # ---- skip undetectable nodes unless ALL_NODES=1 ----
    #   detect==1 means the model's predicted convolved contrast clears the floor.
    #   Undetectable models never enter a real extraction catalogue, so by default
    #   we do not spend RT time on them; set ALL_NODES=1 to compute the full grid.
    if [[ "${ALL_NODES:-0}" == "0" && "$detect" == "0" ]]; then
        [[ "${VERBOSE_SKIP:-0}" == "1" ]] && \
            printf '  [skip ] cSD_%02d/M_%02d/%02d - undetectable (detect=0)\n' "$i" "$j" "$k"
        n_skip=$((n_skip+1)); n_undet=$((n_undet+1))
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

    # ---- skip if already done AND the stored parameters match the catalog ----
    #   The leaf directory is addressed by (i,j,k) only.  If the catalog is
    #   renumbered or re-gridded, a leaf can already hold a COMPLETED model that
    #   was computed with different physics.  write_config above has just
    #   overwritten the config with the new values, so the config would silently
    #   disagree with the images.  Guard against that: stash the parameters that
    #   produced the output and recompute whenever they differ.
    stamp="${leaf_dir}/.mc3d_params"
    want=$(printf '%.6e %.6f %.4f %.6e' "$rho_c" "$T_BE" "$xi_max" "$SD_emb")
    if [[ "$FORCE" == "0" ]]; then
        if [[ -f "${leaf_dir}/${OUTPUT_FILE}" ]] && \
           grep -q "MC3D: DONE" "${leaf_dir}/${LOG_FILE}" 2>/dev/null; then
            have=$(cat "$stamp" 2>/dev/null || echo "")
            if [[ "$have" == "$want" ]]; then
                [[ "${VERBOSE_SKIP:-0}" == "1" ]] && \
                    echo "  [skip ] ${csd}/${md}/${kd} - already done"
                n_skip=$((n_skip+1)); n_done=$((n_done+1))
                continue
            fi
            if [[ -z "$have" ]]; then
                # Legacy leaf from before the guard existed: adopt it, but say so.
                echo "  [adopt] ${csd}/${md}/${kd} - no parameter stamp; assuming it matches"
                printf '%s' "$want" > "$stamp"
                n_skip=$((n_skip+1))
                continue
            fi
            echo "  [STALE] ${csd}/${md}/${kd} - parameters changed, recomputing" >&2
            echo "          had:  $have" >&2
            echo "          want: $want" >&2
        fi
    fi

    # ---- previously failed? skip unless explicitly retrying ----
    if [[ -f "${leaf_dir}/.mc3d_failed" && "${RETRY_FAILED:-0}" == "0" && "$FORCE" == "0" ]]; then
        echo "  [skip] ${csd}/${md}/${kd} — previously FAILED (set RETRY_FAILED=1 to retry)"
        n_skip=$((n_skip+1))
        continue
    fi

    # ---- lock (allows several copies of this script to run in parallel) ----
    try_lock "$leaf_dir" || { echo "  [skip] ${csd}/${md}/${kd} — locked"; continue; }

    if [[ $announced -eq 0 ]]; then
        printf ' Skipped %d node(s): %d already done, %d undetectable.\n' \
               "$n_skip" "$n_done" "$n_undet"
        echo " Computing:"
        announced=1
    fi

    rm -f "${leaf_dir}/.mc3d_failed"

    # ---- write the parameters, now that we know we are actually running ----
    #   This used to run before the skip tests, which meant six awk/grep/mv
    #   passes over the config of every already-completed node -- ~18 forks per
    #   node, silently, before the loop reached anything that needed computing.
    #   With the refined 742-node catalog that preamble became a noticeable
    #   delay.  A node that is skipped keeps the config written when it ran,
    #   which by the stamp test already matches the catalog.
    write_config "$leaf_dir" "$rho_c" "$T_BE" "$xi_max" "$SD_emb"

    echo "  [run ] ${csd}/${md}/${kd}   T_BE=${T_BE}  xi=${xi_max}  M=${M_BE}"
    pushd "$leaf_dir" > /dev/null
    set +e
    $CODE_CMD
    set -e
    popd > /dev/null

    if grep -q "MC3D: DONE" "${leaf_dir}/${LOG_FILE}" 2>/dev/null; then
        echo "  [done] ${csd}/${md}/${kd}"
        printf '%s' "$want" > "${leaf_dir}/.mc3d_params"
        n_run=$((n_run+1))
    else
        echo "  [FAIL] ${csd}/${md}/${kd} — mc3d did not complete" >&2
        n_fail=$((n_fail+1))
        # Do NOT delete the leaf: the log is the only record of why it failed,
        # and deleting it also hides the failure from a later re-run, which
        # would simply recreate the directory and fail again in the same way.
        # Mark it instead, and preserve the log under a distinct name.
        {
            echo "failed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
            echo "params: $want"
            echo "node:   i=${i} j=${j} k=${k}  M_BE=${M_BE}  xi_max=${xi_max}"
        } > "${leaf_dir}/.mc3d_failed"
        [[ -f "${leaf_dir}/${LOG_FILE}" ]] && \
            cp -f "${leaf_dir}/${LOG_FILE}" "${leaf_dir}/mc3d.failed.log" 2>/dev/null || true
        # record for the optional retry pass
        echo "${i} ${j} ${k}" >> "${FAILED_LIST}"
    fi

    release_lock "$leaf_dir"
    leaf_dir=""
done < "$GRID_CATALOG"

echo
echo "Done.  run=${n_run}  skipped=${n_skip}  failed=${n_fail}   (catalog: ${NMODELS} models)"
