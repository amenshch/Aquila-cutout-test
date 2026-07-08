#!/usr/bin/env bash
# run_grid.sh — Fill and compute missing models in a 3D parameter grid.
# Run concurrently on multiple processors; uses lock files to avoid double-runs.
# Must be invoked from the top-level directory containing cSD_01 ... cSD_06.

set -euo pipefail
unset -f find 2>/dev/null || true
# Trap Ctrl+C and other signals — release lock and exit cleanly without deleting anything
trap 'echo; echo "Interrupted."; [[ -n "${leaf_dir:-}" ]] && release_lock "$leaf_dir"; exit 130' INT TERM

THREADS=${1:-4}  #<~ 4 is the default value
NPHOT=${2:-}     # optional: ./run_grid.sh 4 10000000000
FORCE=${3:-0}    # set to 1 to rerun all models: ./run_grid.sh 4 10000000000 1

if [[ "$THREADS" == '' ]]
then
  NPROC=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 2)
  THREADS=$(( NPROC < 4 ? NPROC : 4 ))
  echo "Using ${THREADS} threads (${NPROC} processors available)"
fi

DEBUG=0  # set to 0 for normal run
pause() { [[ $DEBUG -eq 1 ]] && read -rp "  Press Enter to continue..." _ || true; }
ZK_CREATED=0  # set to 1 by ensure_zk when it creates a new directory

# ---------------------------------------------------------------------------
# Grid dimensions
# ---------------------------------------------------------------------------
cSD=('cSD_01' 'cSD_02' 'cSD_03' 'cSD_04' 'cSD_05' 'cSD_06')
tBE=('tBE_01' 'tBE_02' 'tBE_03' 'tBE_04' 'tBE_05' 'tBE_06' 'tBE_07' 'tBE_08' 'tBE_09')
imax=6

TOPDIR="$(pwd)"
OUTPUT_FILE="+mc3d.output"
CONFIG_FILE="+radmc3d.inp"
CODE_CMD="mc3d all"
##CODE_CMD="mock_run"
LOCK_SUFFIX=".lock"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

mock_run() {
    sleep 1
    echo "Done..." >> "${OUTPUT_FILE}"
}

leadzeros() { printf "%02d" "$1"; }

get_besden() {
    awk '/^[[:space:]]*besden[[:space:]]*=/ {
        match($0, /[0-9][0-9.]*[eE][+-][0-9]+/)
        print substr($0, RSTART, RLENGTH)
    }' "$1"
}

set_besden() {
    local cfg="$1" val="$2"
    local formatted
    formatted=$(awk -v v="$val" 'BEGIN { printf "%.5e", v }')
    awk -v new="$formatted" '
        /^[[:space:]]*besden[[:space:]]*=/ {
            sub(/[0-9][0-9.]*[eE][+-][0-9]+/, new)
        }
        { print }
    ' "$cfg" > "${cfg}.tmp" && mv "${cfg}.tmp" "$cfg"
}

get_bestem() {
    awk '/^[[:space:]]*bestem[[:space:]]*=/ {
        match($0, /[0-9]+\.[0-9]+/)
        print substr($0, RSTART, RLENGTH)
    }' "$1"
}

set_bestem() {
    local cfg="$1" val="$2"
    local formatted
    formatted=$(awk -v v="$val" 'BEGIN { printf "%.3f", v }')
    awk -v new="$formatted" '
        /^[[:space:]]*bestem[[:space:]]*=/ {
            sub(/[0-9]+\.[0-9]+/, new)
        }
        { print }
    ' "$cfg" > "${cfg}.tmp" && mv "${cfg}.tmp" "$cfg"
}

set_threads() {
    local cfg="$1"
    awk -v new="$THREADS" '
        /^[[:space:]]*setthreads[[:space:]]*=/ {
            sub(/=[[:space:]]*[0-9]+/, "= " new)
        }
        { print }
    ' "$cfg" > "${cfg}.tmp" && mv "${cfg}.tmp" "$cfg"
}

set_nphot() {
    local cfg="$1" val="$2"
    awk -v new="$val" '
        /^[[:space:]]*nphot[[:space:]]*=/ {
            sub(/=[[:space:]]*[0-9]+/, "= " new)
        }
        { print }
    ' "$cfg" > "${cfg}.tmp" && mv "${cfg}.tmp" "$cfg"
}

half_of()   { awk -v v="$1" 'BEGIN { printf "%.10e", v/2.0 }'; }
double_of() { awk -v v="$1" 'BEGIN { printf "%.10f",  v*2.0 }'; }

discover_jmax() {
    local ref="${1:-}"
    if [[ -z "$ref" ]]; then
        echo 0
        return
    fi
    local jmax=0
    for d in "$ref"/tBE_[0-9][0-9]; do
        [[ -d "$d" ]] || continue
        local num; num=$(basename "$d" | sed 's/tBE_//')
        (( 10#$num > jmax )) && jmax=$((10#$num))
    done
    echo "$jmax"
}

discover_kmax() {
    local ref="${1:-}"
    if [[ -z "$ref" ]]; then
        echo 0
        return
    fi
    local kmax=0
    for d in "$ref"/[0-9][0-9]; do
        [[ -d "$d" ]] || continue
        local num; num=$(basename "$d")
        (( 10#$num > kmax )) && kmax=$((10#$num))
    done
    echo "$kmax"
}

# Atomic lock via mkdir; returns 0 on success, 1 if already locked
try_lock()     { mkdir "${1}${LOCK_SUFFIX}" 2>/dev/null; }
release_lock() { rmdir "${1}${LOCK_SUFFIX}" 2>/dev/null || true; }

# ---------------------------------------------------------------------------
# Directory creation
# ---------------------------------------------------------------------------

# Ensure tBE_j exists inside cSD_i; if not, copy from tBE_(j-1) in same cSD_i and multiply bestem by 2

ensure_tbe() {
    unset -f find 2>/dev/null || true
    local csd_dir="$1" j="$2"
    local tbe_dir="${csd_dir}/${tBE[$((j-1))]}"
    [[ -d "$tbe_dir" ]] && return 0
    local prev_tbe="${csd_dir}/${tBE[$((j-2))]}"
    if [[ ! -d "$prev_tbe" ]]; then
        echo "  ERROR: cannot create ${tBE[$((j-1))]} in $(basename $csd_dir): ${tBE[$((j-2))]} missing" >&2
        return 1
    fi
    echo; echo "  [mkdir] $(basename $csd_dir)/${tBE[$((j-1))]} (copied from ${tBE[$((j-2))]}), bestem doubled)"
    cp -a "$prev_tbe" "$tbe_dir"

# Strip outputs and locks from all copied leaves

    for f in "$tbe_dir"/[0-9][0-9]/"+log.MC3D"; do
        rm -f "$f"
    done
    for f in "$tbe_dir"/[0-9][0-9]/"$OUTPUT_FILE"; do
        rm -f "$f"
    done
    for f in "$tbe_dir"/[0-9][0-9]/*"$LOCK_SUFFIX"; do
        rm -f "$f"
    done
    for cfg in "$tbe_dir"/[0-9][0-9]/"$CONFIG_FILE"; do
        [[ -f "$cfg" ]] || continue
        old_val=$(get_bestem "$cfg")
        new_val=$(awk -v v="$old_val" 'BEGIN { printf "%.10f", v*2.0 }')
        'set_bestem' "$cfg" "$new_val" || true
        'set_threads' "$cfg" || true
        [[ -n "$NPHOT" ]] && 'set_nphot' "$cfg" "$NPHOT" || true
    done
}

# Ensure zk leaf exists inside tBE_j/cSD_i; if not, copy from zk-1 and halve besden

ensure_zk() {
    local tbe_dir="$1" k="$2"
    local zk; zk=$(leadzeros "$k")
    local leaf_dir="${tbe_dir}/${zk}"
    ZK_CREATED=0                     # at the top of ensure_zk, before the existence check:
    [[ -d "$leaf_dir" ]] && return 0
    local prev_zk; prev_zk=$(leadzeros $((k-1)))
    local prev_leaf="${tbe_dir}/${prev_zk}"
    if [[ ! -d "$prev_leaf" ]]
    then
        echo "  ERROR: cannot create ${zk} in ${tbe_dir}: ${prev_zk} missing" >&2
        return 1
    fi
    echo; echo "  [mkdir] ${tbe_dir#$TOPDIR/}/${zk} (copied from ${prev_zk}, besden halved)"
    cp -a "$prev_leaf" "$leaf_dir"
    ZK_CREATED=1                     # after the cp -a line:
    rm -f "${leaf_dir}/+log.MC3D"
    rm -f "${leaf_dir}/${OUTPUT_FILE}"
    rm -f "${leaf_dir}/"*"${LOCK_SUFFIX}" 2>/dev/null || true
    local old_val; old_val=$(get_besden "${leaf_dir}/${CONFIG_FILE}")
    'set_besden' "${leaf_dir}/${CONFIG_FILE}" "$(half_of "$old_val")"
    'set_threads' "${leaf_dir}/${CONFIG_FILE}"
    [[ -n "$NPHOT" ]] && 'set_nphot' "${leaf_dir}/${CONFIG_FILE}" "$NPHOT"
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
echo
echo "        Model grid: i=1..${imax}"
echo " Current directory: ${TOPDIR}"
echo

[[ -n "$THREADS"   ]] && echo " NOTE: SETTHREADS updated: ${THREADS}"
[[ -n "$NPHOT"     ]] && echo " NOTE: NPHOT override: ${NPHOT}"
[[ "$FORCE" == "1" ]] && echo " NOTE: FORCED - recompute existing models in place"
if [[ -n "$NPHOT" ]] || [[ "$FORCE" == "1" ]] || [[ "$THREADS" != "4" ]]; then
    echo
    read -rp "  press Enter to continue..."
fi

for (( i=1; i<=imax; i++ ))
do
    csd_dir="${TOPDIR}/${cSD[$((i-1))]}"
    [[ -d "$csd_dir" ]] || { echo " [skip] ${cSD[$((i-1))]} not found"; continue; }
    
    jmax=$(discover_jmax "$csd_dir")
    (( jmax > 0 )) || { echo " [skip] could not determine jmax for ${cSD[$((i-1))]}"; continue; }
    kmax=0

    for (( j=1; j<=jmax; j++ ))
    do
        'ensure_tbe' "$csd_dir" "$j"; rc=$?
        if [[ $rc -ne 0 ]]; then continue; fi

        tbe_dir="${csd_dir}/${tBE[$((j-1))]}"
        kmax=$(discover_kmax "$tbe_dir")

        (( kmax > 0 )) || { echo " [skip] could not determine kmax for ${tBE[$((j-1))]}"; continue; }

        for (( k=1; k<=kmax; k++ ))
        do
            zk=$(leadzeros "$k")
            leaf_dir="${tbe_dir}/${zk}"

            'ensure_zk' "$tbe_dir" "$k" || continue

# Apply nphot to existing dirs if overriding

            [[ -n "$NPHOT" ]] && 'set_nphot' "${leaf_dir}/${CONFIG_FILE}" "$NPHOT"
            [[ "$FORCE" == "1" ]] && 'set_threads' "${leaf_dir}/${CONFIG_FILE}"
            
# Skip if already computed (unless forced)

            if [[ "$FORCE" == "0" ]]; then
                [[ -f "${leaf_dir}/${OUTPUT_FILE}" ]] && grep -q "MC3D: DONE" "+log.MC3D" 2>/dev/null && continue
            fi

# Skip if another process has it

            'try_lock' "$leaf_dir" || { echo "  [skip] ${cSD[$((i-1))]}/${tBE[$((j-1))]}/${zk} — locked"; continue; }

            echo; echo "  [run] ${cSD[$((i-1))]}/${tBE[$((j-1))]}/${zk}"
            pause

# Change directory

            pushd "$leaf_dir" > /dev/null

##            tmux rename-window "$leaf_dir" 2>/dev/null
##            [[ -n "${TMUX:-}" ]] && tmux rename-window -t "${TMUX_PANE:-}" "$leaf_dir" 2>/dev/null || true

# Execute the code

            set +e
            $CODE_CMD
            exit_code=$?
            set -e

            if ! grep -q "MC3D: DONE" "+log.MC3D" 2>/dev/null; then
                echo "  [failed] ${cSD[$((i-1))]}/${tBE[$((j-1))]}/${zk} — mc3d did not complete successfully" >&2
                popd > /dev/null
                'release_lock' "$leaf_dir"
                if [[ $ZK_CREATED -eq 1 ]]; then
                    echo "  [remove] newly created ${leaf_dir#$TOPDIR/} removed"
                    \rm -rf "$leaf_dir"
                fi
                echo
                exit 1
            fi

            echo "  [done] ${cSD[$((i-1))]}/${tBE[$((j-1))]}/${zk}"

            popd > /dev/null

##            tmux rename-window "$PWD" 2>/dev/null
##            [[ -n "${TMUX:-}" ]] && tmux rename-window -t "${TMUX_PANE:-}" "$leaf_dir" 2>/dev/null || true

            'release_lock' "$leaf_dir"

        done  # k

        echo " FINISHED ${cSD[$((i-1))]}/${tBE[$((j-1))]}"
        pause
        
# Before advancing to j+1, check all k models in current tBE_j are done

        all_done=1
        for (( pk=1; pk<=kmax; pk++ )); do
            pzk=$(leadzeros "$pk")
                 
            if ! grep -q "MC3D: DONE" "+log.MC3D" 2>/dev/null; then
                all_done=0
                break
            fi
        done
        if [[ $all_done -eq 0 ]]; then
            echo " [stop] ${cSD[$((i-1))]}/${tBE[$((j-1))]} not fully computed — stopping j loop"
            break
        fi
    done  # j
done  # i

echo
echo "Done."
