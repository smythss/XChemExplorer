#!/bin/bash
# prepare_mx1_for_xce.sh
#
# Reorganises Australian Synchrotron MX1 autoprocessing output directories
# into the directory layout required by XChemExplorer (XCE).
#
# Input naming convention (MX1 autoprocessing output):
#   <CRYSTAL>_<TIMESTAMP>_<TYPE>/
#     aimless.log
#     <CRYSTAL>_<TIMESTAMP>_<TYPE>_aimless.mtz
#     ...
#   e.g. MPC0022_MB6_0084_20260313-163230_process/
#        MPC0022_MB6_0084_20260313-205140_retrigger/
#
# Crystal name is extracted from the first 3 underscore-delimited fields:
#   MPC0022_MB6_0084
#
# Multiple processing directories per crystal (process, retrigger, etc.) are
# each assigned an incrementing run number (1, 2, ...) ordered by timestamp.
#
# Output layout (XCE-compatible, for a specific target):
#   <DEST>/<crystal>/<run>/mx1_process/output/LogFiles/aimless.log
#   <DEST>/<crystal>/<run>/mx1_process/output/DataFiles/<crystal>.free.mtz
#
# XCE traversal path:
#   DEST/<crystal>           ← sample ID
#         /<run>             ← run directory   (first glob level)
#         /mx1_process       ← proc_code       (second glob level)
#         /output            ← subdir          (third glob level, matches "*")
#         /LogFiles/aimless.log                (matches *aimless.log)
#         /DataFiles/<crystal>.free.mtz        (matches *free.mtz)
#
# In XCE Settings tab, set:
#   Data Collection Directory  →  <DEST>
#   Target                     →  (any specific target name,
#                                  NOT "=== project directory ===")
#
# Usage:
#   bash prepare_mx1_for_xce.sh
#   (The script will prompt for full paths interactively.)

set -euo pipefail

read -rp "Enter full path to source directory (containing MX1 processing output folders): " SOURCE_DIR
read -rp "Enter full path to destination directory (XCE Data Collection Directory):        " DEST_DIR

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if [[ ! -d "${SOURCE_DIR}" ]]; then
    echo "ERROR: Source directory not found: ${SOURCE_DIR}" >&2
    exit 1
fi

mkdir -p "${DEST_DIR}"

echo ""
echo "Source : ${SOURCE_DIR}"
echo "Dest   : ${DEST_DIR}"
echo ""

count=0
skipped=0

# ---------------------------------------------------------------------------
# Collect and sort all processing directories, then group by crystal name.
# Sorting by directory name puts them in timestamp order naturally.
# ---------------------------------------------------------------------------

# Associative array: crystal_name → space-separated list of sorted src dirs
declare -A crystal_dirs

for src in $(ls -d "${SOURCE_DIR}"/*/  2>/dev/null | sort); do
    src="${src%/}"
    basename="${src##*/}"

    # Must contain an underscore-separated structure with at least 4 fields
    # and end in _process or _retrigger (or any _<type>)
    # Crystal name = first 3 underscore-delimited fields
    IFS='_' read -ra fields <<< "${basename}"
    if [[ ${#fields[@]} -lt 4 ]]; then
        continue
    fi

    crystal="${fields[0]}_${fields[1]}_${fields[2]}"

    if [[ -z "${crystal_dirs[${crystal}]+_}" ]]; then
        crystal_dirs["${crystal}"]="${src}"
    else
        crystal_dirs["${crystal}"]+=" ${src}"
    fi
done

# ---------------------------------------------------------------------------
# Process each crystal and its associated run directories
# ---------------------------------------------------------------------------
for crystal in $(echo "${!crystal_dirs[@]}" | tr ' ' '\n' | sort); do
    run_num=1

    for src in ${crystal_dirs["${crystal}"]}; do
        basename="${src##*/}"

        # Skip if no aimless.log (processing failed for this dataset)
        if [[ ! -f "${src}/aimless.log" ]]; then
            echo "  SKIP ${crystal} (${basename}) — no aimless.log found"
            (( skipped++ )) || true
            continue
        fi

        # Find the named aimless MTZ: <basename>_aimless.mtz
        mtz_file="${src}/${basename}_aimless.mtz"
        if [[ ! -f "${mtz_file}" ]]; then
            # Fallback: any *_aimless.mtz in the directory
            mtz_file=$(ls "${src}"/*_aimless.mtz 2>/dev/null | head -1 || true)
        fi

        if [[ -z "${mtz_file}" || ! -f "${mtz_file}" ]]; then
            echo "  SKIP ${crystal} run ${run_num} (${basename}) — no *_aimless.mtz found"
            (( skipped++ )) || true
            (( run_num++ )) || true
            continue
        fi

        # -----------------------------------------------------------------------
        # Build XCE-compatible directory tree (symlinks into original files)
        # <DEST>/<crystal>/<run>/mx1_process/output/LogFiles|DataFiles/
        # -----------------------------------------------------------------------
        log_dir="${DEST_DIR}/${crystal}/${run_num}/mx1_process/output/LogFiles"
        mtz_dir="${DEST_DIR}/${crystal}/${run_num}/mx1_process/output/DataFiles"
        mkdir -p "${log_dir}" "${mtz_dir}"

        # aimless.log → LogFiles/aimless.log   (matches *aimless.log)
        ln -sf "$(realpath "${src}/aimless.log")" "${log_dir}/aimless.log"

        # *_aimless.mtz → DataFiles/<crystal>.free.mtz  (matches *free.mtz)
        ln -sf "$(realpath "${mtz_file}")" "${mtz_dir}/${crystal}.free.mtz"

        echo "  OK  ${crystal}  run ${run_num}  ←  ${basename}"
        (( count++ )) || true
        (( run_num++ )) || true
    done
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "Prepared ${count} dataset(s), skipped ${skipped}."
echo ""
echo "Next steps in XCE:"
echo "  Settings tab → Data Collection Directory = ${DEST_DIR}"
echo "  Settings tab → Target = (your protein target name)"
echo "  Datasets tab → 'Get New Results from Autoprocessing'"
