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

read -rep "Enter full path to source directory (containing MX1 processing output folders): " SOURCE_DIR
read -rep "Enter full path to beamline directory (set as XCE Data Collection Directory):    " DEST_DIR
read -rep "Target name (e.g. cRel, shown in XCE Datasets tab target dropdown):              " TARGET_NAME
read -rep "SMILES library CSV  (e.g. LifeChem...csv) [leave blank to skip]:                 " SMILES_CSV
read -rep "Compound distribution CSV (e.g. MX1-25795...csv) [leave blank to skip]:            " DIST_CSV

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if [[ ! -d "${SOURCE_DIR}" ]]; then
    echo "ERROR: Source directory not found: ${SOURCE_DIR}" >&2
    exit 1
fi

if [[ -z "${TARGET_NAME}" ]]; then
    echo "ERROR: Target name cannot be empty." >&2
    exit 1
fi

PROC_DIR="${DEST_DIR}/processed/${TARGET_NAME}"
mkdir -p "${PROC_DIR}"

echo ""
echo "Source       : ${SOURCE_DIR}"
echo "Beamline dir : ${DEST_DIR}"
echo "Processed dir: ${PROC_DIR}"
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
        log_dir="${PROC_DIR}/${crystal}/${run_num}/mx1_process/output/LogFiles"
        mtz_dir="${PROC_DIR}/${crystal}/${run_num}/mx1_process/output/DataFiles"
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
# Write compound SMILES files (only if both CSVs were provided)
# ---------------------------------------------------------------------------
if [[ -n "${SMILES_CSV:-}" && -n "${DIST_CSV:-}" ]]; then
    echo ""
    echo "Writing compound SMILES files..."
    python3 - "${SMILES_CSV}" "${DIST_CSV}" "${PROC_DIR}" mx1 << 'PYEOF'
import sys, csv, os, re

smiles_csv, dist_csv, dest_dir, mode = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

# SN code → SMILES string
smiles_map = {}
with open(smiles_csv, newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        sn = row.get('CA Sample Number', '').strip()
        smiles = row.get('QCL_SMILES', '').strip()
        if sn and smiles:
            smiles_map[sn] = smiles

# target_directory → SN code → crystal name
crystal_sn = {}  # crystal → SN
with open(dist_csv, newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        sn = row.get('source_directory', '').strip()
        target = row.get('target_directory', '').strip()
        if not sn or not target:
            continue
        fields = target.split('_')
        crystal = '_'.join(fields[:3]) if len(fields) >= 3 else target
        crystal_sn[crystal] = sn

written = 0
missing = 0
for crystal, sn in sorted(crystal_sn.items()):
    smiles = smiles_map.get(sn)
    if not smiles:
        print(f'  WARNING: No SMILES found for {crystal} (SN={sn})')
        missing += 1
        continue
    crystal_dir = os.path.join(dest_dir, crystal)
    if not os.path.isdir(crystal_dir):
        print(f'  SKIP (no dest dir): {crystal}')
        missing += 1
        continue
    smi_path  = os.path.join(crystal_dir, crystal + '.smi')
    cmpd_path = os.path.join(crystal_dir, crystal + '.cmpd')
    with open(smi_path, 'w') as fh:
        fh.write(smiles + '\n')
    with open(cmpd_path, 'w') as fh:
        fh.write(sn + '\n')
    print(f'  SMI  {crystal}  ← {sn}  ({smiles[:30]}...)')
    written += 1

print(f'Wrote {written} SMILES file(s), {missing} missing/skipped.')
PYEOF
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "Prepared ${count} dataset(s), skipped ${skipped}."
echo ""
echo "Next steps in XCE:"
echo "  Settings tab -> Data Collection Directory = ${DEST_DIR}"
echo "  Datasets tab -> Target dropdown            = ${TARGET_NAME}"
echo "  Uncheck 'Read Agamemnon data structure' in Settings tab"
echo "  Datasets tab -> 'Get New Results from Autoprocessing'"
