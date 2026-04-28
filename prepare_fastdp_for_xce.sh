#!/bin/bash
# prepare_fastdp_for_xce.sh
#
# Reorganises Australian Synchrotron MX3 fast_dp output directories into the
# directory layout required by XChemExplorer (XCE).
#
# Input naming convention (MX3 fast_dp output):
#   fast_dp_results_<crystal>-sn<serial>_<run>/
#     aimless.log
#     fast_dp.mtz
#     ...
#
# Output layout (XCE-compatible, for a specific target):
#   <DEST>/<crystal>/<run>/fast_dp/output/LogFiles/aimless.log
#   <DEST>/<crystal>/<run>/fast_dp/output/DataFiles/<crystal>.free.mtz
#
# XCE traversal path:
#   DEST/<crystal>      ← sample ID
#         /<run>        ← run directory   (first glob level)
#         /fast_dp      ← proc_code       (second glob level)
#         /output       ← subdir          (third glob level, matches toParse "*")
#         /LogFiles/aimless.log           (matches *aimless.log)
#         /DataFiles/<crystal>.free.mtz   (matches *free.mtz)
#
# In XCE Settings tab, set:
#   Data Collection Directory  →  <DEST>
#   Target                     →  (any specific target name, NOT "=== project directory ===")
#
# Usage:
#   bash prepare_fastdp_for_xce.sh
#
#   The script will interactively prompt for:
#     SOURCE_DIR  directory containing fast_dp_results_* subdirectories
#     DEST_DIR    XCE Data Collection Directory to write into

set -euo pipefail

read -rp "Enter full path to source directory (containing fast_dp_results_* folders): " SOURCE_DIR
read -rp "Enter full path to destination directory (XCE Data Collection Directory):    " DEST_DIR
read -rp "SMILES library CSV  (e.g. LifeChem...csv) [leave blank to skip]:              " SMILES_CSV
read -rp "Compound distribution CSV (e.g. MX3...csv) [leave blank to skip]:            " DIST_CSV

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if [[ ! -d "${SOURCE_DIR}" ]]; then
    echo "ERROR: Source directory not found: ${SOURCE_DIR}" >&2
    exit 1
fi

mkdir -p "${DEST_DIR}"

echo "Source : ${SOURCE_DIR}"
echo "Dest   : ${DEST_DIR}"
echo ""

count=0
skipped=0

# ---------------------------------------------------------------------------
# Process each fast_dp results directory
# ---------------------------------------------------------------------------
for src in "${SOURCE_DIR}"/fast_dp_results_*/; do
    [[ -d "${src}" ]] || continue
    src="${src%/}"
    basename="${src##*/}"   # e.g. fast_dp_results_mpc-0020-7-sn02731467_1

    # -----------------------------------------------------------------------
    # Extract crystal name and run number from directory name.
    # Pattern: fast_dp_results_<crystal>-sn<serial>_<run>
    # The crystal name is everything after the prefix up to -sn<digits>_<run>.
    # -----------------------------------------------------------------------
    stripped="${basename#fast_dp_results_}"   # e.g. mpc-0020-7-sn02731467_1

    crystal=$(echo "${stripped}" | sed -E 's/-sn[0-9]+_[0-9]+$//')
    run_num=$(echo "${stripped}"  | sed -E 's/.*_([0-9]+)$/\1/')

    if [[ -z "${crystal}" || "${crystal}" == "${stripped}" ]]; then
        echo "WARNING: Could not parse crystal name from '${basename}' — skipping." >&2
        (( skipped++ )) || true
        continue
    fi

    # -----------------------------------------------------------------------
    # Build the XCE-compatible directory tree (symlinks into original files).
    # Structure: <DEST>/<crystal>/<run>/fast_dp/output/LogFiles|DataFiles/
    # -----------------------------------------------------------------------
    log_dir="${DEST_DIR}/${crystal}/${run_num}/fast_dp/output/LogFiles"
    mtz_dir="${DEST_DIR}/${crystal}/${run_num}/fast_dp/output/DataFiles"
    mkdir -p "${log_dir}" "${mtz_dir}"

    # aimless.log → LogFiles/aimless.log  (matches *aimless.log)
    if [[ -f "${src}/aimless.log" ]]; then
        ln -sf "$(realpath "${src}/aimless.log")" "${log_dir}/aimless.log"
    else
        echo "WARNING: ${crystal} run ${run_num}: aimless.log not found in ${src} — skipping." >&2
        (( skipped++ )) || true
        continue
    fi

    # fast_dp.mtz → DataFiles/<crystal>.free.mtz  (matches *free.mtz)
    if [[ -f "${src}/fast_dp.mtz" ]]; then
        ln -sf "$(realpath "${src}/fast_dp.mtz")" "${mtz_dir}/${crystal}.free.mtz"
    else
        echo "WARNING: ${crystal} run ${run_num}: fast_dp.mtz not found in ${src}." >&2
    fi

    echo "  OK  ${crystal}  (run ${run_num})  ←  ${basename}"
    (( count++ )) || true
done

# ---------------------------------------------------------------------------
# Write compound SMILES files (only if both CSVs were provided)
# ---------------------------------------------------------------------------
if [[ -n "${SMILES_CSV:-}" && -n "${DIST_CSV:-}" ]]; then
    echo ""
    echo "Writing compound SMILES files..."
    python3 - "${SMILES_CSV}" "${DIST_CSV}" "${DEST_DIR}" mx3 << 'PYEOF'
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
        # Strip fast_dp_results_ prefix and -sn<digits>_<digits> suffix
        stripped = re.sub(r'^fast_dp_results_', '', target)
        crystal = re.sub(r'-sn\d+_\d+$', '', stripped)
        if crystal:
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
    smi_path = os.path.join(crystal_dir, crystal + '.smi')
    with open(smi_path, 'w') as fh:
        fh.write(smiles + '\n')
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
echo "  Settings tab → Data Collection Directory = ${DEST_DIR}"
echo "  Settings tab → Target = (your protein target name)"
echo "  Datasets tab → 'Get New Results from Autoprocessing'"
