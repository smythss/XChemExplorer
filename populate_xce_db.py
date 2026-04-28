#!/usr/bin/env python3
"""
populate_xce_db.py  —  Directly populate an XChemExplorer SQLite database by
walking a processed/<target>/<crystal>/<run>/<proc_code>/output/ directory tree
and parsing aimless.log files.

Bypasses XCE's built-in autoprocessing scan, which has hard-coded DLS path
assumptions. Works with Australian Synchrotron MX1 and MX3 (fast_dp) data.

Prerequisite — create the SQLite file first via XCE:
    1. Launch XCE and set the Project Directory (Settings tab).
    2. XCE creates <project_dir>/<project_name>.sqlite automatically on startup.
    3. Close XCE (or keep it closed while this script runs).

Usage:
    python3 populate_xce_db.py

    The script will prompt for each required path interactively.
    Tab completion is supported.

Compound / SMILES data:
    If you ran prepare_fastdp_for_xce.sh or prepare_mx1_for_xce.sh first,
    compound information is already stored as <crystal>.smi and <crystal>.cmpd
    files inside each crystal directory. This script reads those files directly
    — you do NOT need to supply the library or distribution CSVs.

    The CSV prompts below are only a fallback for cases where the prep scripts
    were not run (or were run without SMILES data). Leave them blank if the
    .smi / .cmpd files are present.

After running this script:
    1. Open XCE and go to the Datasets tab.
    2. Select your target from the dropdown.
    3. Click "Select Best Autoprocessing Result" (Run button, NOT Status) to
       trigger auto-assignment and populate DataProcessingAutoAssigned.
    4. Do NOT click the Status buttons — they crash at WEHI (CLUSTER_BASTION
       is a DLS-only constant that is not defined here).
"""

import csv
import glob
import math
import os
import re
import readline
import shutil
import sqlite3
from datetime import datetime


# ---------------------------------------------------------------------------
# Space group → lattice + point group (matches XChemUtils.point_group_dict)
# ---------------------------------------------------------------------------

# Maps point group → list of space group strings (no spaces)
_POINT_GROUP_DICT = {
    "P 1":       ["P1"],
    "P 2":       ["P2", "P21", "C2"],
    "P 2 2 2":   ["P222", "P2221", "P21212", "P212121",
                  "C222", "C2221", "F222", "I222", "I212121"],
    "P 4":       ["P4", "P41", "P42", "P43", "I4", "I41"],
    "P 4 2 2":   ["P422", "P4212", "P4122", "P41212", "P4222", "P42212",
                  "P4322", "P43212", "I422", "I4122"],
    "P 3":       ["P3", "P31", "P32", "R3"],
    "P 3 1 2":   ["P312", "P3112", "P3212", "P321", "P3121", "P3221",
                  "R32"],
    "P 6":       ["P6", "P61", "P65", "P62", "P64", "P63"],
    "P 6 2 2":   ["P622", "P6122", "P6522", "P6222", "P6422", "P6322"],
    "P 2 3":     ["P23", "F23", "I23", "P213", "I213"],
    "P 4 3 2":   ["P432", "P4232", "F432", "F4132", "I432",
                  "P4332", "P4132", "I4132"],
}

_LATTICE_MAP = {
    "P1": "triclinic",
    "P2": "monoclinic (primitive)", "P21": "monoclinic (primitive)",
    "C2": "monoclinic (centred)",
    "P222": "orthorhombic", "P2221": "orthorhombic", "P21212": "orthorhombic",
    "P212121": "orthorhombic", "C222": "orthorhombic", "C2221": "orthorhombic",
    "F222": "orthorhombic", "I222": "orthorhombic", "I212121": "orthorhombic",
    "P4": "tetragonal", "P41": "tetragonal", "P42": "tetragonal",
    "P43": "tetragonal", "I4": "tetragonal", "I41": "tetragonal",
    "P422": "tetragonal", "P4212": "tetragonal", "P4122": "tetragonal",
    "P41212": "tetragonal", "P4222": "tetragonal", "P42212": "tetragonal",
    "P4322": "tetragonal", "P43212": "tetragonal", "I422": "tetragonal",
    "I4122": "tetragonal",
    "P3": "trigonal", "P31": "trigonal", "P32": "trigonal",
    "R3": "rhombohedral", "P312": "trigonal", "P3112": "trigonal",
    "P3212": "trigonal", "P321": "trigonal", "P3121": "trigonal",
    "P3221": "trigonal", "R32": "rhombohedral",
    "P6": "hexagonal", "P61": "hexagonal", "P65": "hexagonal",
    "P62": "hexagonal", "P64": "hexagonal", "P63": "hexagonal",
    "P622": "hexagonal", "P6122": "hexagonal", "P6522": "hexagonal",
    "P6222": "hexagonal", "P6422": "hexagonal", "P6322": "hexagonal",
    "P23": "cubic", "F23": "cubic", "I23": "cubic", "P213": "cubic",
    "I213": "cubic", "P432": "cubic", "P4232": "cubic", "F432": "cubic",
    "F4132": "cubic", "I432": "cubic", "P4332": "cubic", "P4132": "cubic",
    "I4132": "cubic",
}

# Reverse map: space group (no spaces) → point group
_SG_TO_PG = {}
for _pg, _sgs in _POINT_GROUP_DICT.items():
    for _sg in _sgs:
        _SG_TO_PG[_sg] = _pg

# Point group → ASU count (for scoring)
_ASU_COUNT = {
    "P 1": 1, "P 2": 2, "P 2 2 2": 4,
    "P 4": 4, "P 4 2 2": 8,
    "P 3": 3, "P 3 1 2": 6,
    "P 6": 6, "P 6 2 2": 12,
    "P 2 3": 12, "P 4 3 2": 24,
}


def sg_to_lattice(sg_no_spaces):
    return _LATTICE_MAP.get(sg_no_spaces, "n/a")


def sg_to_pointgroup(sg_no_spaces):
    return _SG_TO_PG.get(sg_no_spaces, "n/a")


# ---------------------------------------------------------------------------
# Aimless log parser  (mirrors XChemUtils.parse().aimless_logile exactly)
# ---------------------------------------------------------------------------

def parse_aimless_log(logfile_path):
    """Parse an aimless.log and return a dict of DataProcessing* fields."""
    d = {k: "n/a" for k in [
        "DataCollectionWavelength",
        "DataProcessingResolutionLow", "DataProcessingResolutionHigh",
        "DataProcessingResolutionHighOuterShell", "DataProcessingResolutionLowInnerShell",
        "DataProcessingResolutionHigh15sigma", "DataProcessingResolutionHigh20sigma",
        "DataProcessingRmergeOverall", "DataProcessingRmergeLow", "DataProcessingRmergeHigh",
        "DataProcessingIsigOverall", "DataProcessingIsigLow", "DataProcessingIsigHigh",
        "DataProcessingCompletenessOverall", "DataProcessingCompletenessLow", "DataProcessingCompletenessHigh",
        "DataProcessingMultiplicityOverall", "DataProcessingMultiplicityLow", "DataProcessingMultiplicityHigh",
        "DataProcessingCChalfOverall", "DataProcessingCChalfLow", "DataProcessingCChalfHigh",
        "DataProcessingUniqueReflectionsOverall", "DataProcessingUniqueReflectionsLow", "DataProcessingUniqueReflectionsHigh",
        "DataProcessingSpaceGroup", "DataProcessingLattice", "DataProcessingPointGroup",
        "DataProcessingA", "DataProcessingB", "DataProcessingC",
        "DataProcessingAlpha", "DataProcessingBeta", "DataProcessingGamma",
    ]}

    res15_section = False
    res20_section = False

    try:
        with open(logfile_path) as f:
            for line in f:
                p = line.split()

                if "Wavelength" in line and len(p) >= 2:
                    d["DataCollectionWavelength"] = p[1]

                if "Low resolution limit" in line and len(p) == 6:
                    d["DataProcessingResolutionLow"] = p[3]
                    d["DataProcessingResolutionHighOuterShell"] = p[5]

                if "High resolution limit" in line and len(p) == 6:
                    d["DataProcessingResolutionHigh"] = p[3]
                    d["DataProcessingResolutionLowInnerShell"] = p[4]

                if ("Rmerge  (all I+ and I-)" in line or "Rmerge  (all I+ & I-)" in line) and len(p) == 8:
                    d["DataProcessingRmergeOverall"] = p[5]
                    d["DataProcessingRmergeLow"] = p[6]
                    d["DataProcessingRmergeHigh"] = p[7]

                if ("Mean((I)/sd(I))" in line or "Mean(I)/sd(I)" in line) and len(p) == 4:
                    d["DataProcessingIsigOverall"] = p[1]
                    d["DataProcessingIsigLow"] = p[2]
                    d["DataProcessingIsigHigh"] = p[3]

                if line.startswith("Completeness") and len(p) == 4:
                    d["DataProcessingCompletenessOverall"] = p[1]
                    d["DataProcessingCompletenessLow"] = p[2]
                    d["DataProcessingCompletenessHigh"] = p[3]

                if "Completeness (ellipsoidal)" in line and len(p) == 5:
                    d["DataProcessingCompletenessOverall"] = p[2]
                    d["DataProcessingCompletenessLow"] = p[3]
                    d["DataProcessingCompletenessHigh"] = p[4]

                if "Multiplicity" in line and len(p) == 4:
                    d["DataProcessingMultiplicityOverall"] = p[1]
                    d["DataProcessingMultiplicityLow"] = p[2]
                    d["DataProcessingMultiplicityHigh"] = p[3]

                if line.startswith("Mn(I) half-set correlation CC(1/2)") and len(p) == 7:
                    d["DataProcessingCChalfOverall"] = p[4]
                    d["DataProcessingCChalfLow"] = p[5]
                    d["DataProcessingCChalfHigh"] = p[6]

                if line.startswith("     CC(1/2)") and len(p) == 4:
                    d["DataProcessingCChalfOverall"] = p[1]
                    d["DataProcessingCChalfLow"] = p[2]
                    d["DataProcessingCChalfHigh"] = p[3]

                if "Total number unique" in line and len(p) == 6:
                    d["DataProcessingUniqueReflectionsOverall"] = p[3]

                if (line.startswith("Average unit cell:") or
                        line.startswith("  Unit cell parameters")) and len(p) == 9:
                    d["DataProcessingA"] = str(int(float(p[3])))
                    d["DataProcessingB"] = str(int(float(p[4])))
                    d["DataProcessingC"] = str(int(float(p[5])))
                    d["DataProcessingAlpha"] = str(int(float(p[6])))
                    d["DataProcessingBeta"] = str(int(float(p[7])))
                    d["DataProcessingGamma"] = str(int(float(p[8])))

                if (line.startswith("Space group:") or line.startswith("  Spacegroup name")):
                    if "Laue" in line:
                        continue
                    if "Spacegroup name" in line:
                        sg_raw = line.replace("  Spacegroup name", "").strip().rstrip()
                        sg = sg_raw.replace(" ", "")
                    else:
                        sg_raw = line.replace("Space group: ", "").rstrip()
                        sg = sg_raw.strip()
                    d["DataProcessingSpaceGroup"] = sg_raw.strip()
                    d["DataProcessingLattice"] = sg_to_lattice(sg)
                    d["DataProcessingPointGroup"] = sg_to_pointgroup(sg)

                if line.startswith("Estimates of resolution limits: overall"):
                    res15_section = True
                    res20_section = True

                if res15_section and "from Mn(I/sd)" in line and len(p) >= 7:
                    if "1.5" in p[3]:
                        d["DataProcessingResolutionHigh15sigma"] = p[6].rstrip(",")
                        res15_section = False

                if res20_section and "from Mn(I/sd)" in line and len(p) >= 7:
                    if "2.0" in p[3]:
                        d["DataProcessingResolutionHigh20sigma"] = p[6].rstrip(",")
                        res20_section = False

    except Exception as e:
        print(f"  WARNING: could not parse {logfile_path}: {e}")

    return d


def calc_unit_cell_volume(a, b, c, alpha_deg, beta_deg, gamma_deg, lattice):
    """Matches XChemUtils.calc_unitcell_volume_from_logfile (angles in degrees)."""
    al = math.radians(alpha_deg)
    be = math.radians(beta_deg)
    ga = math.radians(gamma_deg)
    if lattice == "triclinic":
        return a * b * c * math.sqrt(
            1 - math.cos(al)**2 - math.cos(be)**2 - math.cos(ga)**2
            + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    if "monoclinic" in lattice:
        return round(a * b * c * math.sin(be), 1)
    if lattice in ("orthorhombic", "tetragonal", "cubic"):
        return round(a * b * c, 1)
    if lattice in ("hexagonal", "rhombohedral"):
        return round(a * b * c * math.sin(math.radians(60)), 1)
    return a * b * c  # fallback


def compute_score_and_volume(d):
    try:
        a = float(d["DataProcessingA"])
        b = float(d["DataProcessingB"])
        c = float(d["DataProcessingC"])
        al = float(d["DataProcessingAlpha"])
        be = float(d["DataProcessingBeta"])
        ga = float(d["DataProcessingGamma"])
        vol = calc_unit_cell_volume(a, b, c, al, be, ga, d["DataProcessingLattice"])
        boost = _ASU_COUNT.get(d["DataProcessingPointGroup"], 1)
        score = (float(d["DataProcessingUniqueReflectionsOverall"])
                 * float(d["DataProcessingCompletenessOverall"])
                 * boost
                 * float(d["DataProcessingIsigOverall"])) / vol
        return str(round(score, 3)), str(round(vol, 1))
    except (ValueError, KeyError, ZeroDivisionError, TypeError):
        return "0.0", "n/a"


def compute_alert(d):
    try:
        res = float(d["DataProcessingResolutionHigh"])
        rm = float(d["DataProcessingRmergeLow"])
        if res > 3.5 or rm > 0.1:
            return "#FF0000"
        if 2.5 < res <= 3.5 or 0.05 < rm <= 0.1:
            return "#FF9900"
        return "#00FF00"
    except (ValueError, TypeError):
        return "#FF0000"


# ---------------------------------------------------------------------------
# Compound / SMILES helpers
# ---------------------------------------------------------------------------

def load_smiles_map(smiles_csv):
    """Return dict: SN_code -> SMILES  from LifeChem-style library CSV.

    Expects columns 'CA Sample Number' and 'QCL_SMILES'.
    """
    smap = {}
    with open(smiles_csv, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sn     = row.get("CA Sample Number", "").strip()
            smiles = row.get("QCL_SMILES", "").strip()
            if sn and smiles:
                smap[sn] = smiles
    return smap


def load_crystal_sn_map(dist_csv, beamline):
    """Return dict: crystal_name -> SN_code from a distribution CSV.

    Supports two formats:
      MX3/ASMX  – columns 'Name' (e.g. MPC-0019-1-SN02731452) and 'Compound' (SN code)
      MX1       – columns 'source_directory' (SN) and 'target_directory' (dir/crystal name)
    """
    crystal_sn = {}
    with open(dist_csv, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            if beamline == "mx3":
                sn     = row.get("Compound", "").strip()
                target = row.get("Name",     "").strip()
                if sn and target:
                    stripped = re.sub(r"^fast_dp_results_", "", target)
                    parts = stripped.split("-")
                    crystal = "-".join(parts[:3]).lower() if len(parts) >= 3 else stripped.lower()
                    crystal_sn[crystal] = sn
            else:  # mx1
                sn     = row.get("source_directory", "").strip()
                target = row.get("target_directory", "").strip()
                if sn and target:
                    parts = target.split("_")
                    crystal = "_".join(parts[:3]) if len(parts) >= 3 else target
                    crystal_sn[crystal] = sn
    return crystal_sn


# ---------------------------------------------------------------------------
# SQLite helpers
# ---------------------------------------------------------------------------

_REQUIRED_COLUMNS = [
    ("CrystalName", "TEXT"), ("ProteinName", "TEXT"),
    ("DataCollectionVisit", "TEXT"), ("DataCollectionRun", "TEXT"),
    ("DataCollectionSubdir", "TEXT"), ("DataCollectionBeamline", "TEXT"),
    ("DataCollectionOutcome", "TEXT"), ("DataCollectionDate", "TEXT"),
    ("DataCollectionWavelength", "TEXT"),
    ("DataProcessingProgram", "TEXT"), ("DataProcessingSpaceGroup", "TEXT"),
    ("DataProcessingUnitCell", "TEXT"), ("DataProcessingAutoAssigned", "TEXT"),
    ("DataProcessingA", "TEXT"), ("DataProcessingB", "TEXT"),
    ("DataProcessingC", "TEXT"), ("DataProcessingAlpha", "TEXT"),
    ("DataProcessingBeta", "TEXT"), ("DataProcessingGamma", "TEXT"),
    ("DataProcessingResolutionOverall", "TEXT"),
    ("DataProcessingResolutionLow", "TEXT"),
    ("DataProcessingResolutionLowInnerShell", "TEXT"),
    ("DataProcessingResolutionHigh", "TEXT"),
    ("DataProcessingResolutionHigh15sigma", "TEXT"),
    ("DataProcessingResolutionHigh20sigma", "TEXT"),
    ("DataProcessingResolutionHighOuterShell", "TEXT"),
    ("DataProcessingRmergeOverall", "TEXT"), ("DataProcessingRmergeLow", "TEXT"),
    ("DataProcessingRmergeHigh", "TEXT"),
    ("DataProcessingIsigOverall", "TEXT"), ("DataProcessingIsigLow", "TEXT"),
    ("DataProcessingIsigHigh", "TEXT"),
    ("DataProcessingCompletenessOverall", "TEXT"),
    ("DataProcessingCompletenessLow", "TEXT"),
    ("DataProcessingCompletenessHigh", "TEXT"),
    ("DataProcessingMultiplicityOverall", "TEXT"),
    ("DataProcessingMultiplicityLow", "TEXT"),
    ("DataProcessingMultiplicityHigh", "TEXT"),
    ("DataProcessingCChalfOverall", "TEXT"), ("DataProcessingCChalfLow", "TEXT"),
    ("DataProcessingCChalfHigh", "TEXT"),
    ("DataProcessingPathToLogfile", "TEXT"), ("DataProcessingPathToMTZfile", "TEXT"),
    ("DataProcessingLOGfileName", "TEXT"), ("DataProcessingMTZfileName", "TEXT"),
    ("DataProcessingDirectoryOriginal", "TEXT"),
    ("DataProcessingUniqueReflectionsOverall", "TEXT"),
    ("DataProcessingUniqueReflectionsLow", "TEXT"),
    ("DataProcessingUniqueReflectionsHigh", "TEXT"),
    ("DataProcessingLattice", "TEXT"), ("DataProcessingPointGroup", "TEXT"),
    ("DataProcessingUnitCellVolume", "TEXT"),
    ("DataProcessingAlert", "TEXT"), ("DataProcessingScore", "TEXT"),
    ("DataProcessingStatus", "TEXT"),
    ("LastUpdated", "TEXT"), ("LastUpdated_by", "TEXT"),
]


def ensure_columns(conn):
    """Add any missing columns to collectionTable (safe — never removes)."""
    cur = conn.execute("PRAGMA table_info(collectionTable)")
    existing = {row[1] for row in cur.fetchall()}
    for col, coltype in _REQUIRED_COLUMNS:
        if col not in existing:
            conn.execute(f"ALTER TABLE collectionTable ADD COLUMN {col} {coltype}")
    conn.commit()


def upsert(conn, xtal, visit, run, proc_code, data):
    data["LastUpdated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    data["LastUpdated_by"] = os.environ.get("USER", "populate_xce_db")

    cur = conn.execute(
        "SELECT ID FROM collectionTable "
        "WHERE CrystalName=? AND DataCollectionRun=? AND DataCollectionSubdir=?",
        (xtal, run, proc_code),
    )
    row = cur.fetchone()

    if row is None:
        data.update({
            "CrystalName": xtal,
            "DataCollectionVisit": visit,
            "DataCollectionRun": run,
            "DataCollectionSubdir": proc_code,
        })
        cols = ", ".join(data.keys())
        placeholders = ", ".join("?" * len(data))
        conn.execute(
            f"INSERT INTO collectionTable ({cols}) VALUES ({placeholders})",
            list(data.values()),
        )
        conn.commit()
        return "INSERT"
    else:
        set_clause = ", ".join(f"{k}=?" for k in data)
        conn.execute(
            f"UPDATE collectionTable SET {set_clause} "
            f"WHERE CrystalName=? AND DataCollectionRun=? AND DataCollectionSubdir=?",
            list(data.values()) + [xtal, run, proc_code],
        )
        conn.commit()
        return "UPDATE"


def upsert_main_table(conn, xtal, compound_code, smiles, protein_name):
    """Insert or update CompoundCode and CompoundSMILES in mainTable."""
    cur = conn.execute("SELECT ID FROM mainTable WHERE CrystalName=?", (xtal,))
    row = cur.fetchone()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    user = os.environ.get("USER", "populate_xce_db")
    if row is None:
        conn.execute(
            "INSERT INTO mainTable (CrystalName, CompoundCode, CompoundSMILES, "
            "ProteinName, LastUpdated, LastUpdated_by) VALUES (?,?,?,?,?,?)",
            (xtal, compound_code, smiles, protein_name, now, user),
        )
    else:
        conn.execute(
            "UPDATE mainTable SET CompoundCode=?, CompoundSMILES=?, "
            "LastUpdated=?, LastUpdated_by=? WHERE CrystalName=?",
            (compound_code, smiles, now, user, xtal),
        )
    conn.commit()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _prompt(msg):
    """Read a path from stdin with tab completion."""
    readline.set_completer_delims(" \t\n;")
    readline.parse_and_bind("tab: complete")
    return input(msg).strip()


def main():
    print(__doc__)

    db_path     = os.path.realpath(_prompt("XCE .sqlite file path                             : "))
    proc_dir    = os.path.realpath(_prompt("processed/<target>/ dir (contains crystal subdirs) : "))
    project_dir = os.path.realpath(_prompt("XCE Project Directory                             : "))
    target      = _prompt(                 "Target / ProteinName (e.g. Bax)                   : ")
    smiles_csv  = _prompt(                 "SMILES library CSV   [blank if .smi files present]  : ")
    dist_csv    = _prompt(                 "Distribution CSV     [blank if .cmpd files present] : ")
    beamline    = _prompt(                 "Beamline mx1/mx3 (only needed if CSVs above given)  : ") or "mx1"

    if not os.path.isfile(db_path):
        raise SystemExit(f"ERROR: .sqlite file not found: {db_path}")
    if not os.path.isdir(proc_dir):
        raise SystemExit(f"ERROR: processed dir not found: {proc_dir}")
    if not os.path.isdir(project_dir):
        raise SystemExit(f"ERROR: project dir not found: {project_dir}")
    if not target:
        raise SystemExit("ERROR: target name cannot be empty.")
    if smiles_csv and not os.path.isfile(smiles_csv):
        raise SystemExit(f"ERROR: SMILES CSV not found: {smiles_csv}")
    if dist_csv and not os.path.isfile(dist_csv):
        raise SystemExit(f"ERROR: Distribution CSV not found: {dist_csv}")

    # Build compound lookup if CSVs supplied
    crystal_to_smiles = {}   # crystal_name -> (compound_code, smiles)
    if smiles_csv and dist_csv:
        print("\nLoading compound/SMILES data...")
        smiles_map  = load_smiles_map(smiles_csv)
        crystal_sn  = load_crystal_sn_map(dist_csv, beamline.lower())
        for crystal, sn in crystal_sn.items():
            smiles = smiles_map.get(sn)
            if smiles:
                crystal_to_smiles[crystal] = (sn, smiles)
            else:
                print(f"  WARNING: no SMILES for {crystal} (SN={sn})")
        print(f"  {len(crystal_to_smiles)} crystals with SMILES, "
              f"{len(crystal_sn) - len(crystal_to_smiles)} missing")
        print()

    # Derive visit to match XCE's getVisitAndBeamline logic for non-DLS paths:
    # processedDir = .../ctd-retry/processed/ctd  → visit = 'ctd-retry'  ([-3])
    parts = proc_dir.rstrip("/").split("/")
    visit = parts[-3] if len(parts) >= 3 else "local"

    print(f"Visit (derived): {visit}")
    print(f"Processed dir  : {proc_dir}")
    print(f"Project dir    : {project_dir}")
    print(f"Database       : {db_path}")
    print()

    conn = sqlite3.connect(db_path)
    ensure_columns(conn)

    args_target = target
    n_insert = n_update = n_skip = 0

    for crystal_dir in sorted(glob.glob(os.path.join(proc_dir, "*"))):
        if not os.path.isdir(crystal_dir):
            continue
        xtal = os.path.basename(crystal_dir)

        for run_dir in sorted(glob.glob(os.path.join(crystal_dir, "*"))):
            if not os.path.isdir(run_dir):
                continue
            run = os.path.basename(run_dir)

            for code_dir in sorted(glob.glob(os.path.join(run_dir, "*"))):
                if not os.path.isdir(code_dir) or os.path.islink(code_dir):
                    continue
                proc_code = os.path.basename(code_dir)

                # Locate aimless.log
                log_matches = glob.glob(
                    os.path.join(code_dir, "output", "LogFiles", "*aimless.log"))
                if not log_matches:
                    print(f"  SKIP  {xtal}/{run}/{proc_code}  (no aimless.log)")
                    n_skip += 1
                    continue
                src_log = log_matches[0]

                # Locate *free.mtz
                mtz_matches = glob.glob(
                    os.path.join(code_dir, "output", "DataFiles", "*free.mtz"))
                if not mtz_matches:
                    print(f"  SKIP  {xtal}/{run}/{proc_code}  (no *free.mtz)")
                    n_skip += 1
                    continue
                src_mtz = mtz_matches[0]

                # Determine program from proc_code
                if "fast_dp" in proc_code:
                    program = "fast_dp"
                elif "mx1" in proc_code:
                    program = "mx1_autoproc"
                else:
                    program = proc_code

                # Replicate XCE's autoprocessing subdir naming:
                # <visit>-<run><DataProcessingProgram>_<DataCollectionSubdir>
                # XCE's checkExistingFiles() reconstructs this path from DB fields,
                # so autoproc_tag MUST equal DataProcessingProgram (not "unknown").
                ap_subdir = f"{visit}-{run}{program}_{proc_code}"
                dest_dir = os.path.join(project_dir, xtal, "autoprocessing", ap_subdir)
                os.makedirs(dest_dir, exist_ok=True)

                # Copy log and mtz into project autoprocessing dir
                dest_log = os.path.join(dest_dir, os.path.basename(src_log))
                dest_mtz = os.path.join(dest_dir, os.path.basename(src_mtz))
                if not os.path.exists(dest_log):
                    shutil.copy2(src_log, dest_log)
                if not os.path.exists(dest_mtz):
                    shutil.copy2(src_mtz, dest_mtz)

                # Convenience symlinks <xtal>.log / <xtal>.mtz (XCE creates these too)
                for link_name, target_name in [
                    (xtal + ".log", os.path.basename(dest_log)),
                    (xtal + ".mtz", os.path.basename(dest_mtz)),
                ]:
                    link_path = os.path.join(dest_dir, link_name)
                    if not os.path.exists(link_path):
                        os.symlink(target_name, link_path)

                # Parse stats
                stats = parse_aimless_log(dest_log)
                score, vol = compute_score_and_volume(stats)
                alert = compute_alert(stats)

                a, b, c = stats["DataProcessingA"], stats["DataProcessingB"], stats["DataProcessingC"]
                al, be, ga = stats["DataProcessingAlpha"], stats["DataProcessingBeta"], stats["DataProcessingGamma"]
                unit_cell = f"{a} {b} {c} {al} {be} {ga}"
                res_overall = f"{stats['DataProcessingResolutionLow']} - {stats['DataProcessingResolutionHigh']}"
                timestamp = datetime.fromtimestamp(os.path.getmtime(run_dir)).strftime(
                    "%Y-%m-%d %H:%M:%S")

                data = {
                    "ProteinName":                         args_target,
                    "DataCollectionBeamline":              "unknown",
                    "DataCollectionOutcome":               "success",
                    "DataCollectionDate":                  timestamp,
                    "DataCollectionWavelength":            stats["DataCollectionWavelength"],
                    "DataProcessingProgram":               program,
                    "DataProcessingSpaceGroup":            stats["DataProcessingSpaceGroup"],
                    "DataProcessingUnitCell":              unit_cell,
                    "DataProcessingA":                     a,
                    "DataProcessingB":                     b,
                    "DataProcessingC":                     c,
                    "DataProcessingAlpha":                 al,
                    "DataProcessingBeta":                  be,
                    "DataProcessingGamma":                 ga,
                    "DataProcessingResolutionOverall":     res_overall,
                    "DataProcessingResolutionLow":         stats["DataProcessingResolutionLow"],
                    "DataProcessingResolutionLowInnerShell": stats["DataProcessingResolutionLowInnerShell"],
                    "DataProcessingResolutionHigh":        stats["DataProcessingResolutionHigh"],
                    "DataProcessingResolutionHigh15sigma": stats["DataProcessingResolutionHigh15sigma"],
                    "DataProcessingResolutionHigh20sigma": stats["DataProcessingResolutionHigh20sigma"],
                    "DataProcessingResolutionHighOuterShell": stats["DataProcessingResolutionHighOuterShell"],
                    "DataProcessingRmergeOverall":         stats["DataProcessingRmergeOverall"],
                    "DataProcessingRmergeLow":             stats["DataProcessingRmergeLow"],
                    "DataProcessingRmergeHigh":            stats["DataProcessingRmergeHigh"],
                    "DataProcessingIsigOverall":           stats["DataProcessingIsigOverall"],
                    "DataProcessingIsigLow":               stats["DataProcessingIsigLow"],
                    "DataProcessingIsigHigh":              stats["DataProcessingIsigHigh"],
                    "DataProcessingCompletenessOverall":   stats["DataProcessingCompletenessOverall"],
                    "DataProcessingCompletenessLow":       stats["DataProcessingCompletenessLow"],
                    "DataProcessingCompletenessHigh":      stats["DataProcessingCompletenessHigh"],
                    "DataProcessingMultiplicityOverall":   stats["DataProcessingMultiplicityOverall"],
                    "DataProcessingMultiplicityLow":       stats["DataProcessingMultiplicityLow"],
                    "DataProcessingMultiplicityHigh":      stats["DataProcessingMultiplicityHigh"],
                    "DataProcessingCChalfOverall":         stats["DataProcessingCChalfOverall"],
                    "DataProcessingCChalfLow":             stats["DataProcessingCChalfLow"],
                    "DataProcessingCChalfHigh":            stats["DataProcessingCChalfHigh"],
                    "DataProcessingPathToLogfile":         dest_log,
                    "DataProcessingPathToMTZfile":         dest_mtz,
                    "DataProcessingLOGfileName":           os.path.basename(dest_log),
                    "DataProcessingMTZfileName":           os.path.basename(dest_mtz),
                    "DataProcessingDirectoryOriginal":     os.path.join(code_dir, "output"),
                    "DataProcessingUniqueReflectionsOverall": stats["DataProcessingUniqueReflectionsOverall"],
                    "DataProcessingUniqueReflectionsLow":  stats["DataProcessingUniqueReflectionsLow"],
                    "DataProcessingUniqueReflectionsHigh": stats["DataProcessingUniqueReflectionsHigh"],
                    "DataProcessingLattice":               stats["DataProcessingLattice"],
                    "DataProcessingPointGroup":            stats["DataProcessingPointGroup"],
                    "DataProcessingUnitCellVolume":        vol,
                    "DataProcessingAlert":                 alert,
                    "DataProcessingScore":                 score,
                    "DataProcessingStatus":                "running",
                }

                action = upsert(conn, xtal, visit, run, proc_code, data)
                res_str = stats["DataProcessingResolutionHigh"]
                sg_str  = stats["DataProcessingSpaceGroup"]

                # Populate mainTable with compound info.
                # Priority 1: .smi / .cmpd files written by prepare_*_for_xce.sh
                # Priority 2: CSV lookup (if CSVs were supplied at prompt)
                compound_note = ""
                smi_file  = os.path.join(proc_dir, xtal, xtal + ".smi")
                cmpd_file = os.path.join(proc_dir, xtal, xtal + ".cmpd")
                if os.path.isfile(smi_file) and os.path.isfile(cmpd_file):
                    with open(smi_file)  as f: file_smiles = f.read().strip()
                    with open(cmpd_file) as f: file_cmpd   = f.read().strip()
                    if file_smiles and file_cmpd:
                        upsert_main_table(conn, xtal, file_cmpd, file_smiles, args_target)
                        compound_note = f"  cmpd={file_cmpd} (from file)"
                elif crystal_to_smiles:
                    # Fall back to CSV lookup
                    entry = crystal_to_smiles.get(xtal) or crystal_to_smiles.get(xtal.lower())
                    if entry:
                        compound_code, smiles = entry
                        upsert_main_table(conn, xtal, compound_code, smiles, args_target)
                        compound_note = f"  cmpd={compound_code} (from CSV)"
                    else:
                        compound_note = "  cmpd=NOT_FOUND"

                print(f"  {action:6s}  {xtal}  run={run}  {proc_code}  "
                      f"res={res_str}  SG={sg_str}{compound_note}")
                if action == "INSERT":
                    n_insert += 1
                else:
                    n_update += 1

    conn.close()
    print()
    print(f"Done: {n_insert} inserted, {n_update} updated, {n_skip} skipped.")
    print()
    print("Next steps in XCE:")
    print("  1. Open XCE -> Datasets tab -> select target from dropdown.")
    print("  2. Click 'Select Best Autoprocessing Result' (Run button, NOT Status).")
    print("  3. Do NOT click Status buttons — they crash at WEHI (CLUSTER_BASTION undefined).")


if __name__ == "__main__":
    main()
