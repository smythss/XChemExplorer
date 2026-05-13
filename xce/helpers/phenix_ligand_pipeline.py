#!/usr/bin/env python3
"""
phenix_ligand_pipeline.py

Per-dataset pipeline for fragment/ligand soaking experiments:

  1. Run phenix.ligand_pipeline with the reference PDB, the dataset MTZ, a FASTA
     sequence file, and *.cif ligand restraints found in the dataset directory.
     R-free flags are sourced from the reference MTZ via
     xray_data.r_free_flags.file_name.  Anisotropic ADP refinement
     (refine.after_mr.adp_type=aniso) is applied automatically when the dataset
     resolution is below the configured cutoff (default 1.5 Å).
     phenix.ligand_pipeline auto-creates a pipeline_N/ subdirectory in cwd.

  2. Post-process outputs into an XCE-compatible Refine_N/ directory structure:

       <dataset_dir>/Refine_N/refine_N.pdb  →  ../pipeline_N/<cmpd>_final.pdb
       <dataset_dir>/Refine_N/refine_N.mtz  →  ../pipeline_N/<cmpd>_final.mtz
       <dataset_dir>/refine.pdb             →  ./Refine_N/refine_N.pdb
       <dataset_dir>/refine.mtz             →  ./Refine_N/refine_N.mtz
       <dataset_dir>/<xtal>.free.mtz        →  <xtal>_rfree_transferred.mtz

     XCE discovers refinement cycles by scanning for Refine_N/ directories and
     reads R-values / resolution from REMARK 3 lines in refine.pdb.

Intended to be submitted as one SLURM job per dataset by
run_ligand_pipeline_batch.sbatch.

Usage:
    python phenix_ligand_pipeline.py \\
        --dataset_dir /path/to/dataset \\
        --ref_mtz     /path/to/reference.mtz \\
        --ref_pdb     /path/to/reference.pdb \\
        --seq_file    /path/to/seq.fa \\
        [--mtz_pattern init.mtz] \\
        [--rfree_label FreeR_flag] \\
        [--rfree_fraction 0.05] \\
        [--nproc 4] \\
        [--phenix_bin /opt/phenix/bin] \\
        [--ligand_pipeline_output_pdb path/to/override.pdb] \\
        [--ligand_pipeline_output_mtz path/to/override.mtz] \\
        [--skip_rfree_transfer] \\
        [--skip_xce_output] \\
        [--dry_run] \\
        [--verbose]
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# ---------- defaults ----------
DEFAULT_MTZ_PATTERN = "init.mtz"
DEFAULT_RFREE_LABEL = "FreeR_flag"
DEFAULT_RFREE_FRACTION = 0.1
DEFAULT_NPROC = 4
DEFAULT_ANISO_CUTOFF = 1.5  # Å — use aniso ADP refinement below this resolution


# ---------- helpers ----------

def find_executable(name: str, bin_dir: Optional[str] = None) -> Optional[str]:
    """Return path to executable, checking bin_dir first then PATH."""
    if bin_dir:
        candidate = Path(bin_dir) / name
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return str(candidate)
    return shutil.which(name)


def run_command(
    cmd: List[str],
    cwd: Optional[str] = None,
    logfile: Optional[str] = None,
    dry_run: bool = False,
) -> Tuple[int, str, str]:
    """Run a subprocess command. Returns (returncode, stdout, stderr)."""
    cmd_str = " ".join(str(c) for c in cmd)
    logging.info("CMD: %s", cmd_str)
    if dry_run:
        return 0, "", ""

    result = subprocess.run(
        [str(c) for c in cmd],
        capture_output=True,
        text=True,
        cwd=cwd,
    )
    stdout = result.stdout or ""
    stderr = result.stderr or ""

    if logfile:
        with open(logfile, "w") as fh:
            fh.write(f"CMD: {cmd_str}\n")
            fh.write(f"\nSTDOUT:\n{stdout}\n")
            fh.write(f"\nSTDERR:\n{stderr}\n")
            fh.write(f"\nReturn code: {result.returncode}\n")

    if stdout:
        logging.debug("STDOUT: %s", stdout[:4000])
    if stderr:
        logging.debug("STDERR: %s", stderr[:2000])

    return result.returncode, stdout, stderr


def get_next_refine_serial(dataset_dir: Path) -> int:
    """Return the next Refine_N serial, mirroring XCE's GetSerial().

    Scans for Refine_<N> directories (ignoring -report suffixes) and returns
    max(N) + 1, or 1 if none exist.
    """
    serials = []
    for item in dataset_dir.iterdir():
        if not item.name.startswith("Refine_"):
            continue
        if item.name.endswith("-report"):
            continue
        try:
            serials.append(int(item.name[item.name.rfind("_") + 1:]))
        except ValueError:
            continue
    return max(serials) + 1 if serials else 1


def find_dataset_mtz(dataset_dir: Path, pattern: str) -> Optional[Path]:
    """Find the dataset MTZ using a glob pattern.

    Returns the first match (symlinks are not resolved; Phenix handles them).
    """
    matches = sorted(dataset_dir.glob(pattern))
    if not matches:
        return None
    if len(matches) > 1:
        logging.warning(
            "Multiple MTZ files matching %r in %s; using: %s",
            pattern, dataset_dir, matches[0].name,
        )
    return matches[0]


def find_cif_files(dataset_dir: Path) -> List[Path]:
    """Return sorted *.cif files in dataset_dir (non-recursive)."""
    return sorted(p for p in dataset_dir.glob("*.cif") if p.exists())


def get_mtz_dmin(mtz_path: Path) -> Optional[float]:
    """Return the high-resolution limit (d_min, Å) of an MTZ file.

    Uses iotbx.mtz directly, which is available inside the Phenix Python
    environment.  Returns None if the value cannot be determined (e.g. outside
    the Phenix environment or if the file is unreadable).
    """
    try:
        from iotbx import mtz as iotbx_mtz  # noqa: PLC0415
        mtz_obj = iotbx_mtz.object(file_name=str(mtz_path))
        arrays = mtz_obj.as_miller_arrays()
        if arrays:
            return float(min(a.d_min() for a in arrays))
    except Exception as exc:
        logging.debug("Could not determine d_min from %s: %s", mtz_path.name, exc)
    return None


# ---------- step 1: R-free flag transfer ----------

def transfer_rfree_flags(
    ref_mtz: Path,
    new_mtz: Path,
    output_mtz: Path,
    rfree_label: str,
    rfree_fraction: float,
    phenix_bin: Optional[str],
    work_dir: Path,
    dry_run: bool,
) -> bool:
    """Ensure the dataset MTZ has R-free flags, writing the result to output_mtz.

    Strategy:
      - If the input MTZ already contains the R-free column (rfree_label), copy
        it to output_mtz unchanged.  Writing a new column with the same name into
        an MTZ that already has it causes a duplicate-column assertion error in
        phenix.reflection_file_converter.
      - If no R-free column is found, generate a fresh set using
        phenix.reflection_file_converter --generate_r_free_flags.

    Returns True on success.
    """
    exe = find_executable("phenix.reflection_file_converter", phenix_bin)
    if not exe:
        logging.error(
            "phenix.reflection_file_converter not found. "
            "Check PATH or pass --phenix_bin."
        )
        return False

    log_path = str(work_dir / "step1_rfree_transfer.log")

    # Detect whether the input MTZ already has the R-free column by listing its
    # arrays (phenix.reflection_file_converter with no output flags prints them).
    has_rfree = False
    if not dry_run:
        probe = subprocess.run(
            [str(exe), str(new_mtz)], capture_output=True, text=True
        )
        # Output lines contain e.g. "Miller array info: /path/file.mtz:FreeR_flag"
        has_rfree = f":{rfree_label}" in probe.stdout
        logging.debug(
            "R-free detection (%s present: %s):\n%s",
            rfree_label, has_rfree, probe.stdout[:500],
        )

    if has_rfree or dry_run:
        # Existing R-free flags — copy MTZ as-is, preserving all columns
        logging.info(
            "Input MTZ already contains %s — copying without modification.",
            rfree_label,
        )
        logging.info("CMD: cp %s %s", new_mtz, output_mtz)
        if not dry_run:
            shutil.copy2(str(new_mtz), str(output_mtz))
    else:
        # No R-free flags found — generate a fresh set
        cmd = [
            exe,
            str(new_mtz),
            "--generate_r_free_flags",
            f"--r_free_flags_fraction={rfree_fraction}",
            f"--output_r_free_label={rfree_label}",
            f"--mtz={output_mtz}",
        ]
        rc, _, _ = run_command(cmd, cwd=str(work_dir), logfile=log_path, dry_run=dry_run)
        if rc != 0:
            logging.error(
                "phenix.reflection_file_converter failed (rc=%d). See %s", rc, log_path
            )
            return False

    logging.info("R-free flags ready → %s", output_mtz.name)
    return True


# ---------- step 2: phenix.ligand_pipeline ----------

def run_ligand_pipeline(
    ref_pdb: Path,
    mtz: Path,
    seq_file: Optional[Path],
    cif_files: List[Path],
    dataset_dir: Path,
    nproc: int,
    phenix_bin: Optional[str],
    dry_run: bool,
    ref_mtz: Optional[Path] = None,
    rfree_label: str = DEFAULT_RFREE_LABEL,
    xray_data_labels: Optional[str] = None,
    aniso_cutoff: float = DEFAULT_ANISO_CUTOFF,
) -> Tuple[bool, Optional[Path]]:
    """Run phenix.ligand_pipeline from dataset_dir as cwd.

    R-free flags are transferred internally by phenix.ligand_pipeline via
    xray_data.r_free_flags.file_name when ref_mtz is provided — no separate
    conversion step is needed.

    Anisotropic ADP refinement (refine.after_mr.adp_type=aniso) is enabled
    automatically when the dataset d_min is below aniso_cutoff (default 1.5 Å).

    phenix.ligand_pipeline auto-creates a pipeline_N/ subdirectory in cwd.
    We snapshot existing pipeline_N/ dirs before the run to detect which new
    one appeared.

    Returns (success, pipeline_output_dir).
    """
    exe = find_executable("phenix.ligand_pipeline", phenix_bin)
    if not exe:
        logging.error(
            "phenix.ligand_pipeline not found. "
            "Check PATH or pass --phenix_bin."
        )
        return False, None

    # Snapshot existing pipeline_N/ dirs before run
    existing_pipeline_dirs = {
        p for p in dataset_dir.glob("pipeline_*/") if p.is_dir()
    }

    # Build command — positional file args followed by Phil keyword arguments.
    cmd = [exe, str(ref_pdb), str(mtz)]
    if seq_file is not None:
        cmd.append(str(seq_file))
    for cif in cif_files:
        cmd.append(str(cif))
    cmd.append(f"runtime.nproc={nproc}")

    # Data column labels (e.g. IMEAN,SIGIMEAN for synchrotron intensity data)
    if xray_data_labels:
        cmd.append(f"xray_data.labels={xray_data_labels}")

    # R-free flags sourced from reference MTZ (no separate transfer step needed)
    if ref_mtz is not None:
        cmd.append(f"xray_data.r_free_flags.file_name={ref_mtz}")
        cmd.append(f"xray_data.r_free_flags.label={rfree_label}")

    # Anisotropic ADP refinement for high-resolution datasets
    dmin = get_mtz_dmin(mtz)
    if dmin is not None:
        logging.info("Dataset resolution: %.3f \u00c5", dmin)
        if dmin < aniso_cutoff:
            logging.info(
                "Resolution %.3f \u00c5 < %.1f \u00c5 cutoff \u2014 "
                "adding refine.after_mr.adp_type=aniso",
                dmin, aniso_cutoff,
            )
            cmd.append("refine.after_mr.adp_type=aniso")
        else:
            logging.info(
                "Resolution %.3f \u00c5 >= %.1f \u00c5 \u2014 using default isotropic ADP refinement",
                dmin, aniso_cutoff,
            )
    else:
        logging.warning(
            "Could not determine dataset resolution from %s; "
            "using default ADP refinement.",
            mtz.name,
        )

    log_path = str(dataset_dir / "step2_ligand_pipeline.log")
    rc, _, _ = run_command(
        cmd, cwd=str(dataset_dir), logfile=log_path, dry_run=dry_run
    )

    if dry_run:
        return True, None

    # Detect newly created pipeline_N/ directory
    all_pipeline_dirs = {
        p for p in dataset_dir.glob("pipeline_*/") if p.is_dir()
    }
    new_dirs = sorted(all_pipeline_dirs - existing_pipeline_dirs)
    if not new_dirs:
        logging.error(
            "Could not find a new pipeline_N/ directory after the run. "
            "Check %s for details.",
            log_path,
        )
        return False, None

    # In the unlikely event multiple new dirs appeared, take the highest-numbered
    pipeline_dir = new_dirs[-1]
    logging.info("phenix.ligand_pipeline output directory: %s", pipeline_dir.name)

    # Check for FAILED flag (ligand fitting failed, but a refined model exists)
    if (pipeline_dir / "FAILED").exists():
        logging.warning(
            "phenix.ligand_pipeline set FAILED flag in %s — "
            "ligand fitting was unsuccessful. "
            "Continuing with the refined model; inspect maps carefully.",
            pipeline_dir.name,
        )

    if rc != 0:
        logging.error(
            "phenix.ligand_pipeline exited with rc=%d. "
            "See %s and %s/pipeline.log",
            rc, log_path, pipeline_dir.name,
        )
        # Still return the directory so the caller can report where outputs are
        return False, pipeline_dir

    return True, pipeline_dir


# ---------- step 3: locate best output PDB/MTZ from pipeline_N ----------

def find_pipeline_outputs(
    pipeline_dir: Path,
    override_pdb: Optional[Path] = None,
    override_mtz: Optional[Path] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Locate the final PDB and MTZ produced by phenix.ligand_pipeline.

    Resolution order:
      1. --ligand_pipeline_output_pdb / --ligand_pipeline_output_mtz overrides
      2. *_final.pdb / *_final.mtz  (standard phenix.ligand_pipeline naming,
         e.g. SN02730832_final.pdb)
      3. Most recently modified *.pdb / *.mtz (fallback)

    Returns (pdb_path, mtz_path); either may be None if not found.
    """
    def _find(
        directory: Path,
        pattern: str,
        ext: str,
        override: Optional[Path],
    ) -> Optional[Path]:
        if override is not None:
            if override.exists():
                return override
            logging.error("Override file not found: %s", override)
            return None
        # Standard naming first
        matches = sorted(directory.glob(pattern))
        if matches:
            if len(matches) > 1:
                logging.warning(
                    "Multiple %s files in %s; using: %s",
                    pattern, directory.name, matches[0].name,
                )
            return matches[0]
        # Fallback: most recently modified file with this extension
        all_files = sorted(
            (p for p in directory.glob(f"*{ext}") if p.is_file()),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if all_files:
            logging.warning(
                "No %s found in %s; using most-recently-modified fallback: %s. "
                "Use --ligand_pipeline_output_%s to specify explicitly.",
                pattern, directory.name, all_files[0].name,
                ext.lstrip("."),
            )
            return all_files[0]
        return None

    pdb = _find(pipeline_dir, "*_final.pdb", ".pdb", override_pdb)
    mtz = _find(pipeline_dir, "*_final.mtz", ".mtz", override_mtz)
    return pdb, mtz


# ---------- step 3b: best-model selection across pipeline stages ----------

_RWORK_RE = re.compile(
    r"R VALUE\s+\(WORKING \+ TEST SET\)\s*:\s*([0-9]+\.[0-9]+)"
)
_RFREE_RE = re.compile(
    r"FREE R VALUE\s*:\s*([0-9]+\.[0-9]+)"
)
# Lines that mention these strings but are NOT the scalar R-free line
_RFREE_SKIP = ("TEST SET SIZE", "TEST SET COUNT")


def parse_rvalues_from_pdb(pdb_path: Path) -> Tuple[Optional[float], Optional[float]]:
    """Parse Rwork and Rfree from REMARK 3 lines in a PDB file.

    Returns (rwork, rfree); either value is None if not found.
    Only reads up to the first ATOM/HETATM record for speed.
    """
    rwork = rfree = None
    try:
        with open(pdb_path) as fh:
            for raw in fh:
                if raw.startswith(("ATOM", "HETATM")):
                    break
                if not raw.startswith("REMARK"):
                    continue
                if rwork is None:
                    m = _RWORK_RE.search(raw)
                    if m:
                        rwork = float(m.group(1))
                if rfree is None:
                    if not any(skip in raw for skip in _RFREE_SKIP):
                        m = _RFREE_RE.search(raw)
                        if m:
                            rfree = float(m.group(1))
                if rwork is not None and rfree is not None:
                    break
    except OSError:
        pass
    return rwork, rfree


def find_companion_mtz(pdb_path: Path) -> Optional[Path]:
    """Find the MTZ file that accompanies a given PDB in the same directory.

    Tries (in order):
      1. Exact stem match:  <stem>.mtz
      2. Longest common prefix among all *.mtz files in the directory
    """
    directory = pdb_path.parent
    stem = pdb_path.stem

    exact = directory / f"{stem}.mtz"
    if exact.exists():
        return exact

    candidates = list(directory.glob("*.mtz"))
    if not candidates:
        return None

    # Pick the MTZ whose stem shares the most leading characters with the PDB stem
    candidates.sort(
        key=lambda m: len(os.path.commonprefix([stem, m.stem])),
        reverse=True,
    )
    return candidates[0]


def select_best_pipeline_output(
    pipeline_dir: Path,
) -> Tuple[Optional[Path], Optional[Path], Optional[float], Optional[float]]:
    """Scan all PDB files under pipeline_dir and return the pair with lowest Rfree.

    Walks the directory tree recursively.  For each PDB that contains REMARK 3
    R-values, looks for a companion MTZ in the same directory.  Returns the
    (pdb, mtz, rwork, rfree) tuple with the lowest Rfree value, or (None, None,
    None, None) if no valid pair is found.

    This is particularly useful when phenix.ligand_pipeline fails to place the
    ligand: the intermediate refine_0/ model (protein + waters only) typically
    has far better R-values than the forced-placement final model.
    """
    best_pdb: Optional[Path] = None
    best_mtz: Optional[Path] = None
    best_rwork: Optional[float] = None
    best_rfree: Optional[float] = None

    candidates: List[Tuple[float, float, Path, Path]] = []

    for pdb in sorted(pipeline_dir.rglob("*.pdb")):
        if not pdb.is_file():
            continue
        rwork, rfree = parse_rvalues_from_pdb(pdb)
        if rfree is None:
            # No REMARK 3 refinement statistics — skip (e.g. ligand-only PDBs,
            # prepare_models outputs, LigandFit intermediate files)
            logging.debug("No R-values in %s — skipping", pdb.relative_to(pipeline_dir))
            continue
        mtz = find_companion_mtz(pdb)
        if mtz is None:
            logging.debug(
                "No companion MTZ for %s — skipping",
                pdb.relative_to(pipeline_dir),
            )
            continue
        logging.debug(
            "Candidate: %s  Rwork=%.4f Rfree=%.4f  MTZ=%s",
            pdb.relative_to(pipeline_dir), rwork or 0.0, rfree,
            mtz.relative_to(pipeline_dir),
        )
        candidates.append((rfree, rwork or 1.0, pdb, mtz))

    if not candidates:
        return None, None, None, None

    # Sort by Rfree ascending; break ties with Rwork
    candidates.sort(key=lambda t: (t[0], t[1]))
    best_rfree, best_rwork, best_pdb, best_mtz = candidates[0]

    logging.info(
        "Best model selected: %s  Rwork=%.4f  Rfree=%.4f",
        best_pdb.relative_to(pipeline_dir), best_rwork, best_rfree,
    )
    if len(candidates) > 1:
        logging.info("All candidates (ascending Rfree):")
        for rfree, rwork, pdb, mtz in candidates:
            marker = " ← selected" if pdb == best_pdb else ""
            logging.info(
                "  Rfree=%.4f  Rwork=%.4f  %s%s",
                rfree, rwork, pdb.relative_to(pipeline_dir), marker,
            )

    return best_pdb, best_mtz, best_rwork, best_rfree


# ---------- step 4: XCE-compatible post-processing ----------

def make_symlink(target: Path, link_path: Path, dry_run: bool) -> None:
    """Create or replace a symlink at link_path pointing to target.

    Uses a relative target so symlinks remain valid if the project directory
    is moved.
    """
    try:
        rel_target = os.path.relpath(str(target), str(link_path.parent))
    except ValueError:
        rel_target = str(target)

    if dry_run:
        logging.info("SYMLINK (dry-run): %s → %s", link_path.name, rel_target)
        return

    if link_path.is_symlink() or link_path.exists():
        link_path.unlink()
    link_path.symlink_to(rel_target)
    logging.info("Symlink: %s → %s", link_path.name, rel_target)


def create_xce_output(
    dataset_dir: Path,
    xtal_name: str,
    pipeline_pdb: Path,
    pipeline_mtz: Path,
    rfree_mtz: Path,
    serial: int,
    dry_run: bool,
) -> bool:
    """Build the Refine_N/ structure and root-level symlinks expected by XCE.

    Creates:
      Refine_N/refine_N.pdb   →  ../pipeline_N/<compound>_final.pdb
      Refine_N/refine_N.mtz   →  ../pipeline_N/<compound>_final.mtz
      refine.pdb              →  ./Refine_N/refine_N.pdb
      refine.mtz              →  ./Refine_N/refine_N.mtz
      <xtal>.free.mtz         →  <xtal>_rfree_transferred.mtz
        (updated to reflect the transferred flags so XCE refinement cycles
         use the same consistent R-free set as the reference)
    """
    refine_dir = dataset_dir / f"Refine_{serial}"
    refine_pdb = refine_dir / f"refine_{serial}.pdb"
    refine_mtz = refine_dir / f"refine_{serial}.mtz"

    if not dry_run:
        refine_dir.mkdir(parents=True, exist_ok=True)

    # Refine_N/refine_N.pdb → pipeline output PDB
    make_symlink(pipeline_pdb, refine_pdb, dry_run)
    # Refine_N/refine_N.mtz → pipeline output MTZ
    make_symlink(pipeline_mtz, refine_mtz, dry_run)

    # dataset_root/refine.pdb → ./Refine_N/refine_N.pdb  (XCE reads this)
    make_symlink(refine_pdb, dataset_dir / "refine.pdb", dry_run)
    # dataset_root/refine.mtz → ./Refine_N/refine_N.mtz  (XCE reads this)
    make_symlink(refine_mtz, dataset_dir / "refine.mtz", dry_run)

    # dataset_root/<xtal>.free.mtz → rfree-transferred MTZ
    # Updating this ensures subsequent XCE refinement cycles use the consistent
    # R-free set from the reference.  Log a clear warning so the user knows.
    free_mtz_link = dataset_dir / f"{xtal_name}.free.mtz"
    if (free_mtz_link.is_symlink() or free_mtz_link.exists()) and not dry_run:
        old_target = (
            os.readlink(str(free_mtz_link))
            if free_mtz_link.is_symlink()
            else str(free_mtz_link)
        )
        logging.warning(
            "%s.free.mtz already exists (→ %s); "
            "replacing with R-free-transferred MTZ.",
            xtal_name, old_target,
        )
    make_symlink(rfree_mtz, free_mtz_link, dry_run)

    return True


# ---------- argument parsing ----------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group("required arguments")
    req.add_argument(
        "--dataset_dir", required=True,
        help="Dataset directory containing the MTZ file and *.cif restraints.",
    )
    req.add_argument(
        "--ref_mtz", default=None,
        help=(
            "Reference MTZ file for R-free flag transfer (optional). "
            "If omitted and the dataset MTZ already contains R-free flags "
            "they are preserved as-is."
        ),
    )
    req.add_argument(
        "--ref_pdb", required=True,
        help="Reference PDB used as starting model for phenix.ligand_pipeline.",
    )
    req.add_argument(
        "--seq_file", default=None,
        help=(
            "FASTA sequence file (seq.fa) for phenix.ligand_pipeline (optional). "
            "If omitted, phenix.ligand_pipeline runs without a sequence file."
        ),
    )

    opt = p.add_argument_group("optional arguments")
    opt.add_argument(
        "--mtz_pattern", default=DEFAULT_MTZ_PATTERN,
        help=(
            f"Glob pattern to find the dataset MTZ within --dataset_dir "
            f"(default: {DEFAULT_MTZ_PATTERN}). "
            f"'init.mtz' is the standard XCE symlink to the dimple output MTZ."
        ),
    )
    opt.add_argument(
        "--rfree_label", default=DEFAULT_RFREE_LABEL,
        help=f"R-free column label in reference MTZ (default: {DEFAULT_RFREE_LABEL}).",
    )
    opt.add_argument(
        "--rfree_fraction", type=float, default=DEFAULT_RFREE_FRACTION,
        help=(
            f"Fraction used when assigning R-free flags to reflections absent "
            f"from the reference (default: {DEFAULT_RFREE_FRACTION})."
        ),
    )
    opt.add_argument(
        "--nproc", type=int, default=DEFAULT_NPROC,
        help=f"Processors for phenix.ligand_pipeline (default: {DEFAULT_NPROC}).",
    )
    opt.add_argument(
        "--phenix_bin", default=None,
        help=(
            "Directory containing Phenix executables (optional). "
            "Used for both phenix.reflection_file_converter and "
            "phenix.ligand_pipeline."
        ),
    )
    opt.add_argument(
        "--ligand_pipeline_output_pdb", default=None,
        metavar="PATH",
        help=(
            "Override auto-detection: explicit path to the PDB produced by "
            "phenix.ligand_pipeline to use as refine_N.pdb. "
            "Use this if the standard *_final.pdb pattern does not match."
        ),
    )
    opt.add_argument(
        "--ligand_pipeline_output_mtz", default=None,
        metavar="PATH",
        help=(
            "Override auto-detection: explicit path to the MTZ produced by "
            "phenix.ligand_pipeline to use as refine_N.mtz. "
            "Use this if the standard *_final.mtz pattern does not match."
        ),
    )
    opt.add_argument(
        "--select_best_model", action="store_true",
        help=(
            "Scan all PDB/MTZ pairs produced inside the pipeline_N/ output "
            "directory and select the one with the lowest R-free, rather than "
            "always using *_final.pdb. Useful when ligand fitting fails: the "
            "intermediate refine_0/ model (protein + waters, no forced ligand) "
            "typically has far better statistics and is more suitable for "
            "PanDDA / XCE analysis. When used together with "
            "--ligand_pipeline_output_pdb/mtz those overrides take priority."
        ),
    )
    opt.add_argument(
        "--xray_data_labels", default=None,
        metavar="LABELS",
        help=(
            "Reflection data column labels passed to phenix.ligand_pipeline as "
            "xray_data.labels (e.g. 'IMEAN,SIGIMEAN' for synchrotron intensity "
            "data, 'F,SIGF' for amplitudes). Omit to let phenix.ligand_pipeline "
            "auto-detect."
        ),
    )
    opt.add_argument(
        "--aniso_cutoff", type=float, default=DEFAULT_ANISO_CUTOFF,
        metavar="ANGSTROM",
        help=(
            f"Resolution cutoff in \u00c5 below which anisotropic ADP refinement "
            f"(refine.after_mr.adp_type=aniso) is enabled automatically. "
            f"Default: {DEFAULT_ANISO_CUTOFF} \u00c5."
        ),
    )
    opt.add_argument(
        "--skip_rfree_transfer", action="store_true",
        help=(
            "Deprecated / no-op. R-free flags are now handled internally by "
            "phenix.ligand_pipeline via xray_data.r_free_flags.file_name. "
            "Accepted for backward compatibility with existing batch scripts."
        ),
    )
    opt.add_argument(
        "--skip_xce_output", action="store_true",
        help=(
            "Skip creation of Refine_N/ structure and root-level symlinks. "
            "For non-XCE workflows."
        ),
    )
    opt.add_argument(
        "--dry_run", action="store_true",
        help="Print commands without executing them.",
    )
    opt.add_argument(
        "--verbose", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    return p.parse_args()


# ---------- main ----------

def main() -> int:
    args = parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    dataset_dir = Path(args.dataset_dir).resolve()
    ref_mtz     = Path(args.ref_mtz).resolve() if args.ref_mtz else None
    ref_pdb     = Path(args.ref_pdb).resolve()
    seq_file    = Path(args.seq_file).resolve() if args.seq_file else None
    xtal_name   = dataset_dir.name

    logging.info("=" * 60)
    logging.info("Dataset       : %s", dataset_dir)
    logging.info("Reference MTZ : %s", ref_mtz or "(not provided)")
    logging.info("Reference PDB : %s", ref_pdb)
    logging.info("Sequence file : %s", seq_file or "(not provided)")
    logging.info("=" * 60)

    # --- Validate required inputs ---
    errors = 0
    for path, label in [
        (dataset_dir, "--dataset_dir"),
        (ref_pdb,     "--ref_pdb"),
    ]:
        if not path.exists():
            logging.error("Not found: %s (%s)", path, label)
            errors += 1
    for path, label in [
        (ref_mtz,  "--ref_mtz"),
        (seq_file, "--seq_file"),
    ]:
        if path is not None and not path.exists():
            logging.error("Not found: %s (%s)", path, label)
            errors += 1
    if errors:
        return 1

    # --- Locate dataset MTZ ---
    dataset_mtz = find_dataset_mtz(dataset_dir, args.mtz_pattern)
    if not dataset_mtz:
        logging.error(
            "No MTZ matching %r found in %s", args.mtz_pattern, dataset_dir
        )
        return 1
    logging.info("Dataset MTZ   : %s", dataset_mtz.name)

    # --- Locate CIF ligand restraints ---
    cif_files = find_cif_files(dataset_dir)
    if cif_files:
        logging.info("CIF files     : %s", ", ".join(c.name for c in cif_files))
    else:
        logging.warning(
            "No *.cif files found in %s — running phenix.ligand_pipeline "
            "without explicit ligand restraints.",
            dataset_dir,
        )

    # --- Determine next Refine_N serial (mirrors XCE's GetSerial()) ---
    serial = get_next_refine_serial(dataset_dir)
    logging.info("Refine serial : %d  (→ Refine_%d/)", serial, serial)

    # --- Step 1: R-free flags are handled internally by phenix.ligand_pipeline ---
    # When --ref_mtz is provided, flags are transferred via
    # xray_data.r_free_flags.file_name; no separate conversion step is needed.

    # --- Step 2: phenix.ligand_pipeline ---
    logging.info("Running phenix.ligand_pipeline...")
    success, pipeline_dir = run_ligand_pipeline(
        ref_pdb=ref_pdb,
        mtz=dataset_mtz,
        seq_file=seq_file,
        cif_files=cif_files,
        dataset_dir=dataset_dir,
        nproc=args.nproc,
        phenix_bin=args.phenix_bin,
        dry_run=args.dry_run,
        ref_mtz=ref_mtz,
        rfree_label=args.rfree_label,
        xray_data_labels=args.xray_data_labels,
        aniso_cutoff=args.aniso_cutoff,
    )
    if not success:
        logging.error("phenix.ligand_pipeline failed. Aborting XCE post-processing.")
        return 1

    # --- Skip XCE output if requested ---
    if args.skip_xce_output:
        logging.info("Skipping XCE output post-processing (--skip_xce_output).")
        logging.info("Pipeline complete for: %s", dataset_dir)
        return 0

    if args.dry_run:
        logging.info(
            "(dry-run) Would create Refine_%d/ structure and root symlinks.",
            serial,
        )
        return 0

    # --- Step 3: locate best output PDB/MTZ ---
    # Overrides (explicit paths) always take priority over auto-detection.
    # If --select_best_model is set we scan all pipeline stages for the lowest
    # Rfree pair; otherwise we use the standard *_final.pdb/*_final.mtz.
    override_pdb = (
        Path(args.ligand_pipeline_output_pdb) if args.ligand_pipeline_output_pdb else None
    )
    override_mtz = (
        Path(args.ligand_pipeline_output_mtz) if args.ligand_pipeline_output_mtz else None
    )
    if override_pdb is not None or override_mtz is not None:
        # Explicit overrides — skip all auto-detection logic
        pipeline_pdb, pipeline_mtz = find_pipeline_outputs(
            pipeline_dir, override_pdb=override_pdb, override_mtz=override_mtz
        )
    elif args.select_best_model:
        logging.info(
            "--select_best_model: scanning pipeline stages for lowest Rfree..."
        )
        pipeline_pdb, pipeline_mtz, sel_rwork, sel_rfree = \
            select_best_pipeline_output(pipeline_dir)
        if pipeline_pdb is None:
            logging.warning(
                "Best-model selection found no valid PDB/MTZ pairs with "
                "REMARK 3 R-values in %s; falling back to *_final.pdb detection.",
                pipeline_dir.name,
            )
            pipeline_pdb, pipeline_mtz = find_pipeline_outputs(pipeline_dir)
    else:
        pipeline_pdb, pipeline_mtz = find_pipeline_outputs(pipeline_dir)

    if not pipeline_pdb:
        logging.error(
            "Could not find output PDB in %s. "
            "Use --ligand_pipeline_output_pdb to specify it explicitly, "
            "or try --select_best_model.",
            pipeline_dir,
        )
        return 1
    if not pipeline_mtz:
        logging.error(
            "Could not find output MTZ in %s. "
            "Use --ligand_pipeline_output_mtz to specify it explicitly, "
            "or try --select_best_model.",
            pipeline_dir,
        )
        return 1

    logging.info("Pipeline PDB  : %s", pipeline_pdb.relative_to(dataset_dir))
    logging.info("Pipeline MTZ  : %s", pipeline_mtz.relative_to(dataset_dir))

    # --- Step 4: XCE-compatible output ---
    logging.info("Creating XCE-compatible Refine_%d/ structure...", serial)
    ok = create_xce_output(
        dataset_dir=dataset_dir,
        xtal_name=xtal_name,
        pipeline_pdb=pipeline_pdb,
        pipeline_mtz=pipeline_mtz,
        rfree_mtz=dataset_mtz,
        serial=serial,
        dry_run=args.dry_run,
    )
    if not ok:
        return 1

    logging.info("=" * 60)
    logging.info("Pipeline complete for: %s", xtal_name)
    logging.info(
        "XCE can read: refine.pdb → Refine_%d/refine_%d.pdb",
        serial, serial,
    )
    logging.info(
        "              refine.mtz → Refine_%d/refine_%d.mtz",
        serial, serial,
    )
    logging.info("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
