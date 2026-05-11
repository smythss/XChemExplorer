[![Build Status](https://travis-ci.org/xchem/XChemExplorer.svg?branch=master)](https://travis-ci.org/xchem/XChemExplorer)
<a href="https://codeclimate.com/github/xchem/XChemExplorer/"><img src="https://codeclimate.com/github/xchem/XChemExplorer/badges/gpa.svg" /></a>
<a href="https://codeclimate.com/github/xchem/XChemExplorer/"><img src="https://codeclimate.com/github/xchem/XChemExplorer/badges/issue_count.svg" /></a>
[![HitCount](http://hits.dwyl.io/xchem/XChemExplorer.svg)](http://hits.dwyl.io/xchem/XChemExplorer)

# XChemExplorer (XCE)
<i> "The XChemExplorer graphical workflow tool for routine or large-scale protein-ligand structure determination." </i>
<p>Acta Crystallogr D Struct Biol. 2017 Mar 1;73(Pt 3):267-278. (https://doi.org/10.1107/S2059798316020234)</p> 

## Scope 

XChemExplorer (XCE) is a data-management and workflow tool to support large-scale simultaneous analysis of protein-ligand complexes during structure-based ligand discovery (SBLD). 

The user interfaces of established crystallographic software packages such as CCP4 [Winn et al. (2011), Acta Cryst. D67, 235-242] or PHENIX [Adams et al. (2010), Acta Cryst. D66, 213-221] have entrenched the paradigm that a 'project' is concerned with solving one structure. This does not hold for SBLD, where many almost identical structures need to be solved and analysed quickly in one batch of work. Functionality to track progress and annotate structures is essential. 

XCE provides an intuitive graphical user interface which guides the user from data processing, initial map calculation, ligand identification and refinement up until data dissemination. It provides multiple entry points depending on the need of each project, enables batch processing of multiple data sets and records metadata, progress and annotations in an SQLite database. 

## Requirements
Operating Systems:
- Linux
- Mac OSX

Prerequisites:
- CCP4 version 7.0 (or higher)
- PHENIX (optional, but recommended)

<b>Please note:</b> The recommended installation process, described below, includes an install of CCP4. This is so that we know that our code works with the version of CCP4, rather than your system version, making user support much easier

### Windows users:
Potential solutions:

a) Partition your hard drive and install a light-weight version of linux, such as Ubuntu: https://www.ubuntu.com/download/desktop

b) Install Ubuntu (or other) on a USB drive, and boot from that: https://tutorials.ubuntu.com/tutorial/tutorial-create-a-usb-stick-on-windows#0 (NB: This will require a huge USB stick or external hard drive - ccp4 is a large package)

c) VirtualBox - emulate a linux environment on your Windows desktop: https://www.virtualbox.org

## Installation
1. Clone the github repository onto your machine with:
```
git clone https://github.com/xchem/XChemExplorer
```

2. Change directory into the repository:
```
cd XChemExplorer/
```

3. Run the test_build.sh script (this is currently only included for bash, but you can view the steps within the script and modify it as you please for other shells):
```
./test_build.sh
```

4. To execute, run the XChemExplorer_local.sh script
```
./XChemExplorer_local.sh
```

(<i>We recommend you add an alias to your bash profile to do this</i>)
```
alias xce="<full_path_to_local_git_repository>/XChemExplorer_local.sh"
```

## Remote running of XCE and model building on Diamond’s filesystem
NX latency is a killer when trying to do anything graphical, like model building during pandda.inspect or during refinement.  So instead use Filesystem in Userspace (FUSE), which mounts Diamond’s disk from your computer using the ssh protocol.  All you need is root privileges, or bribe your IT team.

1. Setup fuse (requires admin rights)
Install FUSE (http://github.com/libfuse/libfuse)

Rocky Linux:
```
sudo dnf install fuse-sshfs
```
See https://gist.github.com/cstroe/e83681e3510b43e3f618 for details.  FUSE is also available for Mac.  Avoid building from source unless you really have to.

Create a /dls mount point and give your user ownership:
```
sudo mkdir /dls
sudo chown <yourUID>:<yourGID) /dls
```

2. Mount the Diamond filesystem

When you need it, run this from your own user account:
```
sshfs -o reconnect <your_fed_id>@nx.diamond.ac.uk:/dls /dls
e.g. sshfs -o reconnect zqr16691@nx.diamond.ac.uk:/dls /dls
```
(Recommended):  make your ssh client to keep the link alive by editing ~/.ssh/config and adding these lines:
```
ServerAliveInterval 15
ServerAliveCountMax 3
```
(The “reconnect” option is meant to do this too, but it slows down and eventually drops it anyway.)

3. Fire up XCE
```
cd <your-labxchem-visit-dir>
```
e.g. cd /dls/labxchem/data/2017/lb18145-3/processing
```
/dls/science/groups/i04-1/software/XChem/xce
```
(See Note 1 below.)

4. When you’re done, unmount the file system
```
fusermount -u /dls
```

### Notes:
1. Launching XCE this way avoids insane installation issues, because it runs the software directly off Diamond’s drives, i.e. it’s exactly what you use if physically at Diamond or through NX.  This is therefore supported.
The trade-offs are currently (but we’re trying to fix):

a. It takes a very long time (up to 10 minutes in the UK and up to an hour(!) in North America) to launch XCE and pandda.inspect.  After that they are fully responsive.

b. The xce-coot plugin does not always work, due to a library mismatch.

2. IF YOU KNOW WHAT YOU’RE DOING and really hate (!) the lag on startup:

a. install XCE locally: https://github.com/xchem/XChemExplorer

b. update your PanDDA installation (the ccp4 version is outdated):

```
ccp4-python -m pip uninstall panddas
ccp4-python -m pip install pip –upgrade
ccp4-python -m pip install numpy –upgrade
ccp4-python -m pip install panddas
```
(Note, this is not supported – you really are on your own!  Recommended only for geeks.)

3. If your IT team doesn’t want to help you, try setting up a virtual machine – then you have root permissions and can do whatever you like.  Performance ought to be okay.


## Running XCE on a Diamond OnDemand GPU Desktop session

XCE can be launched directly from a [Diamond OnDemand](https://ondemand.diamond.ac.uk/) GPU Desktop interactive session.  OnDemand desktop sessions use TurboVNC to deliver the graphical environment, which does not support the X11 MIT-SHM shared-memory extension.  The `XChemExplorer_dls` launch script sets `QT_X11_NO_MITSHM=1` to handle this automatically.

1. Start a **GPU Desktop** interactive session from the OnDemand portal.
2. Open a terminal inside the desktop.
3. Navigate to your visit data directory and launch XCE:
```bash
cd /dls/labxchem/data/<year>/<visit>/processing
/dls/science/groups/i04-1/software/XChem/XChemExplorer_dls
```

No additional configuration is required; the `QT_X11_NO_MITSHM=1` environment variable is set by the launch script so that the Qt GUI renders correctly inside the VNC-backed desktop.

## Installation on WEHI Milton HPC

XChemExplorer can be run on the [WEHI Milton HPC](https://wehieduau.sharepoint.com/sites/rc2/SitePages/using-milton.aspx) using an [Apptainer](https://apptainer.org/) (formerly Singularity) container.  CCP4 and PHENIX are **not** bundled in the container image; they are bind-mounted from the HPC filesystem at runtime.

### Prerequisites

| Requirement | Notes |
|---|---|
| **Apptainer / Singularity** | Available on Milton via `module load apptainer` (or equivalent). |
| **CCP4 7.1+** | Default path used by the launch script: `/stornext/System/data/software/rhel/9/base/structbio/ccp4/ccp4-7.1/` |
| **PHENIX** *(optional, recommended)* | Default path: `/stornext/System/data/software/rhel/9/base/structbio/phenix/1.21.2-5419/` |
| **MOGUL** *(optional)* | Requires a CCDC licence. Set `BDG_TOOL_MOGUL` to the full path of the `mogul` executable. |
| **sbatch / qsub** | Milton uses SLURM natively. `sbatch` is available on all login and compute nodes. |

### 1. Clone the repository

```bash
git clone https://github.com/xchem/XChemExplorer
cd XChemExplorer/
```

### 2. Build the container image

Run the following command **once** from the repository root.  This requires Apptainer to be available:

```bash
module load apptainer   # adjust module name to match your Milton environment
apptainer build xce_wehi.sif Singularity.def
```

The resulting `xce_wehi.sif` image bundles the XChemExplorer source and all required Python/GUI libraries.  It does **not** include CCP4 or PHENIX.  Rebuild the image after any update to the XChemExplorer source tree.

### 3. Configure software paths (if necessary)

The `XChemExplorer_wehi` launch script uses the following defaults:

| Variable | Default | Purpose |
|---|---|---|
| `CCP4` | `/stornext/System/data/software/rhel/9/base/structbio/ccp4/ccp4-7.1/` | CCP4 installation directory |
| `PHENIX` | `/stornext/System/data/software/rhel/9/base/structbio/phenix/1.21.2-5419/` | PHENIX installation directory |
| `BDG_TOOL_MOGUL` | *(unset)* | Full path to the `mogul` executable |

If the software is installed elsewhere, export the relevant environment variable **before** running the launch script, for example:

```bash
export CCP4=/path/to/ccp4
export PHENIX=/path/to/phenix           # optional
export BDG_TOOL_MOGUL=/path/to/mogul    # optional
```

### 4. Launch XChemExplorer

XCE requires a graphical display.  Start an interactive session with X11 forwarding before launching it:

```bash
sinteractive --x11 --mem=8G --cpus-per-task=2
```

Then launch XCE:

```bash
/path/to/XChemExplorer/XChemExplorer_wehi
```

Add an alias to your `~/.bashrc` for convenience:

```bash
alias xce="/path/to/XChemExplorer/XChemExplorer_wehi"
```

The script will:
1. Verify that Apptainer/Singularity and the container image (`xce_wehi.sif`) are present.
2. Bind-mount CCP4 (and optionally PHENIX/MOGUL), the SLURM binaries, and `/vast/` into the container.
3. Start XChemExplorer via `ccp4-python -m xce` inside the container.

---

## WEHI Workflow Tutorial: From Synchrotron Data to XCE Database

This section walks through the complete workflow for taking raw autoprocessing output from the Australian Synchrotron (MX1 or MX3 beamlines) and loading it into XChemExplorer on Milton.

### Workflow overview

```
Synchrotron autoprocessing output
          │
          ▼
 prepare_mx1_for_xce.sh          (MX1 beamline)
 prepare_fastdp_for_xce.sh       (MX3 beamline, fast_dp)
          │
          │  creates XCE-compatible directory tree with symlinks
          ▼
 <beamline_dir>/processed/<target>/<crystal>/<run>/<proc_code>/output/
          │
          ▼
 populate_xce_db.py
          │
          │  reads aimless.log files, writes stats into SQLite DB,
          │  copies log/MTZ files into project autoprocessing dirs
          ▼
 <project_dir>/<project_name>.sqlite
          │
          ▼
 XChemExplorer GUI
   Settings tab  →  point at project dir and beamline dir
   Datasets tab  →  select target, click "Select Best Autoprocessing Result"
   Maps tab      →  run Dimple / PIPEDREAM to calculate electron density maps
   PanDDA tab    →  run pandda.analyse to find ligand binding events
   Reference tab →  build reference model in Coot, run initial refinement
```

---

### Step 1 — Prepare your directory structure

Choose a top-level **project directory** on `/vast/scratch/` for all XCE output (maps, refinements, database):

```bash
mkdir -p /vast/scratch/users/$USER/my_project
```

Choose a separate **beamline directory** that will hold the reorganised autoprocessing output.  This can be the same directory or a subdirectory of it:

```bash
mkdir -p /vast/scratch/users/$USER/my_project/beamline
```

---

### Step 2 — Prepare MX3 (fast_dp) data

If your data were processed by `fast_dp` on the MX3 beamline, use `prepare_fastdp_for_xce.sh`.

**Expected input directory layout:**

```
<SOURCE_DIR>/
  fast_dp_results_mpc-0019-1-sn02731452_1/
    aimless.log
    fast_dp.mtz
    ...
  fast_dp_results_mpc-0019-2-sn02731453_1/
    ...
```

The crystal name is extracted as the first three dash-separated fields (e.g. `mpc-0019-1`).  The run number is the integer after the final underscore.

**Run the script:**

```bash
bash /path/to/XChemExplorer/prepare_fastdp_for_xce.sh
```

The script will prompt for:

| Prompt | Example value |
|---|---|
| Source directory (contains `fast_dp_results_*`) | `/stornext/Home/data/allstaff/smyth.s/mx3_visit/` |
| Beamline directory (XCE Data Collection Directory) | `/vast/scratch/users/$USER/my_project/beamline` |
| Target name | `Bax` |
| SMILES library CSV *(optional)* | `/path/to/LifeChem_library.csv` |
| Compound distribution CSV *(optional)* | `/path/to/MX3_distribution.csv` |

**Output layout written to `<beamline_dir>/processed/<target>/`:**

```
processed/Bax/
  mpc-0019-1/
    1/
      fast_dp/
        output/
          LogFiles/
            aimless.log        ← symlink to original
          DataFiles/
            mpc-0019-1.free.mtz ← symlink to fast_dp.mtz
    mpc-0019-1.smi             ← SMILES string (if CSVs provided)
    mpc-0019-1.cmpd            ← compound SN code (if CSVs provided)
  mpc-0019-2/
    ...
```

**SMILES / compound data:**  If you supply both a SMILES library CSV (columns `CA Sample Number`, `QCL_SMILES`) and a compound distribution CSV (columns `Name`, `Compound`) the script will write per-crystal `.smi` and `.cmpd` files.  These are read automatically by `populate_xce_db.py` — you do not need to re-supply the CSVs later.

---

### Step 3 — Prepare MX1 data

If your data were processed by the MX1 autoprocessing pipeline, use `prepare_mx1_for_xce.sh`.

**Expected input directory layout:**

```
<SOURCE_DIR>/
  MPC0022_MB6_0084_20260313-163230_process/
    aimless.log
    MPC0022_MB6_0084_20260313-163230_process_aimless.mtz
    ...
  MPC0022_MB6_0084_20260313-205140_retrigger/
    aimless.log
    ...
  MPC0022_MB6_0085_20260313-163400_process/
    ...
```

The crystal name is the first three underscore-delimited fields (e.g. `MPC0022_MB6_0084`).  Multiple processing directories for the same crystal (e.g. `_process`, `_retrigger`) are each assigned an incrementing run number.

**Run the script:**

```bash
bash /path/to/XChemExplorer/prepare_mx1_for_xce.sh
```

The script will prompt for:

| Prompt | Example value |
|---|---|
| Source directory (contains MX1 output folders) | `/stornext/Home/data/allstaff/smyth.s/mx1_visit/` |
| Beamline directory (XCE Data Collection Directory) | `/vast/scratch/users/$USER/my_project/beamline` |
| Target name | `Bax` |
| SMILES library CSV *(optional)* | `/path/to/LifeChem_library.csv` |
| Compound distribution CSV *(optional)* | `/path/to/MX1_distribution.csv` |

The output layout is identical to the MX3 case above, except the processing code directory is named `mx1_process` instead of `fast_dp`.

---

### Step 4 — Initialise the XCE project and create the database

XCE creates its SQLite database automatically when you first point it at a new project directory.

1. Launch XCE:

   ```bash
   xce
   ```

2. In the **Settings** tab:
   - Set **Project Directory** to `/vast/scratch/users/$USER/my_project`
   - Set **Data Collection Directory** to `/vast/scratch/users/$USER/my_project/beamline`
   - Set **Datasets Imported from Beamline** → uncheck **Read Agamemnon data structure** (this is a Diamond-specific path format; leave it unchecked for Australian Synchrotron data)
   - Click **Save Settings**

3. XCE will create `<project_dir>/<project_name>.sqlite`.  Close XCE.

---

### Step 5 — Populate the database with autoprocessing statistics

The standard XCE "Get New Results from Autoprocessing" button in the Datasets tab uses Diamond Light Source path conventions and will not find data in the layout created by the prepare scripts.  Use `populate_xce_db.py` to bypass this and directly populate the database:

```bash
python3 /path/to/XChemExplorer/populate_xce_db.py
```

The script will prompt for:

| Prompt | Example value |
|---|---|
| XCE `.sqlite` file path | `/vast/scratch/users/$USER/my_project/my_project.sqlite` |
| `processed/<target>/` directory | `/vast/scratch/users/$USER/my_project/beamline/processed/Bax` |
| XCE Project Directory | `/vast/scratch/users/$USER/my_project` |
| Target / ProteinName | `Bax` |
| SMILES library CSV *(blank if `.smi` files already present)* | *(leave blank)* |
| Distribution CSV *(blank if `.cmpd` files already present)* | *(leave blank)* |
| Beamline `mx1`/`mx3` *(only needed if CSVs above given)* | *(leave blank)* |

For each crystal the script will:
- Parse the `aimless.log` to extract resolution, Rmerge, completeness, space group, and unit cell.
- Compute a data-quality score and alert colour (green/amber/red).
- Copy the log and MTZ into `<project_dir>/<crystal>/autoprocessing/` (matching the path XCE expects when it reads the database back).
- `INSERT` or `UPDATE` the `collectionTable` row.
- Populate `mainTable` with compound code and SMILES if `.smi`/`.cmpd` files are present.

Sample output:

```
  INSERT  mpc-0019-1  run=1  fast_dp  res=1.85  SG=P 21 21 21  cmpd=SN02731452 (from file)
  INSERT  mpc-0019-2  run=1  fast_dp  res=2.10  SG=P 21 21 21  cmpd=SN02731453 (from file)
  ...
Done: 48 inserted, 0 updated, 2 skipped.
```

---

### Step 6 — View datasets in XCE and select best autoprocessing result

1. Launch XCE and go to the **Datasets** tab.
2. Select your target from the dropdown (e.g. `Bax`).
3. The table should now show all crystals with their processing statistics and quality alerts.
4. Click **Select Best Autoprocessing Result** (the **Run** button, **not** the Status button).  This assigns one processing result per crystal to be used for map calculation.

> **Important:** Do **not** click the Status buttons at this stage.  These attempt to contact a DLS-specific cluster bastion node (`CLUSTER_BASTION`) which does not exist at WEHI and will raise a `socket.gaierror`.

---

### Step 7 — Calculate initial electron density maps

With autoprocessing results assigned, XCE can run Dimple (or PIPEDREAM) to calculate difference Fourier maps and perform initial refinement against the reference model.

1. In the **Maps** tab, ensure a reference PDB is set (Settings tab → **Reference PDB** field).
2. Select all crystals (Ctrl+A or the select-all checkbox).
3. Click **Run Dimple** (or **Run PIPEDREAM**).

XCE will generate and submit SLURM batch scripts for each crystal.  Jobs are submitted to the `regular` partition by default.  Monitor progress with:

```bash
squeue -u $USER
```

When Dimple finishes, each crystal directory will contain an `autoprocessing/` subdirectory with the initial MTZ and map files used by subsequent steps.

---

### Step 8 — Run PanDDA analysis

PanDDA statistically analyses the ensemble of electron density maps to identify ligand binding events above the ground state.

1. Go to the **PanDDA** tab.
2. Set the **PanDDA Input Data Directory** to the project directory.
3. Set the **PanDDA Output Directory** (e.g. `<project_dir>/pandda`).
4. Choose which datasets to include (use the filter to select datasets with good map quality).
5. Click **Run PanDDA**.

XCE generates a `pandda.sh` SLURM batch script and submits it.  A typical PanDDA run over 50–100 datasets takes 2–8 hours depending on resolution and ASU size.

After completion, PanDDA outputs include:
- `pandda/analyses/pandda_analyse_events.csv` — table of binding events with Z-scores and BDC values
- `pandda/processed_datasets/<crystal>/` — per-dataset event maps and modelled structures
- `pandda/processed_datasets/<crystal>/<crystal>-pandda-input.mtz` — input MTZ
- `pandda/processed_datasets/<crystal>/<crystal>-z_map.native.mtz` — Z-map
- `pandda/processed_datasets/<crystal>/<crystal>-event_<N>_map.native.mtz` — event map(s)

---

### Step 9 — Build and inspect ligand models with pandda.inspect

After PanDDA analysis, use `pandda.inspect` to review events, build ligand models, and annotate results.

1. In the **PanDDA** tab, click **Open pandda.inspect**.
2. For each event, pandda.inspect opens Coot pre-loaded with the event map.
3. Fit the ligand using Coot's ligand fitting tools or drag the pre-placed ligand into density.
4. Mark each event as a **Hit**, **Interesting**, or **Reject** using the buttons in the pandda.inspect panel.
5. Save and move to the next event.

---

### Step 10 — Select the ground state reference model

Before batch refinement, XCE needs a reference model (apo / ground state structure):

1. In the **Reference** tab (or via the menu), identify a high-resolution apo dataset.
2. Run `select_ground_state_dataset.py` (or manually create symlinks) to populate `<project_dir>/reference/`:
   - `reference.pdb` → `<dataset>-aligned-structure.pdb` (the pandda-aligned version)
   - `reference.mtz` → `<dataset>-pandda-input.mtz`
3. In the **Settings** tab, set **Reference PDB** to `<project_dir>/reference/reference.pdb`.

---

### Step 11 — Build the reference model and run initial refinement

With a reference model set, use the XCE Coot reference model builder to refine the ground state structure before batch refinement of all datasets.

1. In the **Reference** tab, click **Open Coot Reference Model Builder**.
2. In the Coot GUI, make any geometry corrections to the reference model.
3. Click **REFINE** to submit an initial restrained refinement job to SLURM.  The job is submitted directly via `sbatch` to the `regular` partition — no SSH or REST API token is required.
4. Once the refinement completes (check with `squeue -u $USER`), click **Re-load Refined Model** to review the result in Coot.
5. If satisfied, click **Accept** to register the refined reference model.

---

### Step 12 — Batch refinement

With the reference model accepted, return to the **Datasets** tab and run batch refinement:

1. Select all hit datasets (those with pandda events marked as **Hit**).
2. Click **Run Refinement**.

XCE submits one SLURM job per crystal.  Each job runs `refmac5` (with or without PHENIX) against the event map, using the reference model as a starting point and the fitted ligand coordinates from pandda.inspect.

Monitor jobs:

```bash
squeue -u $USER
```

When refinement finishes, the results appear in the **Refinement** tab, showing Rfree, Rwork, and geometry statistics for each crystal.

---

### WEHI-specific notes and known limitations

| Topic | Notes |
|---|---|
| **CLUSTER_BASTION** | XCE contains hard-coded DLS constants (`CLUSTER_BASTION`, `CLUSTER_HOST`) used for SLURM REST API token authentication. At WEHI these are undefined and any button that calls `get_token()` will raise `socket.gaierror`. The REFINE button in the reference model builder has been patched to use direct `sbatch` instead. Avoid clicking **Status** buttons in the Datasets tab. |
| **PanDDA version** | The bundled `pandda.analyse` (PanDDA 1.x) requires CCP4 7.0's Python 2.7 environment. CCP4 7.1 ships a newer Python 2.7 with a different pandas ABI. The `pandda.sh` script sources CCP4 7.0 for the `pandda.analyse` call and restores CCP4 7.1 afterwards. |
| **PyMOL** | `pandda.analyse` calls PyMOL to generate event-map images. Load it in the `pandda.sh` header via `module load pymol` or set `BDG_TOOL_PYMOL` to the full path of the PyMOL executable. |
| **`/vast/` filesystem** | All scratch I/O should use `/vast/scratch/users/$USER/`. The launch script bind-mounts `/vast/` into the container automatically. |
| **Display / X11** | You need a working X11 session. Use `sinteractive --x11` or an NX / remote-desktop session. |
| **Rebuilding the image** | Re-run `apptainer build xce_wehi.sif Singularity.def` whenever you update the XChemExplorer source. |

---

## License

XChemExplorer is licensed under the MIT license.
