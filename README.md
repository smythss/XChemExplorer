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
Ubuntu:
```
sudo apt-get install sshfs
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


## Installation on WEHI Milton HPC

XChemExplorer can be run on the [WEHI Milton HPC](https://wehieduau.sharepoint.com/sites/rc2/SitePages/using-milton.aspx) using an [Apptainer](https://apptainer.org/) (formerly Singularity) container.  CCP4 and PHENIX are **not** bundled in the container image; they are bind-mounted from the HPC filesystem at runtime.

### Prerequisites

| Requirement | Notes |
|---|---|
| **Apptainer / Singularity** | Available on Milton via `module load apptainer` (or equivalent). |
| **CCP4 8.0+** | Must be installed on the HPC filesystem. The default path used by the launch script is `/stornext/System/data/apps/ccp4/ccp4-8.0`. |
| **PHENIX** *(optional, recommended)* | Default path: `/stornext/System/data/apps/phenix/phenix-1.21`. |
| **MOGUL** *(optional)* | Requires a CCDC licence. Set `BDG_TOOL_MOGUL` to the full path of the `mogul` executable. |
| **qsub** | On Milton (SLURM), install the `slurm-torque` compatibility package so that a `qsub` wrapper is present in `PATH`. |

### 1. Clone the repository

```bash
git clone https://github.com/xchem/XChemExplorer
cd XChemExplorer/
```

### 2. Build the container image

Run the following command **once** from the repository root.  This requires Apptainer to be available (load the module first if needed):

```bash
module load apptainer   # adjust module name to match your Milton environment
apptainer build xce_wehi.sif Singularity.def
```

The resulting `xce_wehi.sif` image bundles the XChemExplorer source and all required Python/GUI libraries.  It does **not** include CCP4 or PHENIX.

### 3. Configure software paths (if necessary)

The `XChemExplorer_wehi` launch script uses the following defaults:

| Variable | Default | Purpose |
|---|---|---|
| `CCP4` | `/stornext/System/data/apps/ccp4/ccp4-8.0` | CCP4 installation directory |
| `PHENIX` | `/stornext/System/data/apps/phenix/phenix-1.21` | PHENIX installation directory |
| `BDG_TOOL_MOGUL` | *(unset)* | Full path to the `mogul` executable |

If the software is installed elsewhere, export the relevant environment variable **before** running the launch script, for example:

```bash
export CCP4=/path/to/ccp4
export PHENIX=/path/to/phenix           # optional
export BDG_TOOL_MOGUL=/path/to/mogul    # optional
```

### 4. Launch XChemExplorer

```bash
./XChemExplorer_wehi
```

The script will:
1. Verify that Apptainer/Singularity and the container image (`xce_wehi.sif`) are present.
2. Bind-mount CCP4 (and optionally PHENIX/MOGUL) from the host into the container.
3. Start XChemExplorer via `ccp4-python -m xce` inside the container.

You can add an alias to your shell profile for convenience:

```bash
alias xce="<full_path_to_repository>/XChemExplorer_wehi"
```

### Notes

- **Display / X11**: You need a working X11 session.  On Milton, use an interactive job with X11 forwarding (e.g. `sinteractive --x11`) or an NX/remote-desktop session.
- **SLURM / qsub**: XCE submits processing jobs via `qsub`.  On Milton this is provided by the `slurm-torque` compatibility layer.  Confirm `qsub` is in your `PATH` before launching XCE.
- **Rebuilding the image**: Re-run `apptainer build xce_wehi.sif Singularity.def` whenever you update the XChemExplorer source.

---

## License

XChemExplorer is licensed under the MIT license.
