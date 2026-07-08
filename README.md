# HOLE Pore Radius Analysis with MDAnalysis

A Python script for calculating the pore/channel radius profile along a membrane protein channel across a molecular dynamics (MD) trajectory, using the [HOLE program](http://www.holeprogram.org/) via the [`mdahole2`](https://github.com/MDAnalysis/hole2-mdakit) MDAnalysis plugin.

The script runs HOLE frame-by-frame on a trajectory, collects the pore radius profile at each frame, computes the average radius profile (with standard deviation) along the channel (z) axis, and reports the minimum pore radius per frame — a commonly used metric for identifying channel gating/constriction points (e.g. in ion channels, aquaporins, and other membrane pores).

**Author:** Kunal Rai
**Date:** 17 Jan 2025

---

## Table of Contents

- [What This Script Does](#what-this-script-does)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
- [Configuration — What You Need to Edit](#configuration--what-you-need-to-edit)
- [Usage](#usage)
- [Output Files](#output-files)
- [Understanding the Fixes / Patches in the Script](#understanding-the-fixes--patches-in-the-script)
- [Downstream Plotting Examples](#downstream-plotting-examples)
- [Troubleshooting](#troubleshooting)
- [Citing HOLE / MDAnalysis](#citing-hole--mdanalysis)
- [License](#license)

---

## What This Script Does

1. Loads a topology (`.psf`) and trajectory (`.dcd`) into an `MDAnalysis.Universe`.
2. Strips bond information and patches the PDB writer so HOLE (which internally writes temporary PDB frames) doesn't choke on broken `CONECT` records or bad unit-cell dimensions.
3. Runs `HoleAnalysis` from `mdahole2` on a selected atom group (the pore-lining/protein atoms) for every frame in the trajectory.
4. Collects the raw (reaction coordinate, radius) pairs for every frame into a long-format table (`hole_3d_data.txt`) — useful for 3D surface plots of radius vs. z-coordinate vs. frame/time.
5. Bins the pooled data along the channel axis (default: -70 Å to 70 Å, 200 bins) and computes the **mean radius** and **standard deviation** per bin, producing a 1D averaged pore radius profile (`hole_profile_data.txt`).
6. Extracts the **minimum radius per frame** (the tightest constriction point of the channel at each timestep) and saves it (`minimum_radii.txt`).

---

## Requirements

- Python 3.8+
- [HOLE2](http://www.holeprogram.org/) (compiled binary, e.g. `hole2/exe/hole`)
- Python packages:
  - `MDAnalysis`
  - `mdahole2` (MDAnalysis HOLE2 kit / plugin)
  - `numpy`
  - `pandas`
  - `matplotlib`

### Install Python dependencies

\```bash
pip install MDAnalysis mdahole2 numpy pandas matplotlib
\```

### Install HOLE2

HOLE2 is a separate, standalone Fortran/C program and is **not** installed via pip. Download and compile it from the [HOLE program website](http://www.holeprogram.org/) or from the [hole2 GitHub mirror](https://github.com/osmart/hole2), then note the full path to the compiled `hole` executable (e.g. `/home/username/softwares/hole2/exe/hole`). You will need this path for the script (see [Configuration](#configuration--what-you-need-to-edit)).

---

## Input Files

| File | Description |
|---|---|
| `1u19_popc_wi.psf` | Topology file (CHARMM/NAMD-style PSF) describing the system (protein + POPC membrane + water/ions). |
| `md.dcd` | MD trajectory file in DCD format, matching the topology. |

> Any topology/trajectory format combination supported by MDAnalysis (e.g. `.gro`/`.xtc`, `.pdb`/`.dcd`, `.prmtop`/`.nc`) will work — just change the `mda.Universe(...)` call accordingly.

---

## Configuration — What You Need to Edit

Before running, update the following lines in the script for your own system:

1. **Topology and trajectory paths:**
   \```python
   u = mda.Universe('../../1u19_popc_wi.psf', '../md.dcd')
   \```
   Replace with the paths to your own topology (first argument) and trajectory (second argument).

2. **Atom selection for the pore-lining group and HOLE executable path:**
   \```python
   ha = HoleAnalysis(u, select='index 506:5162',
                      executable='/home/krai/softwares/hole2/exe/hole')
   \```
   - `select='index 506:5162'` — an MDAnalysis atom selection string identifying the protein (or channel-forming) atoms HOLE should use to trace the pore. Adjust the index range (or use a more descriptive selection like `protein` or `segid PROA`) to match your system.
   - `executable=...` — the full path to your local HOLE binary.

3. **Approximate channel length (z-range) and bin count:**
   \```python
   coord_range = (-70, 70)  # mention your roughly channel z length
   bins = 200
   \```
   Set `coord_range` to roughly span your channel's length along the pore axis (in Å), centered near 0. Increase/decrease `bins` for finer/coarser resolution of the averaged profile.

4. **(Optional) Random seed:**
   \```python
   ha.run(random_seed=31415)
   \```
   HOLE uses a Monte Carlo simulated-annealing search internally; fixing the seed makes runs reproducible. Change or remove as needed.

---

## Usage

1. Ensure the folder structure matches the relative paths in the script, or edit the paths directly. 
2. Activate the Python environment with the required packages installed.
3. Run the script:
   \```bash
   python hole_radius_analysis.py
   \```
4. Depending on trajectory length and number of frames, HOLE is invoked once per frame, so runtime scales linearly with the number of frames analyzed. For long trajectories, consider first testing on a subset of frames (e.g. `u.trajectory[::10]` or slicing before passing to `HoleAnalysis`) to confirm settings before running the full analysis.

---

## Output Files

| File | Description |
|---|---|
| `hole_3d_data.txt` | Tab-separated long-format table with columns `Frame`, `Coordinate`, `Radius` — one row per (frame, pore-axis point) pair. Suitable for 3D surface/heatmap plots of radius vs. z vs. time. |
| `hole_profile_data.txt` | Tab-separated table with columns `Coordinate`, `Mean_Radius`, `Std_Dev`, `Upper_Bound`, `Lower_Bound` — the trajectory-averaged pore radius profile along the channel axis, binned every `(coord_range[1]-coord_range[0])/bins` Å. |
| `minimum_radii.txt` | Tab-separated table with columns `Frame`, `Minimum_Radius` — the smallest pore radius (i.e., the tightest constriction) found at each frame. |

---

## Understanding the Fixes / Patches in the Script

The script includes two workarounds needed to run HOLE cleanly through MDAnalysis on some systems:

**Fix 1 — Remove bond information**
\```python
u._topology.bonds = Bonds([])
\```
HOLE analysis writes temporary PDB files per frame. If the topology's bond records don't map cleanly to the selected atom subset, MDAnalysis can raise indexing errors when writing `CONECT` records. Clearing the bond table avoids this, since HOLE itself does not need bond connectivity — only atomic coordinates and radii.

**Fix 2 — Patch the PDB writer to drop invalid unit cell dimensions**
\```python
original_write_timestep = PDBmod.PDBWriter._write_timestep
def patched_write_timestep(self, ts, **kwargs):
    ts = ts.copy()
    ts.dimensions = None
    original_write_timestep(self, ts, **kwargs)
PDBmod.PDBWriter._write_timestep = patched_write_timestep
\```
Some trajectories carry malformed or missing unit-cell/box dimensions, which can cause the PDB `CRYST1` record to be written incorrectly (or crash the writer). This monkey-patch temporarily blanks out `ts.dimensions` before each frame is written to a PDB, then restores the original writer behavior (`PDBmod.PDBWriter._write_timestep = original_write_timestep`) immediately after the HOLE run completes, so it doesn't affect any other PDB writing elsewhere in your workflow.

> **Note:** These are pragmatic workarounds for known edge cases in `mdahole2`/MDAnalysis PDB writing. If your system's topology has valid bonds and box dimensions, you can try removing/commenting out these patches — but keep them if you encounter `CONECT`-record or dimension-related errors.

---

## Downstream Plotting Examples

The script imports `matplotlib.pyplot` but doesn't generate plots directly — the intent is for you to load the saved `.txt` files afterward for visualization. Example snippets:

**1D averaged pore radius profile with standard deviation band:**
\```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('hole_profile_data.txt', sep='\t')

plt.figure(figsize=(6, 5))
plt.plot(df['Coordinate'], df['Mean_Radius'], color='navy', label='Mean radius')
plt.fill_between(df['Coordinate'], df['Lower_Bound'], df['Upper_Bound'],
                  alpha=0.3, color='navy', label='±1 SD')
plt.xlabel('Pore axis coordinate (Å)')
plt.ylabel('Radius (Å)')
plt.legend()
plt.tight_layout()
plt.savefig('pore_radius_profile.png', dpi=300)
\```

**Minimum radius vs. frame (gating over time):**
\```python
df_min = pd.read_csv('minimum_radii.txt', sep='\t')
plt.figure(figsize=(6, 4))
plt.plot(df_min['Frame'], df_min['Minimum_Radius'], color='crimson')
plt.xlabel('Frame')
plt.ylabel('Minimum pore radius (Å)')
plt.tight_layout()
plt.savefig('minimum_radius_vs_frame.png', dpi=300)
\```

**3D surface of radius vs. coordinate vs. frame:**
\```python
df3d = pd.read_csv('hole_3d_data.txt', sep='\t')
pivot = df3d.pivot_table(index='Frame', columns='Coordinate', values='Radius')

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(pivot.columns, pivot.index)
ax.plot_surface(X, Y, pivot.values, cmap='viridis')
ax.set_xlabel('Coordinate (Å)')
ax.set_ylabel('Frame')
ax.set_zlabel('Radius (Å)')
plt.tight_layout()
plt.savefig('pore_radius_3d.png', dpi=300)
\```

---

## Troubleshooting

- **`FileNotFoundError` / HOLE executable not found:** Double-check the `executable=` path points to a valid, executable HOLE2 binary (`chmod +x` it if needed).
- **Empty or `NaN`-heavy profile bins:** Your `coord_range` may not match the actual span of the channel in your system. Inspect `hole_3d_data.txt` to see the real range of `Coordinate` values and adjust `coord_range`/`bins` accordingly.
- **HOLE fails on specific frames:** This is often caused by an unsuitable `select=` atom group (e.g. atoms far outside the pore, or too few atoms to define the pore lining). Try using a more targeted selection such as pore-lining residues (`resid ... and name CA`) instead of a broad index range.
- **Errors related to `CONECT` records or `CRYST1`/box dimensions:** These are exactly what Fix 1 and Fix 2 address — make sure both patches remain active if you see these errors.
- **Long runtime:** HOLE is called once per frame; consider analyzing every Nth frame of the trajectory for exploratory analysis before committing to the full trajectory.

---

## License

Add your preferred license here (e.g. MIT, GPL-3.0) so others know how they may reuse this script.
