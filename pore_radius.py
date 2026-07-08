###########Python code for calculating radius using HOLE program and MDanalysis###########
###########Author: Kunal Rai: Date 17 Jan 2025###############
import MDAnalysis as mda
from mdahole2.analysis import HoleAnalysis
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd
from MDAnalysis.core.topologyattrs import Bonds
import MDAnalysis.coordinates.PDB as PDBmod

warnings.filterwarnings('ignore')

# Load trajectory
u = mda.Universe('../../1u19_popc_wi.psf', '../md.dcd')   # Mention you input files here; Topology file first and then your trajectory file

# Fix 1: Remove bond information to prevent CONECT record indexing error
u._topology.bonds = Bonds([])

# Fix 2: Patch PDB writer to strip broken unit cell info before writing
original_write_timestep = PDBmod.PDBWriter._write_timestep

def patched_write_timestep(self, ts, **kwargs):
    ts = ts.copy()
    ts.dimensions = None
    original_write_timestep(self, ts, **kwargs)

PDBmod.PDBWriter._write_timestep = patched_write_timestep

# Initialize HOLE analysis
ha = HoleAnalysis(u, select='index 506:5162',
                        executable='/home/krai/softwares/hole2/exe/hole') #Mention here you hole programm path and you selection of protein in index

# Run analysis
ha.run(random_seed=31415)

# Restore original PDB writer after run
PDBmod.PDBWriter._write_timestep = original_write_timestep

# Get flat data for profile analysis
flat_data = ha.gather(flat=True)
coord_range = (-70, 70) #mention your roughly channel z length 
bins = 200

# Create bins for coordinates
coord_bins = np.linspace(coord_range[0], coord_range[1], bins+1)
midpoints = 0.5 * (coord_bins[1:] + coord_bins[:-1])

# Save 3D plot data
profiles = ha.results.profiles
rxn_coords = []
radii_list = []
frames = []

# Collect data for 3D plot
for frame in range(len(profiles)):
    profile = profiles[frame]
    for point in profile:
        rxn_coords.append(point[0])
        radii_list.append(point[1])
        frames.append(frame)

df_3d = pd.DataFrame({
    'Frame': frames,
    'Coordinate': rxn_coords,
    'Radius': radii_list
})

# Save with tab separator
df_3d.to_csv('hole_3d_data.txt', sep='\t', index=False)

# Calculate binned statistics for profile data
df = pd.DataFrame({
    'Coordinate': rxn_coords,
    'Radius': radii_list
})

# Group by coordinate bins and calculate statistics
binned_stats = []
for i in range(len(coord_bins)-1):
    mask = (df['Coordinate'] >= coord_bins[i]) & (df['Coordinate'] < coord_bins[i+1])
    radii_in_bin = df.loc[mask, 'Radius']
    if len(radii_in_bin) > 0:
        mean_radius = radii_in_bin.mean()
        std_radius = radii_in_bin.std()
    else:
        mean_radius = np.nan
        std_radius = np.nan
    binned_stats.append({
        'Coordinate': midpoints[i],
        'Mean_Radius': mean_radius,
        'Std_Dev': std_radius
    })

# Create and save profile data
df_profile = pd.DataFrame(binned_stats)
df_profile['Upper_Bound'] = df_profile['Mean_Radius'] + df_profile['Std_Dev']
df_profile['Lower_Bound'] = df_profile['Mean_Radius'] - df_profile['Std_Dev']

# Save with tab separator
df_profile.to_csv('hole_profile_data.txt', sep='\t', index=False)

# Calculate and save minimum radii
min_radii = ha.min_radius()
with open('minimum_radii.txt', 'w') as f:
    f.write("Frame\tMinimum_Radius\n")
    for frame, min_radius in min_radii:
        f.write(f"{int(frame)}\t{min_radius:.3f}\n")
