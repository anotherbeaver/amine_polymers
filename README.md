# amine_polymers
Sketching and simulating generic polymer chains with functional "amines" as sites for reversible dynamic crosslinking.

# Running Polymer Simulation

Broadly, a molecule is defined and provided to the 

## LAMMPS input script
Calling input script requires defining the following variables:

| Variable               | Description                                           |
|------------------------|-------------------------------------------------------|
| N                      | Number of  polymers                                   |
| BOX_SIZE               | Size of the simulation box                            |
| patch_bond_r0          | Equilibrium bond length between patches and backbone  |
| patch_gauss_A_strength | Strength of Gaussian attraction between amine patches |
| patch_gauss_B_width    | Width of Gaussian attraction between amine patches    |
| output_traj            | Set to 1 to output trajectory file                    |
| output_traj_file       | Name of trajectory file                               |
| output_amine           | Set to 1 to output amine positions for bond analysis  |
| output_amine_file      | Name of amine positions file                          |
| output_msd             | Set to 1 to output mean squared displacement data     |
| output_msd_file        | Name of msd data file                                 |
| output_press           | Set to 1 to output pressure_xy per atom data          |
| output_press_file      | Name of pressure data file                            |

Here's an example call:
```
mpirun -np 12 --oversubscribe "$LAMMPS" -in "$INPUT" \
    -var N 300 \
    -var BOX_SIZE 40 \
    -var patch_bond_r0 0.45 \
    -var patch_gauss_A_strength 50 \
    -var patch_gauss_B_width 13.8504255125 \
    -var output_amine "$LOGDIR/patch_n_${n}_${scale}_${length}.dat" \
    -var output_traj 0 \
    -var output_msd 0 \
    -var output_amine 1 \
    -var output_amine_file "$LOGDIR/patch_locations.dat" \

```

In total there are 3 optional output files and 1 default output. The amine positions are tracked and optionally output to determine reversible bonds. The mean squared displacement is calculated by LAMMPS and optionally output to a file. The trajectory of all particles is optionally output, mostly for sanity checks and visualization. Several thermo values are recorded, as well as program output, in a log file. 
## Molecule file
Note that the molecule file has already been generated and saved as molecule.data at this point. 
If there are parameters that need to be changed in the molecule file, change them with the script `generate_molecule.py` before running this script. The following parameters must be defined:

| Variable          | Default Value           | Description                     |
|-------------------|-------------------------|---------------------------------|
| filename          |`'data/polymer.molecule'`| Output filename                 |
| chain_length      | `25`                    | Length of the polymer chain     |
| amine_spacing     | `5`                     | Spacing between amine patches   |
| backbone_to_patch | `0.45`                  | Distance from backbone to patch |

Please note that the filename provided is the one that is currently read by the input file. There is currently no reason to change this as only one input molecule file is needed at any time. Any other filename will NOT be read by the LAMMPS input file.

# Postprocessing

Several post processing options are available. The reversible bonds can be found and bonds involving more than two amines can be collected. The lifetimes can be calculated and plotted. 

## find_reversible_bonds.py
The reversible bonds are calculated with an input position/time file. Please note that the traj file must be generated for this to be calculated, which requires `output_amine` to be true while running the simulation. 

| Variable    | Default Value                      | Description                         |
|-------------|------------------------------------|-------------------------------------|
| traj        | `data/amine_positions.lammpstrj`   | Input LAMMPS trajectory file        |
| output      | `data/amine_bonds.csv`             | Output CSV file for bond data       |
| min_dist    | `0.6`                              | Minimum distance to consider a bond |
| box         | `(40, 40, 40)`                     | Simulation box dimensions (x y z)   |
| max_workers | `4`                                | Maximum number of parallel workers  |

The output is a csv file with columns `["timestep","id1","id2","x1","y1","z1","x2","y2","z2","distance"]`. 

## calc_bond_lifetime.py
The lifetimes of bonds are calculated by gathering consecutive timesteps where bonds exist. 

| Argument    | Default Value                | Description                         |
|-------------|------------------------------|-------------------------------------|
| csv_file    | `data/amine_bonds.csv`       | Input CSV file containing bond data |
| output      | `data/bond_lifetimes.csv`    | Output CSV file for bond lifetimes  |

The output is a csv file with columns `['id1', 'id2', 'lifetime']`. The order of rows does not matter; each entry represents a bond that occured at some point in the simulation.

## find_multiplets.py
A deceivingly named file, this locates and counts the number of neighbours involved in a bond. We ideally want only 1 neighbour in each bond, but there are some parameters where we get more. 

| Argument | Default Value               | Description                    |
|----------|-----------------------------|--------------------------------|
| output   | `data/amine_multiplets.csv` | Output CSV file for multiplets |
| input    | `data/amine_bonds.csv`      | Input CSV file with bond data  |

The output is formatted as a csv with timestep and center (atom id) columns. There are 10 additional columns of atom ids and x, y, z, which allow up to 10 neighbours to be tracked. This is meant to be quite sparse, as we don't generally have much more than 3 neighbours, and those should also be rare. There may be unexpected behaviour if more than 10 neighbours are found. 

TODO: add $G(t)$ post processing file
