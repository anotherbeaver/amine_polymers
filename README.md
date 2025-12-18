# amine_polymers
Sketching and simulating generic polymer chains with functional "amines" as sites for reversible dynamic crosslinking. The name is very much not relevant anymore, I did just name it that on a whim when starting the repo -- other functional groups are also allowed, we won't be picky.

# Table of Contents

- [Installing LAMMPS](#installing-lammps)
- [Running polymer simulation](#running-polymer-simulation)
  - [LAMMPS input script](#lammps-input-script)
  - [Molecule file](#molecule-file)
- [Postprocessing](#postprocessing)
  - [find_reversible_bonds.py](#find_reversible_bondspy)
  - [calc_bond_lifetime.py](#calc_bond_lifetimepy)
  - [find_multiplets.py](#find_multipletspy)
  - [TODO: $G(t) post processing](#todo-gt-post-processing)
- [State of the project](#state-of-the-project)


# Installing LAMMPS
The normal process of installing LAMMPS is perfect, but note that there are a couple flags to set when compiling with `cmake`:

`cmake ../cmake   -D PKG_MOLECULE=ON   -D PKG_EXTRA-PAIR=ON -D PKG_EXTRA-FIX=ON`

the equivalent flags can be found for other installation methods by looking in the LAMMPS documentation

# Running polymer simulation

Broadly, a molecule file is defined and provided to the simulation. Some simulation and bonding parameters are defined, and the simulation will output bonding site trajectories, autocorrelated stress data, or mean squared displacement of all the beads. Previous versions would have saved other data, though those are not hard to add.

Some scripts are provided to postprocess output data. Bonding site trajectories can be used to determine bond lengths, durations, and errors where valency is violated. Please note that the input formats of these scripts are meant to be linked to each other, but that there may still remain errors -- read through all file inputs before running the scripts.

## LAMMPS input script
Calling input scripts require defining the following variables:

| Variable               | Description                                           |
|------------------------|-------------------------------------------------------|
| N                      | Number of  polymers (DEFUNCT, CURRENTLY USE DENSITY)  |
| BOX_SIZE               | Size of the simulation box (DEFUNCT, CURRENTLY USE DENSITY)    |
| patch_bond_r0          | Equilibrium bond length between patches and backbone  |
| patch_gauss_A_strength | Strength of Gaussian attraction between amine patches |
| patch_gauss_B_width    | Width of Gaussian attraction between amine patches    |
| output_traj            | Set to 1 to output all beads trajectory file          |
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

Please note that there are several versions available, including in `\archive`. These represent the history of the input file, and may be run differently than above and may certainly have some errors. Use those with discretion and initial review. Some variables are no longer used and will not do anything to the simulation.


## Molecule file
Note that the molecule file has likely been generated and saved as `molecule.data` at this point. 
If there are parameters that need to be changed in the molecule file, change them in `generate_molecule.py` before running this script. The following parameters must be defined:

| Variable          | Default Value           | Description                     |
|-------------------|-------------------------|---------------------------------|
| filename          |`'data/polymer.molecule'`| Output filename                 |
| chain_length      | `25`                    | Length of the polymer chain     |
| amine_spacing     | `5`                     | Spacing between amine patches   |
| backbone_to_patch | `0.45`                  | Distance from backbone to patch |

Please note that the filename provided is the one that is currently read by the input file. There is currently no reason to change this as only one input molecule file is needed at any time. Any other filename will NOT be read by the LAMMPS input file without changing the code.

# Postprocessing

Several post processing options are available. The reversible bonds can be found and bonds involving more than two patches (violating valency) can be collected. The lifetimes can be calculated and plotted. 

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

## $G(t)$ postprocessing
Although calculating $G(t)$ with stress autocorrelation is a possibility (see in `\archive\`), the best method is [on-the-fly multiple-tau autocorrelation](https://doi.org/10.1063/1.3491098), which is built-in to LAMMPS. This has been used extensively in papers on related simulations, which are worth checking out for some intuition on simulation conditions and length. 

The first part of the $G(t)$ graph is the result of 

# State of the project
Date added: 17. December 2025

At this moment, the code is complete and able to calculate storage and loss moduli data for full simulations with reversible crosslinking. Comparison has been done for the bond density parameter (number of bonds/chain length) for a model system of density $\rho=0.45$, $n=50000$ beads, at $T=1.0$ and Gaussian well depth $A=50$. 

Note that most literature about similar simulations uses $\rho=0.85$, a lower density was used to align with computational resources -- it is certainly possible to increase this density, especially to match the dynamics of these papers. 

Here's some relevant reading for now:

- [Relaxation Dynamics of Entangled Linear Polymer Melts via Molecular Dynamics Simulations (Alireza F. Behbahani & Friederike Schmid)](https://doi.org/10.1021/acs.macromol.4c02168)
   - Follows very closely the process that we do, though addresses chain length and does not have any connection to reversible cross-linking. Very good to read through to understand the simulation + postprocessing pipeline. 
-  [Stress relaxation in tunable gels (Raffaelli, C., & Ellenbroek, W. G.)](https://doi.org/10.1039/D1SM00091H)
   - This is where we found the Gaussian potential method originally. Note the paper erroneously writes it as an exponential.
- [Three-body potential for simulating bond swaps in molecular dynamics (Sciortino, F.)](https://doi.org/10.1140/epje/i2017-11496-5)
   - Another method explored, using a three body potential that prevents valency violation. The lack of bonding energy and possible repulsion in the three body case was not good though.
- [Equilibrium and non-equilibrium molecular dynamics approaches for the linear viscoelasticity of polymer melts (Adeyemi, O., Zhu, S., Xi, L.)](https://doi.org/10.1063/5.0090540)
   - Details how $G(t), G'(T), G''(t)$ are retrived with equilibrium molecular dynamics as is done in this project. Discusses the difference and equivalence with non-equilibrium methods.

On the software side, quite a few data files that were evaluated through Jupyter notebooks over the last few months are not around anymore -- either they can be regenerated or at this stage are no longer needed, if they were addressing old problems. As such, the state of most of the notebooks is left as is, and are not intended to be run again.

## What's next?
