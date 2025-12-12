#!/bin/bash

: '
Calling input script requires defining the following variables:
- N: number of polymers
- BOX_SIZE: size of the simulation box
- patch_bond_r0: equilibrium bond length between patches and backbone
- patch_gauss_A_strength: strength of Gaussian attraction between amine patches
- patch_gauss_B_width: width of Gaussian attraction between amine patches

- output_traj: set to 1 to output trajectory file
- output_traj_file: name of trajectory file
- output_amine: set to 1 to output amine positions for bond analysis
- output_amine_file: name of amine positions file
- output_msd: set to 1 to output mean squared displacement data
- output_msd_file: name of msd data file
```
example call:
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
Note that the molecule file has already been generated and saved as molecule.data at this point. 
If there are parameters that need to be changed in the molecule file, change them with the script generate_molecule_no_dummy.py before running this script.


'

LAMMPS=../lammps/build/lmp
INPUT=in.amine_polymers_sh
LOGDIR=logs
DATADIR=data

mkdir -p "$LOGDIR"

patch_bond_r0=0.45
patch_gauss_A_strength=50
patch_gauss_B_width=13.8504255125

output_traj=0
output_traj_file="$DATADIR/patch_trajectory.lammpstrj"
output_msd=0
output_msd_file="$LOGDIR/patch_msd.dat"
output_amine=1
output_amine_file="$LOGDIR/patch_locations.dat"
output_press=0
output_press_file="$LOGDIR/press_${patch_gauss_A_strength}_${P}.dat"


echo "Generate long polymer molecule file."
python3 generate_molecule.py --chain_length 100 --amine_spacing 1 --filename data/polymer_long.molecule

echo "Generate short polymer molecule file."
python3 generate_molecule.py --chain_length 28 --amine_spacing 1 --filename data/polymer_short.molecule



spacing=(15 10 5 4 3)
seeds=(100 200 300 400 500)


for j in "${seeds[@]}"; do
        for i in "${spacing[@]}"; do
                echo "Run spacing $i with seed $j"

                echo "Generate... other... polymer molecule file."
                python3 generate_molecule.py --chain_length 30 --amine_spacing $i --filename data/polymer.molecule
                #Single run
                mpirun -np 13 --oversubscribe $LAMMPS -in in.amine_polymers_sh_short \
                        -var patch_bond_r0 $patch_bond_r0 \
                        -var patch_gauss_A_strength $patch_gauss_A_strength \
                        -var patch_gauss_B_width $patch_gauss_B_width \
                        -var output_traj $output_traj \
                        -var output_traj_file $output_traj_file \
                        -var output_msd $output_msd \
                        -var output_msd_file $output_msd_file \
                        -var output_amine $output_amine \
                        -var output_amine_file $output_amine_file \
                        -var output_press $output_press \
                        -var output_press_file "$LOGDIR/press_${patch_gauss_A_strength}_spacing${i}_seed${j}.dat" \
                        -var random_seed $j
        done
done