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

mkdir -p "$LOGDIR"

# scales=(50.0 60.0 70.0 80.0 90.0 150.0)
# scales=(50.0 60.0 70.0)
# scales=(49.0 99.0 149.0 200.0 300.0 400.0)
# lengths=(0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7)
# lengths=(0.6 0.65 0.7)

# scales=(0 40 50 60)
# scales=(0)
scales=(65 70 75 80)
# scales=(50)
patch_gauss_B_width=13.8504255125
patch_bond_r0=0.45



# lengths=(0.2 0.3 0.4 0.5 0.6 0.7)
# lengths=(0.2 0.4 0.6)
# lengths=(0.2 0.3)
lengths=(0.45)

# n_patch=(1 2 3 4 5)
n_patch=(5)


# for n in "${n_patch[@]}"; do
#   echo "Generating molecule with $n amine patches."
#   python3 generate_molecule_no_dummy.py --amine_spacing $n
#   for scale in "${scales[@]}"; do
#     for length in "${lengths[@]}"; do
#       echo "Launching hbond_energy_scale=$scale and hbond_distance_cutoff=$length"
#       mpirun -np 12 --oversubscribe "$LAMMPS" -in "$INPUT" \
#         -var patch_gauss_A_strength "$scale" \
#         -var patch_bond_r0 "$length" \
#         -var output_file "$LOGDIR/patch_n_${n}_${scale}_${length}.dat"
#     done
#   done
# done


# wait
# echo "All LAMMPS jobs finished."



# Process each energy scale sequentially
for scale in "${scales[@]}"; do
  echo "========================================="
  echo "Processing energy scale: $scale"
  echo "========================================="
  
  # Step 1: Run LAMMPS simulation
  echo "Step 1/4: Running LAMMPS simulation for scale $scale"
  mpirun -np 12 --oversubscribe \
    "$LAMMPS" -in "$INPUT" \
    -var patch_gauss_A_strength "$scale" \
    -var output_file "$LOGDIR/energy_${scale}.lammpstrj"
  
  echo "LAMMPS simulation completed successfully for scale $scale"
  
  # Step 2: Find bonds
  echo "Step 2/4: Finding bonds for scale $scale"
  python3 find_reversible_bonds.py --traj "logs/energy_${scale}.lammpstrj" --output "data/amine_bonds_${scale}.csv"
  
  echo "Bond finding completed successfully for scale $scale"
  
  # Step 3: Find multiplets
  echo "Step 3/4: Finding multiplets for scale $scale"
  output_csv="data/amine_multiplets_${scale}.csv"
  python3 find_multiplets.py --input "data/amine_bonds_${scale}.csv" --output "$output_csv"
  
  echo "Multiplet finding completed successfully for scale $scale"

  echo "Step 4/4: Calculating bond lifetimes for scale $scale"
  python3 calc_bond_lifetime.py --csv_file "data/amine_bonds_${scale}.csv" --output "data/bond_lifetimes_${scale}.csv"
  echo "Bond lifetime calculation completed successfully for scale $scale"

  echo "All processing completed for energy scale $scale"
  echo ""
done

echo "========================================="
echo "All energy scales processed successfully!"
echo "========================================="

