#!/bin/bash
# Original location: <PROJECT_ROOT>/run.sh
# Old LAMMPS script for running on the box upstairs -- lots of 
# legacy and leftover stuff in here, but might be good to read.
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

# N=131
# P=87
# N=2000 # total number of polymers
# P=0 # number of long polymers
# BOX_SIZE=43.088
# patch_bond_r0=0.45
# patch_gauss_A_strength=0
# patch_gauss_B_width=13.8504255125
# temp=1.0


N=30 # total number of polymers
P=0 # number of long polymers
BOX_SIZE=7.06
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

# P_values=(0 50 100 150 200 250 300)
P_values=(0)
# patch_A_strengths=(0 40 50 60 100 150)
patch_A_strengths=(60)


echo "Generate long polymer molecule file."
python3 generate_molecule.py --chain_length 100 --amine_spacing 1 --filename data/polymer_long.molecule

echo "Generate short polymer molecule file."
python3 generate_molecule.py --chain_length 28 --amine_spacing 1 --filename data/polymer_short.molecule



# for scale in "${patch_A_strengths[@]}"; do
#   for P in "${P_values[@]}"; do
# #     N=$((300))
#     echo "Running LAMMPS with N=$N, P=$P, patch_gauss_A_strength=$scale"
    
#     mpirun -np 12 --oversubscribe ../lammps/build/lmp -in in.amine_polymers_sh \
#           -var N $N \
#           -var P $P \
#           -var BOX_SIZE $BOX_SIZE \
#           -var patch_bond_r0 $patch_bond_r0 \
#           -var patch_gauss_A_strength $scale \
#           -var patch_gauss_B_width $patch_gauss_B_width \
#           -var output_traj 0 \
#           -var output_traj_file "$LOGDIR/patch_trajectory_P${P}.lammpstrj" \
#           -var output_msd 0 \
#           -var output_msd_file "$LOGDIR/patch_msd_P${P}.dat" \
#           -var output_amine 0 \
#           -var output_amine_file "$LOGDIR/patch_locations_P${P}_scale${scale}.dat" \
#           -var output_press 1 \
#           -var output_press_file "$LOGDIR/press_${scale}_${P}_long.dat" \
#           -var random_seed $((scale + 1))
#   done
# done

# spacing=(1 2 3 5 10 15 20 29)
# seeds=(200 300 400 500)

spacing=(15 10 5 4 3)
seeds=(100 200 300 400 500)

# spacing=(15)
# seeds=(500)


for j in "${seeds[@]}"; do
        for i in "${spacing[@]}"; do
                echo "Run spacing $i with seed $j"

                echo "Generate... other... polymer molecule file."
                python3 generate_molecule.py --chain_length 30 --amine_spacing $i --filename data/polymer.molecule
                #Single run
                mpirun -np 13 --oversubscribe $LAMMPS -in in.amine_polymers_sh_short \
                        -var N $N \
                        -var P $P \
                        -var BOX_SIZE $BOX_SIZE \
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

# spacing=(3 4 5 7 10 15)
# seeds=(100 200 300 400 500)

# for j in "${seeds[@]}"; do
#         for i in "${spacing[@]}"; do
#                 echo "Run spacing $i with seed $j"

#                 echo "Generate... other... polymer molecule file."
#                 python3 generate_molecule.py --chain_length 30 --amine_spacing $i --filename data/polymer.molecule
#                 #Single run
#                 mpirun -np 13 --oversubscribe $LAMMPS -in in.amine_polymers_sh_short \
#                         -var N $N \
#                         -var P $P \
#                         -var BOX_SIZE $BOX_SIZE \
#                         -var patch_bond_r0 $patch_bond_r0 \
#                         -var patch_gauss_A_strength $patch_gauss_A_strength \
#                         -var patch_gauss_B_width $patch_gauss_B_width \
#                         -var output_traj $output_traj \
#                         -var output_traj_file $output_traj_file \
#                         -var output_msd $output_msd \
#                         -var output_msd_file $output_msd_file \
#                         -var output_amine $output_amine \
#                         -var output_amine_file $output_amine_file \
#                         -var output_press $output_press \
#                         -var output_press_file "$LOGDIR/press_${patch_gauss_A_strength}_spacing${i}_seed${j}.dat" \
#                         -var random_seed $j
#         done
# done

i=100

# # Single run
# python3 generate_molecule.py --chain_length 30 --amine_spacing 6 --filename data/polymer.molecule

# mpirun -np 12 --oversubscribe $LAMMPS -in in.amine_polymers_sh_short \
#         -var N $N \
#         -var P $P \
#         -var BOX_SIZE $BOX_SIZE \
#         -var patch_bond_r0 $patch_bond_r0 \
#         -var patch_gauss_A_strength $patch_gauss_A_strength \
#         -var patch_gauss_B_width $patch_gauss_B_width \
#         -var output_traj $output_traj \
#         -var output_traj_file $output_traj_file \
#         -var output_msd $output_msd \
#         -var output_msd_file $output_msd_file \
#         -var output_amine $output_amine \
#         -var output_amine_file $output_amine_file \
#         -var output_press $output_press \
#         -var output_press_file "$LOGDIR/press_test_freq_low_temp.dat" \
#         -var random_seed $i

# n_len_values=(2 4 8 16 64)

# for n_len in "${n_len_values[@]}"; do
#         mpirun -np 12 --oversubscribe $LAMMPS -in in.amine_polymers_sh_short \
#                 -var N $N \
#                 -var P $P \
#                 -var n_len $n_len \
#                 -var BOX_SIZE $BOX_SIZE \
#                 -var patch_bond_r0 $patch_bond_r0 \
#                 -var patch_gauss_A_strength $patch_gauss_A_strength \
#                 -var patch_gauss_B_width $patch_gauss_B_width \
#                 -var output_traj $output_traj \
#                 -var output_traj_file $output_traj_file \
#                 -var output_msd $output_msd \
#                 -var output_msd_file $output_msd_file \
#                 -var output_amine $output_amine \
#                 -var output_amine_file $output_amine_file \
#                 -var output_press 0 \
#                 -var output_press_file "$LOGDIR/n_len_${n_len}.dat" \
#                 -var random_seed $i
# done



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



# # Process each energy scale sequentially
# for scale in "${scales[@]}"; do
#   echo "========================================="
#   echo "Processing energy scale: $scale"
#   echo "========================================="
  
#   # Step 1: Run LAMMPS simulation
#   echo "Step 1/4: Running LAMMPS simulation for scale $scale"
#   mpirun -np 12 --oversubscribe \
#     "$LAMMPS" -in "$INPUT" \
#     -var patch_gauss_A_strength "$scale" \
#     -var output_file "$LOGDIR/energy_${scale}.lammpstrj"
  
#   echo "LAMMPS simulation completed successfully for scale $scale"
  
#   # Step 2: Find bonds
#   echo "Step 2/4: Finding bonds for scale $scale"
#   python3 find_reversible_bonds.py --traj "logs/energy_${scale}.lammpstrj" --output "data/amine_bonds_${scale}.csv"
  
#   echo "Bond finding completed successfully for scale $scale"
  
#   # Step 3: Find multiplets
#   echo "Step 3/4: Finding multiplets for scale $scale"
#   output_csv="data/amine_multiplets_${scale}.csv"
#   python3 find_multiplets.py --input "data/amine_bonds_${scale}.csv" --output "$output_csv"
  
#   echo "Multiplet finding completed successfully for scale $scale"

#   echo "Step 4/4: Calculating bond lifetimes for scale $scale"
#   python3 calc_bond_lifetime.py --csv_file "data/amine_bonds_${scale}.csv" --output "data/bond_lifetimes_${scale}.csv"
#   echo "Bond lifetime calculation completed successfully for scale $scale"

#   echo "All processing completed for energy scale $scale"
#   echo ""
# done

# echo "========================================="
# echo "All energy scales processed successfully!"
# echo "========================================="

