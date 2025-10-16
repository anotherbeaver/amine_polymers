#!/bin/bash


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
# scales=(35 40 45 50 55 60 75 80)
scales=(50)



# lengths=(0.2 0.3 0.4 0.5 0.6 0.7)
# lengths=(0.2 0.4 0.6)
# lengths=(0.2 0.3)
lengths=(0.45)

# n_patch=(1 2 3 4 5)
n_patch=(5)


for n in "${n_patch[@]}"; do
  echo "Generating molecule with $n amine patches."
  python3 generate_molecule_no_dummy.py --amine_spacing $n
  for scale in "${scales[@]}"; do
    for length in "${lengths[@]}"; do
      echo "Launching hbond_energy_scale=$scale and hbond_distance_cutoff=$length"
      mpirun -np 12 --oversubscribe "$LAMMPS" -in "$INPUT" \
        -var patch_gauss_A_strength "$scale" \
        -var patch_bond_r0 "$length" \
        -var output_file "$LOGDIR/patch_n_${n}_${scale}_${length}.dat"
    done
  done
done


wait
echo "All LAMMPS jobs finished."



# # Process each energy scale sequentially
# for scale in "${scales[@]}"; do
#   echo "========================================="
#   echo "Processing energy scale: $scale"
#   echo "========================================="
  
#   # # Step 1: Run LAMMPS simulation
#   # echo "Step 1/4: Running LAMMPS simulation for scale $scale"
#   # mpirun -np 12 --oversubscribe \
#   #   "$LAMMPS" -in "$INPUT" \
#   #   -var patch_gauss_A_strength "$scale" \
#   #   -var output_file "$LOGDIR/energy_${scale}.lammpstrj"
  
#   # echo "LAMMPS simulation completed successfully for scale $scale"
  
#   # # Step 2: Find bonds
#   # echo "Step 2/4: Finding bonds for scale $scale"
#   # python3 find_reversible_bonds.py --traj "logs/energy_${scale}.lammpstrj" --output "data/amine_bonds_${scale}.csv"
  
#   # echo "Bond finding completed successfully for scale $scale"
  
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

