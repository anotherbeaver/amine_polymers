#!/bin/bash

LAMMPS=../lammps/build/lmp
INPUT=in.amine_polymers_sh
LOGDIR=data

mkdir -p "$LOGDIR"

# scales=(50.0 60.0 70.0 80.0 90.0 150.0)
# scales=(50.0 60.0 70.0)
scales=(49.0 99.0 149.0 200.0 300.0 400.0)
# lengths=(0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7)
# lengths=(0.6 0.65 0.7)




# core=0
# for scale in "${scales[@]}"; do
#   echo "Launching hbond_energy_scale=$scale on cores $core and $((core+1))"
#   taskset -c $core,$((core+1)) mpirun -np 2 --bind-to core \
#     "$LAMMPS" -in "$INPUT" \
#     -var patch_gauss_A_strength "$scale" \
#     -var output_file "$LOGDIR/energy_${scale}.lammpstrj" &
#   core=$((core+2))
# done

# wait
# echo "All LAMMPS jobs finished."

# for scale in "${scales[@]}"; do
#   # find the bonds
#   python3 find_reversible_bonds.py --traj "data/energy_${scale}.lammpstrj" --output "data/amine_bonds_${scale}.csv"
# done

wait
echo "All bond finding jobs finished."

for scale in "${scales[@]}"; do
  # find multiple bonds
  output_csv="data/amine_multiplets_${scale}.csv"
  python3 find_multiplets.py --input "data/amine_bonds_${scale}.csv" --output "$output_csv"
  echo "Job with hbond_energy_scale=$scale completed."
done

wait
echo "All jobs finished."

