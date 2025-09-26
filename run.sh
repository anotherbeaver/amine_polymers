#!/bin/bash

LAMMPS=../lammps/build/lmp
INPUT=in.amine_polymers_sh
LOGDIR=logs

mkdir -p "$LOGDIR"

scales=(50.0 60.0 70.0 80.0 90.0 150.0)

core=0
for scale in "${scales[@]}"; do
  echo "Launching hbond_energy_scale=$scale on cores $core and $((core+1))"
  taskset -c $core,$((core+1)) mpirun -np 2 --bind-to core \
    "$LAMMPS" -in "$INPUT" \
    -var patch_gauss_A_strength "$scale" \
    -var output_file "$LOGDIR/energy_${scale}.dat" &
  core=$((core+2))
done

wait
echo "All jobs finished."
