#!/bin/bash

LAMMPS=../lammps/build/lmp
INPUT=in.amine_polymers
LOGDIR=logs

mkdir -p "$LOGDIR"

scales=(0.0 3.0 5.0 10.0 20.0 50.0 100.0)

core=0
for scale in "${scales[@]}"; do
  echo "Launching hbond_energy_scale=$scale on cores $core and $((core+1))"
  taskset -c $core,$((core+1)) mpirun -np 2 --map-by core --bind-to core \
    "$LAMMPS" -in "$INPUT" \
    -var hbond_energy_scale "$scale" \
    -var output_file "$LOGDIR/energy_${scale}.dat" &
  core=$((core+2))
done

wait
echo "All jobs finished."
