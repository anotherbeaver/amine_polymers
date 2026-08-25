#!/bin/bash

LAMMPS=lmp
BINARY2TXT=binary2txt
INPUTEQ=in.polymers_eq
INPUTPROD=in.polymers_prod
LOGDIR=logs
DATADIR=data
DUMPSDIR=data/dumps
AUTOCORRDIR=data/autocorrelation
PRESSDIR=data/pressures
RESTARTSDIR=data/restarts
PATCHESDIR=data/patches

mkdir -p "$LOGDIR" "$DATADIR" "$DUMPSDIR" "$AUTOCORRDIR" "$PRESSDIR" "$RESTARTSDIR" "$PATCHESDIR"

# -----------------------------
# Parameter space
# -----------------------------
N_list=(31) # lengths of chains (axis 1)
N_beads_list=(24800) # number of beads (total system size), each corresponds to a chain length (axis 1)
patch_spacing_list=(5) # spacing between patches (axis 2)
patch_strength_list=(0) # strength of patch attraction (axis 3)
random_seed_list=(12345) # random seeds for Langevin thermostat (axis 4)

# -----------------------------
# Map array index -> parameters
# -----------------------------

# Flatten 3D parameter space into single array index
# index=$((SLURM_ARRAY_TASK_ID)) # this is for slurm, we are just running this script

# temporary, only have one axis rn
for index in "${!patch_strength_list[@]}"; do
  N=${N_list[0]}
  N_beads=${N_beads_list[0]}
  patch_spacing=${patch_spacing_list[0]}
  patch_strength=${patch_strength_list[index]}
  random_seed=${random_seed_list[0]}

  echo "Running N=${N}, N_beads=${N_beads}, patch_spacing=${patch_spacing}, patch_strength=${patch_strength}, random_seed=${random_seed}"

  # -----------------------------
  # Molecule generation
  # -----------------------------
  # python3 generate_molecule.py \
  #     --chain_length ${N} \
  #     --patch_spacing ${patch_spacing} \
  #     --backbone_to_patch 0.25 \
  #     --random 1

  # -----------------------------
  # Run LAMMPS
  # -----------------------------
  export OMP_NUM_THREADS=2
  export OMP_PROC_BIND=close

  mpirun -n 8 ${LAMMPS} -in "${INPUTEQ}" \
      -var chain_length "${N}" \
      -var N_beads "${N_beads}"\
      -var patch_spacing "${patch_spacing}" \
      -var patch_gauss_A_strength "${patch_strength}" \
      -var random_seed "${random_seed}"
      # -var random_patch_spacing 1

  # mpirun -n 8 ${LAMMPS} \
  #   -sf omp \
  #   -pk omp ${OMP_NUM_THREADS} \
  #   -in "${INPUTPROD}" \
  #   -var chain_length "${N}" \
  #   -var N_beads "${N_beads}" \
  #   -var patch_spacing "${patch_spacing}" \
  #   -var patch_gauss_A_strength "${patch_strength}" \
  #   -var random_seed "${random_seed}" 

  echo "Finished at $(date)"

done