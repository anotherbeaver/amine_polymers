#!/bin/bash
#SBATCH --job-name=plain_polymer
#SBATCH --account=st-jrottler-1
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --partition=skylake
#SBATCH --output=logs/test_sticky_polymer_%A_%a.out
#SBATCH --error=logs/test_sticky_polymer_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=azhu13@student.ubc.ca

# this file is for benchmarking purposes

set -e

echo "Job ${SLURM_JOB_ID}"
echo "Started on $(hostname) at $(date)"

# -----------------------------
# Environment
# -----------------------------
module purge
module load gcc openmpi python
cd "$SLURM_SUBMIT_DIR"

# -----------------------------
# Paths
# -----------------------------
LAMMPS=../lammps_build/lammps/build/lmp
INPUTEQ=in.polymers_eq
INPUTPROD=in.polymers_prod
LOGDIR=logs
DATADIR=data
DUMPSDIR=data/dumps
AUTOCORRDIR=data/autocorrelation
PRESSDIR=data/pressures
RESTARTSDIR=data/restarts

mkdir -p "$LOGDIR" "$DATADIR" "$DUMPSDIR" "$AUTOCORRDIR" "$PRESSDIR" "$RESTARTSDIR"

# -----------------------------
# Parameter space
# -----------------------------
# N_list=(30 100 200) # lengths of chains (axis 1)
# N_beads_list=(24000 40000 80000) # number of beads (total system size), each corresponds to a chain length (axis 1)
# patch_spacing=(15 10 5 4 3 2 1) # spacing between patches (axis 2)
N_list=(30) # lengths of chains (axis 1)
N_beads_list=(24000) # number of beads (total system size), each corresponds to a chain length (axis 1)
patch_spacing_list=(5) # spacing between patches (axis 2)
patch_strength_list=(0 20 50 75 100) # strength of patch attraction (axis 3)
random_seed_list=(12345) # random seeds for Langevin thermostat (axis 4)

# -----------------------------
# Map array index -> parameters
# -----------------------------

Flatten 3D parameter space into single array index
index=$((SLURM_ARRAY_TASK_ID))

# temporary, only have one axis rn
# TODO: need to update this part when we have more axes in the parameter space
N=${N_list[index]}
N_beads=${N_beads_list[index]}
patch_spacing=${patch_spacing_list[index]}
patch_strength=${patch_strength_list[index]}
random_seed=${random_seed_list[index]}

echo "Running N=${N}, N_beads=${N_beads}, patch_spacing=${patch_spacing}, patch_strength=${patch_strength}, random_seed=${random_seed}"

# -----------------------------
# Molecule generation
# -----------------------------
python3 generate_molecule.py \
    --chain_length ${N} \
    --patch_spacing ${patch_spacing} \

# -----------------------------
# Run LAMMPS
# -----------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# srun ${LAMMPS} -in "${INPUTEQ}" \
#     -var chain_length "${N}" \
#     -var N_beads "${N_beads}"
srun ${LAMMPS} \
  -sf omp \
  -pk omp ${OMP_NUM_THREADS} \
  -in "${INPUTPROD}" \
  -var chain_length "${N}" \
  -var N_beads "${N_beads}"

echo "Finished at $(date)"
