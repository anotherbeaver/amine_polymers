#!/bin/bash
#SBATCH --job-name=plain_polymer
#SBATCH --account=st-jrottler-1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --partition=skylake
#SBATCH --output=logs/test_plain_polymer_%A_%a.out
#SBATCH --error=logs/test_plain_polymer_%A_%a.err
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
# N_list=(30 100 200) # length of chains
# N_beads_list=(24000 40000 80000)

# -----------------------------
# Map array index -> parameters
# -----------------------------

# Flatten 3D parameter space into single array index
# index=$((SLURM_ARRAY_TASK_ID))
N=29
N_beads=20000

echo "Running N=${N}, N_beads=${N_beads}"

# -----------------------------
# Molecule generation
# -----------------------------
python3 generate_molecule.py \
    --chain_length ${N} \
    --amine_spacing 0 

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
