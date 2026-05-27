#!/bin/bash
#SBATCH --job-name=plain_polymer
#SBATCH --account=st-jrottler-1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=20:00:00
#SBATCH --partition=skylake
#SBATCH --array=0-2
#SBATCH --output=logs/plain_polymer_%A_%a.out
#SBATCH --error=logs/plain_polymer_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=azhu13@student.ubc.ca

set -e

echo "Job ${SLURM_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}"
echo "Started on $(hostname) at $(date)"

# -----------------------------
# Environment
# -----------------------------
module purge
module load gcc openmpi python
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1

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
N_list=(30 100 200) # length of chains
N_beads_list=(24000 40000 80000)

# -----------------------------
# Map array index -> parameters
# -----------------------------

# Flatten 3D parameter space into single array index
index=$((SLURM_ARRAY_TASK_ID))
N=${N_list[$index]}
N_beads=${N_beads_list[$index]}

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

srun ${LAMMPS} -in "${INPUTEQ}" \
    -sf omp \
    -pk omp ${OMP_NUM_THREADS} \
    -var chain_length "${N}" \
    -var N_beads "${N_beads}"
srun ${LAMMPS} \
  -sf omp \
  -pk omp ${OMP_NUM_THREADS} \
  -in "${INPUTPROD}" \
  -var chain_length "${N}" \
  -var N_beads "${N_beads}"

echo "Finished at $(date)"
