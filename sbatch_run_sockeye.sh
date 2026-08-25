#!/bin/bash
#SBATCH --job-name=sticky_polymer_sockeye
#SBATCH --account=st-jrottler-1
#SBATCH --partition=skylake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=0-8
#SBATCH --time=168:00:00
#SBATCH --output=logs/N_31_random_spacing_5_sticky_polymer_%A_%a.out
#SBATCH --error=logs/N_31_random_spacing_5_sticky_polymer_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=azhu13@student.ubc.ca

# this file is for the proper sticky production run

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
BINARY2TXT=../lammps_build/lammps/tools/binary2txt
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
patch_strength_list=(60 70 80 90 100 110 120 130 140) # strength of patch attraction (axis 3)
random_seed_list=(12345) # random seeds for Langevin thermostat (axis 4)

# -----------------------------
# Map array index -> parameters
# -----------------------------

# Flatten 3D parameter space into single array index
index=$((SLURM_ARRAY_TASK_ID))

# temporary, only have one axis rn
# TODO: need to update this part when we have more axes in the parameter space
N=${N_list[0]}
N_beads=${N_beads_list[0]}
patch_spacing=${patch_spacing_list[${index}]}
patch_strength=${patch_strength_list[0]}
random_seed=${random_seed_list[${index}]}

echo "Running N=${N}, N_beads=${N_beads}, patch_spacing=${patch_spacing}, patch_strength=${patch_strength}, random_seed=${random_seed}"

# record parameters to a shared lookup CSV (append-only)
timestamp=$(date +%Y-%m-%dT%H:%M:%S)
LOOKUP_FILE="${LOGDIR}/params_lookup.csv"
if [ ! -f "${LOOKUP_FILE}" ]; then
  printf "job_id,array_id,timestamp,N,N_beads,patch_spacing,patch_strength,random_seed,hostname,submit_dir,lammps_bin,input_file\n" > "${LOOKUP_FILE}"
fi
printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
  "${SLURM_JOB_ID:-}" "${SLURM_ARRAY_TASK_ID:-}" "${timestamp}" "${N}" "${N_beads}" "${patch_spacing}" "${patch_strength}" "${random_seed}" "$(hostname)" "${SLURM_SUBMIT_DIR:-$(pwd)}" "${LAMMPS}" "${INPUTPROD}" >> "${LOOKUP_FILE}"
echo "Appended parameters to ${LOOKUP_FILE}"

# -----------------------------
# Molecule generation
# -----------------------------
python3 generate_molecule.py \
    --chain_length ${N} \
    --patch_spacing ${patch_spacing} \
    --backbone_to_patch 0.25

# -----------------------------
# Run LAMMPS
# -----------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ${LAMMPS} -in "${INPUTEQ}" \
    -var chain_length "${N}" \
    -var N_beads "${N_beads}"\
    -var patch_spacing "${patch_spacing}" \
    -var patch_gauss_A_strength "${patch_strength}" \
    -var random_seed "${random_seed}"
srun ${LAMMPS} \
  -sf omp \
  -pk omp ${OMP_NUM_THREADS} \
  -in "${INPUTPROD}" \
  -var chain_length "${N}" \
  -var N_beads "${N_beads}" \
  -var patch_spacing "${patch_spacing}" \
  -var patch_gauss_A_strength "${patch_strength}" \
  -var random_seed "${random_seed}" 

# srun ${BINARY2TXT} "prodpatch_chainlength_31_patchspacing_${patch_spacing}_N_beads_24800_gaussA_${patch_strength}_r0_0.25_gaussB_10_seed_12345.bin.txt"

echo "Finished at $(date)"


# mpirun -np 6 ../lammps/build/lmp -in "in.polymers_eq" \
#     -var chain_length "30" \
#     -var N_beads "24000" \
#     -var patch_spacing "5" \
#     -var patch_gauss_A_strength "30" \
#     -var random_seed "12345"

# mpirun -np 6 ../lammps/build/lmp \
#   -in "in.polymers_prod" \
#   -var chain_length "30" \
#   -var N_beads "24000" \
#   -var patch_spacing "5" \
#   -var patch_gauss_A_strength "30" \
#   -var random_seed "12345" 
