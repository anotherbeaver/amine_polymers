#!/bin/bash
#SBATCH --job-name=patchy_polymer
#SBATCH --account=st-jrottler-1
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=25:00:00
#SBATCH --mem=16G
#SBATCH --partition=skylake
#SBATCH --array=0-5              # total combinations = spacing * seeds * A_strengths
#SBATCH --output=logs/patchy_polymer_%A_%a.out
#SBATCH --error=logs/patchy_polymer_%A_%a.err
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
INPUT=in.amine_polymers_sh_arc
LOGDIR=logs
DATADIR=data
mkdir -p "$LOGDIR" "$DATADIR"

# -----------------------------
# Fixed parameters
# -----------------------------
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
output_amine=0
output_amine_file="$LOGDIR/patch_locations.dat"
output_press=0
# output_press_file="$LOGDIR/press_${patch_gauss_A_strength}_${P}.dat"


# -----------------------------
# Parameter space
# -----------------------------
spacing=(4)
seeds=(1000)
A_strengths=(50 75 100)   # new parameter

# -----------------------------
# Map array index → parameters
# -----------------------------
n_spacing=${#spacing[@]}
n_seeds=${#seeds[@]}
n_A=${#A_strengths[@]}

# Flatten 3D parameter space into single array index
tmp=$((SLURM_ARRAY_TASK_ID))
seed_index=$(( tmp / (n_spacing*n_A) ))
tmp=$(( tmp % (n_spacing*n_A) ))
spacing_index=$(( tmp / n_A ))
A_index=$(( tmp % n_A ))

seed=${seeds[$seed_index]}
sp=${spacing[$spacing_index]}
patch_gauss_A_strength=${A_strengths[$A_index]}

echo "Running spacing=${sp}, seed=${seed}, A_strength=${patch_gauss_A_strength}"

# -----------------------------
# Molecule generation
# -----------------------------
python3 generate_molecule.py \
    --chain_length 30 \
    --amine_spacing "${sp}" \
    --filename "${DATADIR}/polymer_${seed}_${sp}_${patch_gauss_A_strength}.molecule"

# -----------------------------
# Run LAMMPS
# -----------------------------
srun ${LAMMPS} -in "${INPUT}" \
    -var patch_bond_r0 "${patch_bond_r0}" \
    -var patch_gauss_A_strength "${patch_gauss_A_strength}" \
    -var patch_gauss_B_width "${patch_gauss_B_width}" \
    -var output_traj "${output_traj}" \
    -var output_msd "${output_msd}" \
    -var output_amine "${output_amine}" \
    -var output_press "${output_press}" \
    -var output_press_file "${LOGDIR}/press_${patch_gauss_A_strength}_spacing${sp}_seed${seed}.dat" \
    -var molecule_file "${DATADIR}/polymer_${seed}_${sp}_${patch_gauss_A_strength}.molecule" \
    -var random_seed "${seed}"

echo "Finished at $(date)"
