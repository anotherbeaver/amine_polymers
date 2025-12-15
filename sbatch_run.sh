#!/bin/bash
#SBATCH --job-name=patchy_polymer_baseline
#SBATCH --output=logs/patchy_polymer_%A_%a.out
#SBATCH --error=logs/patchy_polymer_%A_%a.err
#SBATCH --nodes=1
#SBATCH --account=st-jrottler-1
#SBATCH --ntasks=16
#SBATCH --time=15:00:00         # max walltime hh:mm:ss
#SBATCH --output=logs/patchy_polymer_%A.out
#SBATCH --error=logs/patchy_polymer_%A.err
#SBATCH --mem=16G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=azhu13@student.ubc.ca

echo "Job started on $(hostname) at $(date)"

# =============================
# Environment
# =============================
module purge
module load gcc
module load openmpi
module load python

cd "$SLURM_SUBMIT_DIR"

export OMP_NUM_THREADS=1

# =============================
# Paths
# =============================
LAMMPS=../lammps_build/lammps/build/lmp
INPUT=in.amine_polymers_sh_arc
LOGDIR=logs
DATADIR=data

mkdir -p "$LOGDIR" "$DATADIR"

# =============================
# Simulation parameters
# =============================
patch_bond_r0=0.45
patch_gauss_A_strength=50
patch_gauss_B_width=13.8504255125

output_traj=0
output_traj_file="${DATADIR}/patch_trajectory.lammpstrj"

output_msd=0
output_msd_file="${LOGDIR}/patch_msd.dat"

output_amine=1
output_amine_file="${LOGDIR}/patch_locations.dat"

output_press=0

# =============================
# Molecule generation
# =============================
echo "Generating long polymer molecule..."
python3 generate_molecule.py \
    --chain_length 100 \
    --amine_spacing 1 \
    --filename "${DATADIR}/polymer_long.molecule"

echo "Generating short polymer molecule..."
python3 generate_molecule.py \
    --chain_length 28 \
    --amine_spacing 1 \
    --filename "${DATADIR}/polymer_short.molecule"

# =============================
# Parameter sweep
# =============================
spacing=(0 2)
seeds=(1000)

# =============================
# Main loop
# =============================
for seed in "${seeds[@]}"; do
    for sp in "${spacing[@]}"; do

        echo "--------------------------------------"
        echo "Spacing = ${sp}, Seed = ${seed}"
        echo "--------------------------------------"

        echo "Generating polymer molecule (spacing=${sp})..."
        python3 generate_molecule.py \
            --chain_length 30 \
            --amine_spacing "${sp}" \
            --filename "${DATADIR}/polymer.molecule"

        echo "Running LAMMPS with mpirun..."
        mpirun -np ${SLURM_NTASKS} --oversubscribe \
            ${LAMMPS} -in "${INPUT}" \
            -var patch_bond_r0 "${patch_bond_r0}" \
            -var patch_gauss_A_strength "${patch_gauss_A_strength}" \
            -var patch_gauss_B_width "${patch_gauss_B_width}" \
            -var output_traj "${output_traj}" \
            -var output_traj_file "${output_traj_file}" \
            -var output_msd "${output_msd}" \
            -var output_msd_file "${output_msd_file}" \
            -var output_amine "${output_amine}" \
            -var output_amine_file "${output_amine_file}" \
            -var output_press "${output_press}" \
            -var output_press_file "${LOGDIR}/press_A${patch_gauss_A_strength}_spacing${sp}_seed${seed}.dat" \
            -var random_seed "${seed}"

    done
done

echo "Job finished at $(date)"
