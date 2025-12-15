#!/bin/bash
#SBATCH --job-name=patchy_polymer
#SBATCH --nodes=1
#SBATCH --output=logs/patchy_polymer_%A_%a.out
#SBATCH --error=logs/patchy_polymer_%A_%a.err
#SBATCH --ntasks=10
#SBATCH --time=02:00:00
#SBATCH --partition=skylake

# =============================
# Paths and LAMMPS executable
# =============================
LAMMPS=../lammps_build/lammps/build/lmp # Path to LAMMPS executable on HPC
INPUT_SHORT=in.amine_polymers_sh_short
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
output_traj_file="$DATADIR/patch_trajectory.lammpstrj"
output_msd=0
output_msd_file="$LOGDIR/patch_msd.dat"
output_amine=1
output_amine_file="$LOGDIR/patch_locations.dat"
output_press=0

# =============================
# Generate molecules
# =============================
echo "Generate long polymer molecule file."
python3 generate_molecule.py --chain_length 100 --amine_spacing 1 \
    --filename "$DATADIR/polymer_long.molecule"

echo "Generate short polymer molecule file."
python3 generate_molecule.py --chain_length 28 --amine_spacing 1 \
    --filename "$DATADIR/polymer_short.molecule"

# =============================
# Define spacings and seeds
# =============================
spacing=(15)
seeds=(100)

# =============================
# Main loop: generate polymer and run LAMMPS
# =============================
for j in "${seeds[@]}"; do
    for i in "${spacing[@]}"; do
        echo "Run spacing $i with seed $j"

        echo "Generate polymer molecule file for spacing $i."
        python3 generate_molecule.py --chain_length 30 --amine_spacing "$i" \
            --filename "$DATADIR/polymer.molecule"

        # Run LAMMPS using srun (Slurm-managed parallel execution)
        srun $LAMMPS -in "$INPUT_SHORT" \
            -var patch_bond_r0 "$patch_bond_r0" \
            -var patch_gauss_A_strength "$patch_gauss_A_strength" \
            -var patch_gauss_B_width "$patch_gauss_B_width" \
            -var output_traj "$output_traj" \
            -var output_traj_file "$output_traj_file" \
            -var output_msd "$output_msd" \
            -var output_msd_file "$output_msd_file" \
            -var output_amine "$output_amine" \
            -var output_amine_file "$output_amine_file" \
            -var output_press "$output_press" \
            -var output_press_file "$LOGDIR/press_${patch_gauss_A_strength}_spacing${i}_seed${j}.dat" \
            -var random_seed "$j"
    done
done
