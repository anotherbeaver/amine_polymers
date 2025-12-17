# Original location: <PROJECT_ROOT>/run_patchy.py
# Script to try running purely with python with MDMC approach
import numpy as np
import random
import os
from lammps import lammps

# --------------------------
# User parameters
# --------------------------
LAMMPS_EXEC = "../lammps/build/lmp"
N = 200                    # number of molecules
BOX_SIZE = 40              # simulation box
PATCH_BOND_R0 = 0.45       # equilibrium distance for patch bonds
PATCH_GAUSS_A_STRENGTH = 75
PATCH_GAUSS_B_WIDTH = 2.63
STEPS_PER_UPDATE = 10       # tau: steps before updating reversible bonds
TOTAL_STEPS = 200
CUT_OFF = 1.5               # LJ distance cutoff for patch bonds
BOND_TYPE = 2               # harmonic patch bond type
BREAK_PROB = 0.1            # probability to break a bond
BONDS_FILE = "data/current_bonds.npy"  # optional persistence

# --------------------------
# Initialize bonded pairs
# --------------------------
if os.path.exists(BONDS_FILE):
    bonded_pairs = set(map(tuple, np.load(BONDS_FILE)))
else:
    bonded_pairs = set()

# --------------------------
# Initialize LAMMPS
# --------------------------
lmp = lammps()

# --------------------------
# LAMMPS commands
# --------------------------
lmp.command("units lj")
lmp.command("dimension 3")
lmp.command("boundary p p p")
lmp.command("atom_style bond")
lmp.command("log logs/log.amine_polymers")

# create box
lmp.command(f"region box block 0 {BOX_SIZE} 0 {BOX_SIZE} 0 {BOX_SIZE}")
lmp.command("create_box 2 box bond/types 2 extra/bond/per/atom 4 extra/special/per/atom 20")

# masses
lmp.command("mass 1 1.0")
lmp.command("mass 2 1.0")

# molecule
lmp.command("molecule polymer data/polymer.molecule")

# bonds
lmp.command("bond_style harmonic")
lmp.command("bond_coeff 1 50.0 1.0")
lmp.command(f"bond_coeff 2 1000 {PATCH_BOND_R0}")

# non-bonded
lmp.command(f"pair_style hybrid/overlay lj/cut 1.4592006628 gauss 2.0")
lmp.command("pair_coeff 1 1 lj/cut 1.0 1.3")
lmp.command("pair_coeff 1 2 lj/cut 1.0 1.3")
lmp.command("pair_coeff 2 2 lj/cut 0.0 1.3")
lmp.command("pair_modify shift yes")
lmp.command(f"pair_coeff 2 2 gauss {PATCH_GAUSS_A_STRENGTH} {PATCH_GAUSS_B_WIDTH} 2.0")

# create atoms
lmp.command(f"create_atoms 0 random {N} 12345 box mol polymer 12345")

# groups
lmp.command("group monomer type 1")
lmp.command("group amine type 2")

# neighbor list
lmp.command("neighbor 0.3 bin")
lmp.command("neigh_modify every 1 delay 0 check yes")

# energy minimization
lmp.command("min_style cg")
lmp.command("minimize 1.0e-6 1.0e-8 10000 100000")
lmp.command("reset_timestep 0")

# velocities and integrator
lmp.command("timestep 0.005")
lmp.command("velocity all create 1.0 1")
lmp.command(f"fix 1 all nvt temp 1.0 1.0 {0.5*0.005}")

# dumps
lmp.command("dump 2 all atom 10 data/dump.lammpstrj")
lmp.command("dump aminepos amine custom 10 data/amine_positions.lammpstrj id type x y z")


# Bond update function
def update_bonds(lmp):
    global bonded_pairs
    types = lmp.numpy.extract_atom("type")
    x = lmp.numpy.extract_atom("x")
    amines = np.where(types == 2)[0]

    # create new candidate bonds
    new_bonds = []
    for i, idx_i in enumerate(amines):
        for idx_j in amines[i+1:]:
            pair = (min(idx_i+1, idx_j+1), max(idx_i+1, idx_j+1))
            if pair in bonded_pairs:
                continue
            dist = np.linalg.norm(x[idx_i] - x[idx_j])
            if dist <= CUT_OFF:
                new_bonds.append(pair)

    # decide which existing bonds to keep
    bonds_to_keep = [pair for pair in bonded_pairs if random.random() >= BREAK_PROB]

    # update bonded_pairs
    bonded_pairs = set(bonds_to_keep + new_bonds)

    # delete all type-2 bonds
    if bonded_pairs:
        lmp.command(f"delete_bonds all bond {BOND_TYPE} remove special")

    # recreate bonds
    for i, j in bonded_pairs:
        lmp.command(f"create_bonds many {BOND_TYPE} {i} {j}")

    print(f"[Python] created {len(new_bonds)} bonds, kept {len(bonds_to_keep)} existing bonds")


# Main simulation loop
steps_remaining = TOTAL_STEPS
while steps_remaining > 0:
    steps = min(STEPS_PER_UPDATE, steps_remaining)
    lmp.command(f"run {steps}")
    update_bonds(lmp)
    steps_remaining -= steps

# Save bonded pairs
os.makedirs(os.path.dirname(BONDS_FILE), exist_ok=True)
np.save(BONDS_FILE, np.array(list(bonded_pairs)))

print("Simulation finished. Bonded pairs saved.")
