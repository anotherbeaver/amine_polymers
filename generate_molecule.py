chain_length = 25   # number of monomers in a chain
amine_spacing = 5   # spacing between amine groups
backbone_to_dummy = 1.0   # distance from backbone to dummy bead
dummy_to_patch = 0.45     # distance from dummy bead to patch bead

def generate_molecule_file(filename="data/polymer.molecule"):
    # sanity check
    if amine_spacing > chain_length:
        raise ValueError("Amine spacing cannot be greater than chain length.")

    # positions where amine groups are attached
    amine_positions = list(range(amine_spacing, chain_length, amine_spacing))

    # total atoms = backbone + 2 per amine group
    N_atoms = chain_length + len(amine_positions) * 2
    # total bonds = (chain_length - 1) backbone bonds + 2 per amine
    N_bonds = (chain_length - 1) + len(amine_positions)

    with open(filename, "w") as f:
        f.write("LAMMPS Description\n\n")
        f.write(f"{N_atoms} atoms\n")
        f.write(f"{N_bonds} bonds\n\n")

        # --- Coords section ---
        f.write("Coords\n\n")
        # backbone atoms along x-axis
        for i in range(chain_length):
            f.write(f"{i+1} {i*1.0:.3f} 0.0 0.0\n")
        # dummies + patches
        for idx, pos in enumerate(amine_positions):
            dummy_id = chain_length + idx * 2 + 1
            patch_id = chain_length + idx * 2 + 2
            f.write(f"{dummy_id} {(pos+1)*1.0:.3f} {backbone_to_dummy:.3f} 0.0\n")
            f.write(f"{patch_id} {(pos+1)*1.0:.3f} {backbone_to_dummy+dummy_to_patch:.3f} 0.0\n")

        # --- Types section ---
        f.write("\nTypes\n\n")
        # backbone atoms = type 1
        for i in range(chain_length):
            f.write(f"{i+1} 1\n")
        # dummy = type 2, patch = type 3
        for idx, pos in enumerate(amine_positions):
            dummy_id = chain_length + idx * 2 + 1
            patch_id = chain_length + idx * 2 + 2
            f.write(f"{dummy_id} 2\n")
            f.write(f"{patch_id} 3\n")

        # --- Molecules section ---
        f.write("\nMolecules\n\n")
        for i in range(N_atoms):
            f.write(f"{i+1} 1\n")

        # --- Bonds section ---
        f.write("\nBonds\n\n")
        bond_id = 1
        # backbone bonds
        for i in range(chain_length - 1):
            f.write(f"{bond_id} 1 {i+1} {i+2}\n")
            bond_id += 1
        # backbone–dummy + dummy–patch bonds
        for idx, pos in enumerate(amine_positions):
            dummy_id = chain_length + idx * 2 + 1
            patch_id = chain_length + idx * 2 + 2
            # backbone–dummy
            f.write(f"{bond_id} 2 {pos+1} {dummy_id}\n")
            bond_id += 1
            # dummy–patch
            # f.write(f"{bond_id} 3 {dummy_id} {patch_id}\n") # no bond type 3, rigid body
            bond_id += 1

if __name__ == "__main__":
    generate_molecule_file()
