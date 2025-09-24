#!/usr/bin/env python3

chain_length = 25   # number of monomers in a chain
amine_spacing = 5   # spacing between amine groups
backbone_to_donor = 1.0   # distance from backbone to donor bead
donor_to_h_dummy = 0.45     # distance from donor bead to h_dummy bead
backbone_to_acceptor = 0.45   # distance from backbone bead to acceptor bead

def generate_molecule_file(filename="data/polymer.molecule"):
    # sanity check
    if amine_spacing > chain_length:
        raise ValueError("Amine spacing cannot be greater than chain length.")

    # positions where amine groups are attached
    amine_positions = list(range(amine_spacing, chain_length, amine_spacing))

    # total atoms = backbone + 3 per amine group
    N_atoms = chain_length + len(amine_positions) * 3
    # total bonds = (chain_length - 1) backbone bonds + 3 per amine
    N_bonds = (chain_length - 1) + len(amine_positions) * 3

    with open(filename, "w") as f:
        f.write("LAMMPS Description\n\n")
        f.write(f"{N_atoms} atoms\n")
        f.write(f"{N_bonds} bonds\n\n")

        # --- Coords section ---
        f.write("Coords\n\n")
        # backbone atoms along x-axis
        for i in range(chain_length):
            f.write(f"{i+1} {i*1.0:.3f} 0.0 0.0\n")
        # donors + h_dummies + acceptors
        for idx, pos in enumerate(amine_positions):
            donor_id = chain_length + idx * 3 + 1
            h_dummy_id = chain_length + idx * 3 + 2
            acceptor_id = chain_length + idx * 3 + 3
            f.write(f"{donor_id} {(pos+1)*1.0:.3f} {backbone_to_donor:.3f} 0.0\n")
            f.write(f"{h_dummy_id} {(pos+1)*1.0:.3f} {backbone_to_donor+donor_to_h_dummy:.3f} 0.0\n")
            f.write(f"{acceptor_id} {(pos+1)*1.0:.3f} 0.0 {backbone_to_acceptor:.3f}\n") # init in z direction

        # --- Types section ---
        f.write("\nTypes\n\n")
        # backbone atoms = type 1
        for i in range(chain_length):
            f.write(f"{i+1} 1\n")
        # donor = type 2, h_dummy = type 3, acceptor = type 4
        for idx, pos in enumerate(amine_positions):
            donor_id = chain_length + idx * 3 + 1
            h_dummy_id = chain_length + idx * 3 + 2
            acceptor_id = chain_length + idx * 3 + 3
            f.write(f"{donor_id} 2\n")
            f.write(f"{h_dummy_id} 3\n")
            f.write(f"{acceptor_id} 4\n")

        # # --- Molecules section ---
        # f.write("\nMolecules\n\n")
        # for i in range(N_atoms):
        #     f.write(f"{i+1} 1\n")

        # --- Bonds section ---
        f.write("\nBonds\n\n")
        bond_id = 1
        # backbone bonds
        for i in range(chain_length - 1):
            f.write(f"{bond_id} 1 {i+1} {i+2}\n")
            bond_id += 1
        # backbone–donor + donor–h_dummy + backbone–acceptor bonds
        for idx, pos in enumerate(amine_positions):
            donor_id = chain_length + idx * 3 + 1
            h_dummy_id = chain_length + idx * 3 + 2
            acceptor_id = chain_length + idx * 3 + 3
            # backbone–donor
            f.write(f"{bond_id} 2 {pos+1} {donor_id}\n")
            bond_id += 1
            # donor–h_dummy
            f.write(f"{bond_id} 3 {donor_id} {h_dummy_id}\n")
            bond_id += 1
            # backbone–acceptor
            f.write(f"{bond_id} 4 {pos+1} {acceptor_id}\n")
            bond_id += 1

if __name__ == "__main__":
    generate_molecule_file()
