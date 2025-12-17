# Original location: <PROJECT_ROOT>/generate_data.py
# Very old, originally used to create locations of boxes before
# molecule files were used
"""
Generate data file for LAMMPS simulation

The output data file will be a replicable module that can be used to create a larger system.

The general structure of the molecule is as follows:
    Backbone: 1 - 1 - 1 - 1 - 1 - 1 - 1 - 1
                      |
                  Dummy: 1
                      |
                  Patch: 2

    Where '1' is a backbone bead, '2' is an amine patch bead, and the dummy bead (also type '1') is used to create a rigid bond angle between the backbone and the patch.

Output:
    data.particles: LAMMPS data file containing single molecule information
"""

chain_length = 25  # number of monomers in a chain
amine_spacing = 10  # spacing between amine groups in a chain
BOX_X = chain_length + 2 # length of monomer-monomer distances with extra padding (1)
BOX_Y = 2.45 # height of one monomer-monomer distance (1) + patch height (.45) + padding (1)
BOX_Z = 1.2 # padding (1.2)


def generate_data_file(filename='data/data.particles'):
    if amine_spacing > chain_length:
        raise ValueError("Amine spacing cannot be greater than chain length.")
    amine_positions = list(range(0, chain_length, amine_spacing))
    if amine_positions[-1] == chain_length - 1:
        print("Warning: Last monomer is an amine patch. This affects chain length.")
    
    with open(filename, 'w') as f:
        f.write("LAMMPS Description\n\n")

        # Counts 
        f.write(f"{chain_length + len(amine_positions) * 2} atoms\n") # each amine position has 2 atoms (patch + dummy)
        f.write(f"{chain_length - 1 + len(amine_positions) * 2} bonds\n") # each amine position has 2 bonds (to dummy and dummy to patch)

        # Types
        f.write("2 atom types\n") # backbone and amine patch
        f.write("2 bond types\n") # backbone bond and patch-dummy bond

        # box is sized for a single chain with padding
        f.write(f"0.0 {BOX_X} xlo xhi\n")
        f.write(f"0.0 {BOX_Y} ylo yhi\n")
        f.write(f"0.0 {BOX_Z} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write("1 1.0\n")    # backbone + dummy mass
        f.write("2 1\n\n")   # patch mass TODO: figure out if we can reduce this without blowing up the simulation

        f.write("Atoms\n\n")
        for i in range(chain_length):
            f.write(f"{i+1} 1 1 {i+1} 0.0 0.0\n")
        for idx, pos in enumerate(amine_positions):
            dummy_id = chain_length + idx * 2 + 1
            patch_id = chain_length + idx * 2 + 2
            f.write(f"{dummy_id} 1 1 {pos + 1} 1.0 0.0\n")
            f.write(f"{patch_id} 2 2 {pos + 1} 1.45 0.0\n")

        f.write("\nBonds\n\n")
        for i in range(chain_length - 1):
            f.write(f"{i+1} 1 {i+1} {i+2}\n")
        for idx, pos in enumerate(amine_positions):
            dummy_id = chain_length + idx * 2 + 1
            patch_id = chain_length + idx * 2 + 2
            f.write(f"{chain_length - 1 + idx * 2 + 1} 1 {pos + 1} {dummy_id}\n")
            f.write(f"{chain_length - 1 + idx * 2 + 2} 2 {dummy_id} {patch_id}\n")

def generate_chain_positions():
    """
    Generate positions for chains of particles.
    
    Returns:
        List of tuples: (particle_id, x, y, z)
    """
    # positions = []
    # for chain_id in range(num_chains):
    #     for i in range(chain_length):
    #         x = chain_id * amine_spacing
    #         y = 0.0
    #         z = i * amine_spacing
    #         positions.append((chain_id * chain_length + i + 1, x, y, z))
    # return positions
    pass

if __name__ == "__main__":
    generate_chain_positions()
    generate_data_file()