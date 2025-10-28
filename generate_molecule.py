"""
Similar to generate_data.py but generates a .molecule file, which should be easier to use to generate random configs

Note that the formatting of individual atoms has been split across the sections for molecule entries
"""




import argparse


def generate_molecule_file(filename='data/polymer.molecule',
                           chain_length=25,
                           amine_spacing=5,
                           backbone_to_patch=0.45):
    
    BOX_X = chain_length + 2 # length of monomer-monomer distances with extra padding (1)
    BOX_Y = 1.45 # patch height (.45) + padding (1)
    BOX_Z = 1.2 # padding (1.2)


    if amine_spacing > chain_length:
        raise ValueError("Amine spacing cannot be greater than chain length.")
    amine_positions = list(range(amine_spacing, chain_length, amine_spacing))
    if amine_positions[-1] == chain_length - 1:
        print("Warning: Last monomer is an amine patch. This affects chain length.")

    N = chain_length + len(amine_positions)  # total atoms
    
    with open(filename, 'w') as f:
        f.write("LAMMPS Description\n\n")

        # Counts 
        f.write(f"{N} atoms\n") # each amine position has 2 atoms (patch + dummy)
        f.write(f"{N - 1} bonds\n") # each amine position has 2 bonds (to dummy and dummy to patch)

        # # Masses
        # f.write("\nMasses\n\n")
        # f.write("1 1.0\n")    # backbone + dummy mass
        # f.write("2 1\n\n")   # patch mass TODO: figure out if we can reduce this without blowing up the simulation

        f.write("Coords\n\n")
        for i in range(chain_length):
            f.write(f"{i+1} {i+1}.0 0.0 0.0\n")
        for idx, pos in enumerate(amine_positions):
            patch_id = chain_length + idx + 1
            f.write(f"{patch_id} {pos + 1} {backbone_to_patch} 0.0\n")
        f.write("\nTypes\n\n")
        for i in range(chain_length):
            f.write(f"{i+1} 1\n")
        for idx, pos in enumerate(amine_positions):
            patch_id = chain_length + idx + 1
            f.write(f"{patch_id} 2\n")


        # f.write("\nMasses\n\n")
        # for i in range(N):
        #     f.write(f"{i+1} 1.0\n")

        f.write("\nBonds\n\n")
        for i in range(chain_length - 1):
            f.write(f"{i+1} 1 {i+1} {i+2}\n")
        for idx, pos in enumerate(amine_positions):
            patch_id = chain_length + idx + 1
            f.write(f"{chain_length - 1 + idx * 2} 2 {pos + 1} {patch_id}\n")

    return filename

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
    # generate_chain_positions()
    args = argparse.ArgumentParser()
    args.add_argument('--filename', type=str, default='data/polymer.molecule', help='Output filename')
    args.add_argument('--chain_length', type=int, default=25, help='Length of the polymer chain')
    args.add_argument('--amine_spacing', type=int, default=5, help='Spacing between amine patches')
    args.add_argument('--backbone_to_patch', type=float, default=0.45, help='Distance from backbone to patch')
    args = args.parse_args()
    print(args)
    fname = generate_molecule_file(
        filename=args.filename,
        chain_length=args.chain_length,
        amine_spacing=args.amine_spacing,
        backbone_to_patch=args.backbone_to_patch
    )
    print(f"Generated molecule file: {fname}")