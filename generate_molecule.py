"""
Similar to generate_data.py but generates a .molecule file, which should be easier to use to generate random configs

Note that the formatting of individual atoms has been split across the sections for molecule entries
"""

import argparse


def generate_molecule_file(filename=None,
                           chain_length=25,
                           patch_spacing=5,
                           backbone_to_patch=0.35):

    if patch_spacing > chain_length:
        print("Warning: Patch spacing cannot be greater than chain length.")
    if patch_spacing == 0:
        patch_positions = []
    else:
        patch_positions = list(range(patch_spacing, chain_length, patch_spacing))
        if patch_positions[-1] == chain_length - 1:
            print("Warning: Last monomer is a patch. This affects chain length (slightly).")

    if filename is None:
        filename = f"data/polymer_chainlength_{chain_length}_patchspacing_{patch_spacing}_backbone2patch_{backbone_to_patch}.molecule"

    N = chain_length + len(patch_positions)  # total atoms

    with open(filename, 'w') as f:
        f.write("LAMMPS Description\n\n")

        # Counts 
        f.write(f"{N} atoms\n") # each patch position has 2 atoms (patch + dummy)
        f.write(f"{N - 1} bonds\n") # each patch position has 2 bonds (to dummy and dummy to patch)

        f.write("Coords\n\n")
        for i in range(chain_length):
            f.write(f"{i+1} {i+1}.0 0.0 0.0\n")
        for idx, pos in enumerate(patch_positions):
            patch_id = chain_length + idx + 1
            f.write(f"{patch_id} {pos + 1} {backbone_to_patch * 0.9} 0.0\n")
        f.write("\nTypes\n\n")
        for i in range(chain_length):
            f.write(f"{i+1} 1\n")
        for idx, pos in enumerate(patch_positions):
            patch_id = chain_length + idx + 1
            f.write(f"{patch_id} 2\n")

        f.write("\nBonds\n\n")
        for i in range(chain_length - 1):
            f.write(f"{i+1} 1 {i+1} {i+2}\n")
        for idx, pos in enumerate(patch_positions):
            patch_id = chain_length + idx + 1
            f.write(f"{chain_length + idx} 2 {pos + 1} {patch_id}\n")

    return filename


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--filename', type=str, default=None, help='Output filename')
    args.add_argument('--chain_length', type=int, default=25, help='Length of the polymer chain')
    args.add_argument('--patch_spacing', type=int, default=5, help='Spacing between sticky patches')
    args.add_argument('--backbone_to_patch', type=float, default=0.35, help='Distance from backbone to patch')
    args = args.parse_args()
    print(args)
    fname = generate_molecule_file(
        filename=args.filename,
        chain_length=args.chain_length,
        patch_spacing=args.patch_spacing,
        backbone_to_patch=args.backbone_to_patch
    )
    print(f"Generated molecule file: {fname}")