"""
Generate data file for LAMMPS simulation

Output:
    data.particles: LAMMPS data file containing particle information
"""

BOX_X, BOX_Y, BOX_Z = 100.0, 100.0, 100.0
chain_length = 10  # number of particles in a chain
amine_spacing = 3  # spacing between amine groups in a chain
num_chains = 10
num_particles = 10

atoms = [] # atoms list (id, x, y, z)

def generate_data_file(filename='data.particles'):
    with open(filename, 'w') as f:
        f.write("LAMMPS Description\n\n")

        # Counts 
        f.write(f"{num_chains} atoms\n")
        f.write(f"{num_chains - 1} bonds\n")

        # Types
        f.write("1 atom types\n")
        f.write("1 bond types\n")

        # 100x100x100 box
        f.write(f"-{BOX_X/2} {BOX_X/2} xlo xhi\n")
        f.write(f"-{BOX_Y/2} {BOX_Y/2} ylo yhi\n")
        f.write(f"-{BOX_Z/2} {BOX_Z/2} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write("* 1.0\n\n") # all masses are 1.0

        f.write("Atoms\n\n")

        for i in range(num_chains):
            f.write(f"{i+1} 1 1 {i*1.0} {0.0} {0.0}\n")
        f.write(f" 1 1 {i*1.0} {0.0} {0.0}\n")
        

        f.write("\nBonds\n\n")
        for i in range(num_chains - 1):
            f.write(f"{i+1} 1 {i+1} {i+2}\n")

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