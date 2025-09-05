"""
Generate data file for LAMMPS simulation

Output:
    data.particles: LAMMPS data file containing particle information
"""

def generate_data_file(filename='data.particles', num_chains=10, num_particles=10):
    with open(filename, 'w') as f:
        f.write("LAMMPS Description\n\n")

        # Counts 
        f.write(f"{num_chains} atoms\n")
        f.write(f"{num_chains - 1} bonds\n")
        # f.write(f"{num_particles - 2} angles\n")
        # f.write(f"0 dihedrals\n")
        # f.write(f"0 impropers\n\n")

        # Types
        f.write("1 atom types\n")
        f.write("1 bond types\n")
        # f.write("1 angle types\n")
        # f.write("0 dihedral types\n")
        # f.write("0 improper types\n\n")

        # 100x100x10 box
        f.write("-50.0 50.0 xlo xhi\n")
        f.write("-50.0 50.0 ylo yhi\n")
        f.write("-50.0 50.0 zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write("1 1.0\n\n")

        f.write("Atoms\n\n")
        for i in range(num_chains):
            f.write(f"{i+1} 1 1 {i*1.0} {0.0} {0.0}\n")

        f.write("\nBonds\n\n")
        for i in range(num_chains - 1):
            f.write(f"{i+1} 1 {i+1} {i+2}\n")
        
        # f.write("\nAngles\n\n")
        # for i in range(num_particles - 2):
        #     f.write(f"{i+1} 1 {i+1} {i+2} {i+3}\n")

if __name__ == "__main__":
    generate_data_file()