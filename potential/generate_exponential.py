"""
Generate an exponential potential file for LAMMPS.

This potential applies ONLY to the amine groups forming transient bonds, superimposed on
the WCA potential.

Form of the potential:
    V(r) = A * exp(-r/(2 * sigma_G^2))

Based on Rafaelli and Ellenbroek, Soft Matter, 2021 10.1039/d1sm00091h

NOTE: THIS IS NOT USED IN THE FINAL SIMULATIONS. INSTEAD, A GAUSSIAN POTENTIAL IS USED.
There was an error in the paper, it should be a Gaussian potential, not an exponential.
"""

import numpy as np

# Parameters for the exponential potential
A = 50 # Tune this parameter by replicating the Arrhenius equation (see paper)
sigma_G = 0.19 # from paper

def exponential_potential(r):
    return - A * np.exp(-r / (2 * sigma_G**2))

def exponential_force(r):
    return -(A / (2 * sigma_G**2)) * np.exp(-r / (2 * sigma_G**2))

def generate_exponential_file(filename='potential/exp.table', r_min=0.00001, r_max=2.0, num_points=500):
    r_values = np.linspace(r_min, r_max, num_points)
    V_values = exponential_potential(r_values)
    F_values = exponential_force(r_values)

    with open(filename, 'w') as f:
        f.write("EXP\n")
        f.write(f"N {num_points} R {r_min} {r_max}\n\n")
        # f.write("hi")
        # f.write("# r (distance)    V(r) (potential energy)\n\n")
        for i, (r, V, F) in enumerate(zip(r_values, V_values, F_values)):
            f.write(f"{i} {r:.15f} {V:.15f} {F:.15f}\n")


if __name__ == "__main__":
    generate_exponential_file()
