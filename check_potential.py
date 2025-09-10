import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


def shift_and_sum(r1, e1, r2, e2, x0):
    # interpolate U2
    f2 = interp1d(r2, e2, kind='linear', bounds_error=False, fill_value=0.0)
    
    # evaluate shifted U2 at r1 positions
    e2_shifted = f2(r1 - x0)  # shift by x0
    
    # sum potentials
    e_total = e1 + e2_shifted
    return e2_shifted, e_total

def plot_summed_potentials():
    pot_1 = "data/pair_coeff_1_2.txt"
    pot_2 = "data/pair_coeff_2_2.txt"

    r1, e1 = [], []
    r2, e2 = [], []

    def read_pot(filename):
        r_vals, e_vals = [], []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip() == "" or line.startswith("#") or line[0].isalpha():
                    continue  # skip comments, headers
                try:
                    idx, r, e, fr = map(float, line.split())
                except ValueError:
                    continue  # skip bad lines
                r_vals.append(r)
                e_vals.append(e)
        return r_vals, e_vals

    r1, e1 = read_pot(pot_1)
    r2, e2 = read_pot(pot_2)

    r1, e1 = np.array(r1), np.array(e1)
    r2, e2 = np.array(r2), np.array(e2)

    if len(r1) != len(r2):
        raise ValueError("Potential tables have different r-grid lengths!")
    
    x0 = 0.9
    e2_shifted, e_total = shift_and_sum(r1, e1, r2, e2, x0)
    # e_total = e1 + e2

    plt.plot(r1, e1, label='Backbone-Patch (WCA)')
    plt.plot(r1, e2, label='Patch-Patch (Gaussian)')
    plt.plot(r1, e_total, label='Sum', linestyle='--')
    plt.xlabel('r')
    plt.ylabel('U(r)')
    plt.ylim(-40, 40)
    plt.legend()
    plt.title('Summed Interaction Potentials')
    # plt.savefig('data/summed_potentials.png', dpi=300)
    plt.show()
    # plt.close()

if __name__ == "__main__":
    plot_summed_potentials()

