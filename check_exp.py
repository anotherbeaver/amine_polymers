import matplotlib.pyplot as plt

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

    if len(r1) != len(r2):
        raise ValueError("Potential tables have different r-grid lengths!")

    e_total = [e1[i] + e2[i] for i in range(len(r1))]

    # plt.plot(r2, e1, label='Backbone-Patch')
    plt.plot(r2, e2, label='Patch-Patch')
    # plt.plot(r2, e_total, label='Sum', linestyle='--')
    plt.xlabel('r')
    plt.ylabel('U(r)')
    # plt.ylim(-5, 5)
    plt.legend()
    plt.title('Summed Interaction Potentials')
    # plt.savefig('data/summed_potentials.png', dpi=300)
    plt.show()
    # plt.close()

if __name__ == "__main__":
    plot_summed_potentials()

