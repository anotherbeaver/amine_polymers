import matplotlib.pyplot as plt

def plot_exp(filename):
    r = []
    V = []
    F = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[5:]:  # Skip header lines
            parts = line.split()
            r.append(float(parts[1]))
            V.append(float(parts[2]))
            F.append(float(parts[3]))

    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(r, V, label='V(r)')
    plt.xlabel('r')
    plt.ylabel('Potential Energy V(r)')
    plt.title('Exponential Potential V(r)')
    plt.grid()
    plt.legend()
    plt.subplot(1, 2, 2)
    plt.plot(r, F, label='F(r)', color='orange')
    plt.xlabel('r')
    plt.ylabel('Force F(r)')
    plt.title('Force F(r)')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_exp('potential/exp.table')    