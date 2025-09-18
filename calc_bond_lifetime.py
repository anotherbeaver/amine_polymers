import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def compute_bond_lifetimes(data, timestep_size=10):
    """
    Compute bond lifetimes from bond formation and breaking data.

    Parameters:
        data (pd.DataFrame): DataFrame containing columns:
            ['timestep', 'id1', 'id2', 'x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'distance']

    Returns:
        dict: Keys are (id1, id2) tuples, values are lists of bond lifetimes (in timesteps)
        dict: Keys are (id1, id2) tuples, values are lists of lists of timesteps when each bond was active
    """
    bond_lifetimes = {}
    bond_timesteps = {}

    for (id1, id2), group in data.groupby(["id1", "id2"]):
        bond_key = (id1, id2)
        timesteps = group['timestep'].values
        norm_timesteps = timesteps // timestep_size  # normalize to start at 0 and convert to index

        # Identify continuous segments of timesteps
        segments = np.split(norm_timesteps, np.where(np.diff(norm_timesteps) != 1)[0] + 1)

        lifetimes = [len(segment) for segment in segments if len(segment) > 0] * timestep_size
        bond_lifetimes[bond_key] = lifetimes
        bond_timesteps[bond_key] = [segment.tolist() for segment in segments if len(segment) > 0]

    return bond_lifetimes, bond_timesteps

def plot_bond_lifetimes(bond_lifetimes, title='Bond Lifetimes Histogram'):
    """
    Plot histogram of bond lifetimes.
    Parameters:
        bond_lifetimes (dict): Keys are (id1, id2) tuples, values are lists of bond lifetimes (in timesteps)
    """

    all_lifetimes = [lifetime for lifetimes in bond_lifetimes.values() for lifetime in lifetimes]

    plt.hist(all_lifetimes, bins=100, density=True)
    plt.xlabel('Bond Lifetime (timesteps)')
    plt.ylabel('Probability Density')
    plt.title(title)
    plt.show()

def save_bond_lifetimes(bond_lifetimes, filename='data/bond_lifetimes.csv'):
    """
    Save bond lifetimes to a CSV file.
    
    Parameters:
        bond_lifetimes (dict): Keys are (id1, id2) tuples, values are lists of bond lifetimes (in timesteps)
        filename (str): Output CSV file path
    """
    rows = []
    for (id1, id2), lifetimes in bond_lifetimes.items():
        for lifetime in lifetimes:
            rows.append([id1, id2, lifetime])
    df = pd.DataFrame(rows, columns=['id1', 'id2', 'lifetime'])
    df.to_csv(filename, index=False)

if __name__ == "__main__":
    # Specify the path to your CSV file
    csv_file_path = 'data/amine_bonds_A.csv'

    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file_path, sep=r'\s+')

    # print(data.head())
    # print(data.columns)

    # sort id1, id2 so that (id1, id2) and (id2, id1) are treated the same
    data[['id1', 'id2']] = pd.DataFrame(
        np.sort(data[['id1', 'id2']].values, axis=1),
        index=data.index
    )

    bond_lifetimes, bond_timesteps = compute_bond_lifetimes(data)

    print(f"Average bond lifetime: {np.mean([lifetime for lifetimes in bond_lifetimes.values() for lifetime in lifetimes]):.2f} timesteps")

    save_bond_lifetimes(bond_lifetimes, filename='data/bond_lifetimes_A.csv')

    # graph histogram of bond lifetimes
    plot_bond_lifetimes(bond_lifetimes, title='Bond Lifetimes Histogram (Type A)')

    csv_file_path = 'data/amine_bonds_C.csv'


    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file_path, sep=r'\s+')

    # print(data.head())
    # print(data.columns)

    # sort id1, id2 so that (id1, id2) and (id2, id1) are treated the same
    data[['id1', 'id2']] = pd.DataFrame(
        np.sort(data[['id1', 'id2']].values, axis=1),
        index=data.index
    )

    bond_lifetimes, bond_timesteps = compute_bond_lifetimes(data)

    print(f"Average bond lifetime: {np.mean([lifetime for lifetimes in bond_lifetimes.values() for lifetime in lifetimes]):.2f} timesteps")

    save_bond_lifetimes(bond_lifetimes, filename='data/bond_lifetimes_C.csv')

    # graph histogram of bond lifetimes
    plot_bond_lifetimes(bond_lifetimes, title='Bond Lifetimes Histogram (Type C)')


