import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from time import time


def process_timestep(snapshot, minimum_distance, box):
    """
    Finds all pairs within a minimum distance and returns it a as a list of lists.

    Parameters:
        snapshot (tuple): (timestep, data) where data is a NxM numpy array
        minimum_distance (float): Distance cutoff for identifying reversible bonds
        box (tuple): Dimensions of the simulation box (x, y, z)

    Returns:
        list of lists: Each inner list contains [timestep, id1, id2, x1, y1, z1, x2, y2, z2, distance]
    """
    ts, data = snapshot
    ids = data[:,0].astype(int)
    coords = data[:,2:5]
    box = np.array(box)

    print("Processing timestep:", ts)

    delta = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    delta -= np.round(delta / box) * box
    dist_matrix = np.linalg.norm(delta, axis=-1)

    a_idx, b_idx = np.triu_indices(len(ids), k=1)
    mask = dist_matrix[a_idx, b_idx] < minimum_distance

    results = []
    for ia, ib in zip(a_idx[mask], b_idx[mask]):
        results.append([
            ts,
            ids[ia], ids[ib],
            coords[ia,0], coords[ia,1], coords[ia,2],
            coords[ib,0], coords[ib,1], coords[ib,2],
            dist_matrix[ia, ib]
        ])
    return results

def find_smallest_dist(snapshot, box):
    """
    Finds the smallest distance between all pairs of amine groups in a snapshot.

    (For testing purposes)

    Parameters:
        snapshot (tuple): (timestep, data) where data is a NxM numpy array
        box (tuple): Dimensions of the simulation box (x, y, z)

    Returns:
        tuple: (timestep, min_distance)
    """ 
    ts, data = snapshot
    print("Processing timestep:", ts)
    ids = data[:,0].astype(int)
    coords = data[:,2:5]
    box = np.array(box)

    delta = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    delta -= np.round(delta / box) * box
    dist_matrix = np.linalg.norm(delta, axis=-1)

    a_idx, b_idx = np.triu_indices(len(ids), k=1)
    min_dist = np.min(dist_matrix[a_idx, b_idx])
    print("Minimum distance at timestep", ts, "is", min_dist)
    return ts, min_dist
    

def compute_reversible_bonds(traj_file, minimum_distance=0.5, box=(40,40,40), max_workers=4):
    """
    Compute pairwise distances between amine groups from a LAMMPS trajectory file and return as a DataFrame.
    
    Parameters:
        traj_file (str): Path to the LAMMPS trajectory file
        minimum_distance (float): Distance cutoff for identifying reversible bonds
        box (tuple): Dimensions of the simulation box (x, y, z)
    
    Returns:
        pd.DataFrame: Columns = timestep, id1, id2, x1, y1, z1, x2, y2, z2, distance
    """
    # read snapshots first
    snapshots = []
    with open(traj_file) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1])
            n_atoms = int(lines[i+3])
            data_start = i + 9
            data_end = data_start + n_atoms
            data = np.loadtxt(lines[data_start:data_end])
            snapshots.append((ts, data))
            i = data_end
        else:
            i += 1

    # parallel processing
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_timestep, snap, minimum_distance, box) for snap in snapshots]
        # futures = [executor.submit(find_smallest_dist, snap, box) for snap in snapshots]
        for future in futures:
            results.extend(future.result())
            # results.append(future.result())

    df = pd.DataFrame(
        results,
        columns=["timestep","id1","id2","x1","y1","z1","x2","y2","z2","distance"]
    )
    # df = pd.DataFrame(
    #     results,
    #     columns=["timestep","min_distance"]
    # )
    return df


if __name__ == "__main__":
    # Example usage
    start = time()
    df = compute_reversible_bonds("data/amine_positions.lammpstrj", minimum_distance=(0.19*2), box=(40,40,40), max_workers=4)
    df.to_csv("data/amine_bonds.csv", index=False, sep=' ')
    print(f"Total time: {time() - start:.2f} seconds")
    print(df.head())
