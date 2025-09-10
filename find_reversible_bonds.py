import numpy as np
import pandas as pd

def compute_reversible_bonds(traj_file, minimum_distance=0.5, box=(100,100,100)):
    """
    Compute pairwise distances between amine groups from a LAMMPS trajectory file and return as a DataFrame.

    Parameters:
        traj_file (str): Path to the LAMMPS trajectory file
        box (tuple): Dimensions of the simulation box (x, y, z)

    Returns:
        pd.DataFrame: DataFrame containing timestep, atom IDs, coordinates, and distances
    """
    from itertools import combinations
    
    with open(traj_file, "r") as f:
        lines = f.readlines()
    
    timesteps, snapshots = [], []
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            ts = int(lines[i+1])
            print("Reading timestep", ts)
            n_atoms = int(lines[i+3])
            box_bounds = [list(map(float, lines[i+5+j].split())) for j in range(3)]
            data_start = i+9
            data_end = data_start + n_atoms
            data = np.loadtxt(lines[data_start:data_end])
            # data columns: id, type, x, y, z
            snapshots.append((ts, data))
            i = data_end
        else:
            i += 1
    
    results = []
    for ts, data in snapshots:
        print("Processing timestep", ts)
        # print("Data head:\n", data[:5])
        print("Total pairs before filtering:", len(data)*(len(data)-1)//2)
        ids = data[:,0].astype(int)
        coords = data[:,2:5]

        # first perform a quick filter based on coordinate differences (x,y,z)
        for (i1, i2) in combinations(range(len(ids)), 2):
            if ids[i1] > ids[i2]: # prevent double counting
                continue
            # only consider distance in x,y,z less than min
            if coords[i1,0] - coords[i2,0] < minimum_distance and \
               coords[i1,1] - coords[i2,1] < minimum_distance and \
               coords[i1,2] - coords[i2,2] < minimum_distance:
                results.append([ts, ids[i1], ids[i2], *coords[i1], *coords[i2]])
        
        print('Found preliminary', len(results), 'pairs')
        
        # only keep results with norm less than minimum_distance (because norm is expensive)
        for row in results:
            dist = np.linalg.norm(np.array(row[3:6]) - np.array(row[6:9]))
            if dist < minimum_distance:
                row.append(dist)
            else:
                row.append(np.nan)

        print('After distance check', len([row for row in results if not np.isnan(row[-1])]), 'pairs remain')

        # clean up results to only keep valid distances
        results = [row for row in results if not np.isnan(row[-1])]
        print(f"Found {len(results)} pairs in timestep {ts}")
    
    df = pd.DataFrame(results, columns=["timestep", "id1", "id2", "x1", "y1", "z1", "x2", "y2", "z2", "distance"])
        
    return df

# Example usage
df = compute_reversible_bonds("data/amine_positions.lammpstrj", box=(40,40,40))
print(df.head())
