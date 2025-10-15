import argparse
import pandas as pd
import csv
from tqdm import tqdm

def find_multiplets_to_csv(data, output_csv):
    """
    Stream groups of atoms with more than 2 members to CSV.
    Each group has the center and all its bonded neighbors.
    """
    # print(data.columns[0])
    # input()
    with open(output_csv, mode="w", newline="") as f:
        writer = csv.writer(f)
        # Header: timestep, atom IDs, coordinates
        writer.writerow(["timestep", "center"] + 
                        ["neighbor_id", "neighbor_x", "neighbor_y", "neighbor_z"] * 10)  # max 10 neighbors for example

        for timestep, bonds in tqdm(data.groupby("timestep"), total=data['timestep'].nunique(), desc="Finding groups"):
            adjacency = {}
            center_coords = {}

            # Build adjacency and store center coords
            for _, row in bonds.iterrows():
                i, j = row['id1'], row['id2']

                adjacency.setdefault(i, []).append((j, row['x2'], row['y2'], row['z2']))
                adjacency.setdefault(j, []).append((i, row['x1'], row['y1'], row['z1']))

                if i not in center_coords:
                    center_coords[i] = (row['x1'], row['y1'], row['z1'])
                if j not in center_coords:
                    center_coords[j] = (row['x2'], row['y2'], row['z2'])

            seen_groups = set()

            for center, neighbors in adjacency.items():
                # if len(neighbors) < 2:
                #     continue

                # Create a sorted tuple of all atom IDs to ensure uniqueness
                group_ids = tuple(sorted([center] + [n[0] for n in neighbors]))
                if group_ids in seen_groups:
                    continue
                seen_groups.add(group_ids)

                row_data = [timestep, center]
                for n_id, nx, ny, nz in neighbors:
                    row_data.extend([n_id, nx, ny, nz])
                writer.writerow(row_data)

# include a command line input to change the name of the output file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find multiplets in bond data")
    parser.add_argument("--output", type=str, default="data/amine_multiplets.csv",
                        help="Output CSV file for multiplets")
    parser.add_argument("--input", type=str, default="data/amine_bonds.csv",
                        help="Input CSV file with bond data")
    args = parser.parse_args()
    output_csv = args.output
    
    # Load bond data
    data = pd.read_csv(args.input, delimiter=' ')

    # Find multiplets and write to CSV
    find_multiplets_to_csv(data, output_csv)
