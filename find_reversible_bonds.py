import argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import time
from scipy.spatial import cKDTree
import os


def snapshot_generator(traj_file):
    """Yields (timestep, data) tuples from a LAMMPS trajectory file without loading it all."""
    with open(traj_file) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                ts = int(f.readline().strip())
                # Skip next header lines
                while not (l := f.readline()).startswith("ITEM: NUMBER OF ATOMS"):
                    continue
                n_atoms = int(f.readline().strip())
                # Skip to atom data
                while not (l := f.readline()).startswith("ITEM: ATOMS"):
                    continue
                # Read atom data block
                data = np.loadtxt([f.readline() for _ in range(n_atoms)])
                yield ts, data


def process_timestep(snapshot, minimum_distance, box):
    """
    Finds all pairs within a minimum distance using cKDTree and periodic boundaries.
    """
    ts, data = snapshot
    coords = data[:, 2:5]
    ids = data[:, 0].astype(int)
    box = np.array(box)

    print(f"Processing timestep {ts}...")

    # Build extended periodic image positions for PBC search
    offsets = np.array([[i, j, k] for i in (-1, 0, 1)
                                  for j in (-1, 0, 1)
                                  for k in (-1, 0, 1)])
    extended_coords = np.vstack([coords + offset * box for offset in offsets])
    extended_ids = np.tile(ids, len(offsets))

    tree = cKDTree(extended_coords)
    pairs = tree.query_pairs(r=minimum_distance)

    results = []
    for ia, ib in pairs:
        id1, id2 = extended_ids[ia], extended_ids[ib]
        if id1 >= id2:
            continue
        p1, p2 = extended_coords[ia], extended_coords[ib]
        # Minimum image convention
        delta = p2 - p1 - np.round((p2 - p1) / box) * box
        dist = np.linalg.norm(delta)
        results.append([ts, id1, id2, *p1, *p2, dist])
    return results


def compute_reversible_bonds(traj_file, output_file, minimum_distance=0.5, box=(40, 40, 40), max_workers=4):
    """
    Stream over snapshots, compute pairwise distances in parallel, and stream results to CSV.
    """
    header = "timestep id1 id2 x1 y1 z1 x2 y2 z2 distance\n"
    with open(output_file, "w") as f:
        f.write(header)

    start_time = time()
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for snapshot in snapshot_generator(traj_file):
            futures.append(executor.submit(process_timestep, snapshot, minimum_distance, box))

        for future in as_completed(futures):
            results = future.result()
            if results:
                df = pd.DataFrame(results)
                df.to_csv(output_file, mode="a", index=False, header=False, sep=" ")
    print(f"Total time: {time() - start_time:.2f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute reversible bonds from LAMMPS trajectory")
    parser.add_argument("--traj", type=str, default="data/amine_positions.lammpstrj",
                        help="Input LAMMPS trajectory file")
    parser.add_argument("--output", type=str, default="data/amine_bonds.csv",
                        help="Output CSV file for bond data")
    parser.add_argument("--min_dist", type=float, default=0.6,
                        help="Minimum distance to consider a bond")
    parser.add_argument("--box", type=float, nargs=3, default=(40, 40, 40),
                        help="Simulation box dimensions (x y z)")
    parser.add_argument("--max_workers", type=int, default=4,
                        help="Maximum number of parallel workers")
    args = parser.parse_args()

    compute_reversible_bonds(
        args.traj,
        args.output,
        minimum_distance=args.min_dist,
        box=args.box,
        max_workers=args.max_workers
    )
