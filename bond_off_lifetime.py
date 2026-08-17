import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time

cutoffs = [0.4, 0.5, 0.6, 0.7]

# As = [60, 70, 80, 90, 100, 110, 120, 130, 140]
As = [140]


spacing = 5

# in_dir = "data/patches"
in_dir = "/mnt/hdd1/patchy_polymers_files"
out_dir = "data"


def get_current_bonds(file_iter, N_patches, pairs=False):
    current_bonds = set() # track current bonded pairs in this frame to compare with bonded_pairs dict
    clusters = pd.read_csv(file_iter, delimiter='\s+', header=None, nrows=N_patches, usecols=[0,5,6], names=['id', 'patch_coord', 'patch_cluster']).to_numpy()
    clusters = clusters[clusters[:, 1] > 0] # filter only consider bonds with at least 1 patch (patch_coord is in column 5)
    bonded_clusters = clusters[clusters[:, 2].argsort()] # sort by patch_cluster


    for i in range(len(bonded_clusters)-1):
        id1, id2 = int(bonded_clusters[i, 0]), int(bonded_clusters[i+1, 0]) # get ids of current and next cluster
        if int(bonded_clusters[i, 2]) == int(bonded_clusters[i+1, 2]): # if same patch_cluster, they are bonded
            if pairs:
                ids = (id1, id2) if id1 < id2 else (id2, id1) # order ids to avoid duplicates
                current_bonds.add(ids) # add to current bonds set
            else:
                current_bonds.add(id1) # add to current bonds set
                current_bonds.add(id2) # add to current bonds set
    
    return current_bonds

def read_frame(f, timestep, pairs=True):
    while True:
        # stop reading if file is empty
        line = f.readline()
        if not line:
            return None, None
        if line.startswith("ITEM: TIMESTEP"): # timestep is on next line
            time = float(f.readline().strip()) * timestep # convert to time units
            # print(f"Processing timestep: {time / timestep:.0f} (time units: {time})")
        elif line.startswith("ITEM: NUMBER OF ATOMS"): # number of entries is on next line
            N_patches = int(f.readline().strip())
        elif line.startswith("ITEM: ATOMS"): # data starts on next line
            current_bonds = get_current_bonds(f, N_patches, pairs=pairs)
            return time, current_bonds


def process_buffered_lifetimes(cutoffs, A, timestep=0.002):
    """
    Process lifetimes from clustering data file. Each frame of the file contains the following sections:
    - "ITEM: TIMESTEP" followed by the timestep on the next line
    - "ITEM: NUMBER OF ATOMS" followed by the number of entries on the next line
    - "ITEM: ATOMS" followed by the data lines for each entry, with 7 values per line: id, mol, x, y, z, patch_coord, patch_cluster

    The bonded_pairs list is a dictionary where the key is a tuple of (id1, id2) with id1 < id2 to avoid duplicates, and the value is the 
    start_time of the bond. 

    The bond_lifetimes list is a list of tuples (id1, id2, start_time, end_time, lifetime) for each bonded pair that has broken.

    If there are any bonded_pairs that have not broken by the end of the dump, we ignore them and do not record a lifetime for them.

    All times for breaking/creating are considered on the timestep we note them, so we should on average get the right time 
    (i.e. both events are noted AFTER they happen)

    Args:
        cutoffs: list of cutoff distances to process
        A: Gaussian potential strength, used to find files
        timestep: timestep size in LJ units (default is 0.002)
    Returns:
        bond_lifetimes: dict of lists of tuples (id1, id2, start_time, end_time, lifetime) for each bonded pair, keys are cutoff pairs (cutoff_on, cutoff_off)
    """
    files = {
        cutoff: open(f"{in_dir}/{cutoff}prodpatch_chainlength_31_patchspacing_{spacing}_N_beads_24800_gaussA_{A}_r0_0.25_gaussB_10_seed_12345.bin.txt", 'r') for cutoff in cutoffs
    }
    cutoff_pairs = [(cutoff_on, cutoff_off) for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off]

    # init empty bond lifetime list for each cutoff. lookup by (cutoff_on, cutoff_off
    bond_lifetimes = {(cutoff_on, cutoff_off): [] for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off}

    bonded_pairs = {(cutoff_on, cutoff_off): {} for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off} # init empty bonded pair dict for each cutoff. lookup by (id1, id2) where id1 < id2 to avoid duplicates. value is start_time of the bond

    frames = {cutoff: None for cutoff in cutoffs}

    while True:
        for cutoff in cutoffs:
            frames[cutoff] = read_frame(files[cutoff], timestep, pairs=True)
            # print(f"Read frame for cutoff {cutoff}: time {frames[cutoff][0]}, number of bonds {len(frames[cutoff][1]) if frames[cutoff][1] is not None else 0}")
            if frames[cutoff][0] is None:
                for f in files.values():
                    f.close() # idk if we really need this, but better to be safe
                return bond_lifetimes # finish if we reach the end of any file
        for cutoff_on, cutoff_off in cutoff_pairs:
            curr_bonded_pairs = bonded_pairs[(cutoff_on, cutoff_off)]
            time_on, current_bonds = frames[cutoff_on]
            for ids in current_bonds:
                if ids not in curr_bonded_pairs: # if new bond, add to bonded_pairs dict with start_time as current timestep
                    curr_bonded_pairs[ids] = time_on

            time_off, current_bonds = frames[cutoff_off]
            assert time_on == time_off, f"Time mismatch between cutoffs: {time_on} != {time_off}" # sanity check that both cutoffs are at the same timestep
            broken_bonds = curr_bonded_pairs.keys() - current_bonds # bonds that are in bonded_pairs but not in current_bonds are broken
            for ids in broken_bonds:
                start_time = curr_bonded_pairs.pop(ids) # remove from bonded_pairs and get start_time
                lifetime = time_off - start_time
                if (lifetime > timestep):
                    bond_lifetimes[(cutoff_on, cutoff_off)].append((*ids, start_time, time_off, lifetime)) # add to bond_lifetimes list
        if time_on % 1.0 < timestep * 10: # print every 1.0 time units
            print(f"Processed frame at time {time_on}")


def process_buffered_offtimes(cutoffs, A, timestep=0.002, N_patches=4000):
    """
    Process offtimes from clustering data file. Each frame of the file contains the following sections:
    - "ITEM: TIMESTEP" followed by the timestep on the next line
    - "ITEM: NUMBER OF ATOMS" followed by the number of entries on the next line
    - "ITEM: ATOMS" followed by the data lines for each entry, with 7 values per line: id, mol, x, y, z, patch_coord, patch_cluster

    The unbonded_patches list is a dictionary where the key is a patch's id and the value is the last timestep when a bond it was part of breaks.

    The bond_offtimes list is a list of tuples (id, start_time, end_time, offtime) for each bonded pair that has formed.

    If there are any bonded_pairs that have not created by the end of the dump, we ignore them and do not record a offtime for them.

    At the start of the simulation, all patches are unbonded, so we do not start any timers for them. Therefore, it takes at least one creation and breaking event 
    for bond offtimes to start being recorded, might take some time to get to a steady state.

    All times for breaking/creating are considered on the timestep we note them, so we should on average get the right time 
    (i.e. both events are noted AFTER they happen).

    Args:
        cutoffs: list of cutoff distances to process
        A: Gaussian potential strength, used to find files
        timestep: timestep size in LJ units (default is 0.002)
        N_patches: number of patches, (default is 4000, from N_b=31, M=5)
    Returns:
        - bond_offtimes: list of tuples (id, start_time, end_time, offtime) for each patch
    """

    all_ids = set(list(range(1, N_patches+1)))

    files = {
        cutoff: open(f"{in_dir}/{cutoff}prodpatch_chainlength_31_patchspacing_{spacing}_N_beads_24800_gaussA_{A}_r0_0.25_gaussB_10_seed_12345.bin.txt", 'r') for cutoff in cutoffs
    }
    cutoff_pairs = [(cutoff_on, cutoff_off) for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off]

    # init empty bond offtime list for each cutoff. lookup by (cutoff_on, cutoff_off
    bond_offtimes = {(cutoff_on, cutoff_off): [] for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off}

    unbonded_ids = {(cutoff_on, cutoff_off): {} for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off} # init empty unbonded ids dict for each cutoff

    frames = {cutoff: None for cutoff in cutoffs}

    while True:
        for cutoff in cutoffs:
            frames[cutoff] = read_frame(files[cutoff], timestep, pairs=False)
            # print(f"Read frame for cutoff {cutoff}: time {frames[cutoff][0]}, number of bonds {len(frames[cutoff][1]) if frames[cutoff][1] is not None else 0}")
            if frames[cutoff][0] is None:
                for f in files.values():
                    f.close() # idk if we really need this, but better to be safe
                return bond_offtimes # finish if we reach the end of any file
        for cutoff_on, cutoff_off in cutoff_pairs:
            curr_unbonded_ids = unbonded_ids[(cutoff_on, cutoff_off)]
            time_on, current_bonds = frames[cutoff_on]
            for id in current_bonds:
                if id in curr_unbonded_ids: # if just formed, stop their timers
                    start_time = curr_unbonded_ids.pop(id)
                    bond_offtimes[(cutoff_on, cutoff_off)].append((id, start_time, time_on, time_on - start_time))

            time_off, current_bonds = frames[cutoff_off]
            assert time_off == time_on, f"Time mismatch between two cutoffs: {time_on} != {time_off}"
            broken_bonds = all_ids.difference(curr_unbonded_ids, current_bonds) # bonds that are in unbonded_ids but not in current_bonds are broken
            for id in broken_bonds: # just broken, set start of offtime and track timers
                unbonded_ids[(cutoff_on, cutoff_off)][int(id)] = time_on
        if time_on % 1.0 < timestep * 10: # print every 1.0 time units
                    print(f"Processed frame at time {time_on}")


def process_and_save_lifetimes(cutoffs, A):
    
    start = time()
    # filename_on = f"{in_dir}/{cutoff_on}prodpatch_chainlength_31_patchspacing_{spacing}_N_beads_24800_gaussA_{A}_r0_0.25_gaussB_10_seed_12345.bin.txt"
    # filename_off = f"{in_dir}/{cutoff_off}prodpatch_chainlength_31_patchspacing_{spacing}_N_beads_24800_gaussA_{A}_r0_0.25_gaussB_10_seed_12345.bin.txt"
    lifetimes = process_buffered_lifetimes(cutoffs, A, 0.002)
    for cutoff_on, cutoff_off in [(cutoff_on, cutoff_off) for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off]:
        lifetimes_df = pd.DataFrame(lifetimes[(cutoff_on, cutoff_off)], columns=["id1", "id2", "start_time", "end_time", "lifetime"])
        lifetimes_df.to_csv(f"{out_dir}/lifetimes/buffered_lifetimes_A{A}_cutoff_on{cutoff_on}_cutoff_off{cutoff_off}_r0_0.25_spacing_{spacing}.csv", index=False)
        print(f"Processed lifetimes for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}, total lifetimes recorded: {len(lifetimes)} in {time() - start:.2f} seconds")



def process_and_save_offtimes(cutoffs, A):
    start = time()
    # filename_on = f"{in_dir}/{cutoff_on}prodpatch_chainlength_31_patchspacing_{spacing}_N_beads_24800_gaussA_{A}_r0_0.25_gaussB_10_seed_12345.bin.txt"
    # filename_off = f"{in_dir}/{cutoff_off}prodpatch_chainlength_31_patchspacing_{spacing}_N_beads_24800_gaussA_{A}_r0_0.25_gaussB_10_seed_12345.bin.txt"
    offtimes = process_buffered_offtimes(cutoffs, A, 0.002)
    for cutoff_on, cutoff_off in [(cutoff_on, cutoff_off) for cutoff_on in cutoffs for cutoff_off in cutoffs if cutoff_on <= cutoff_off]:
        offtimes_df = pd.DataFrame(offtimes[(cutoff_on, cutoff_off)], columns=["id", "start_time", "end_time", "offtime"])
        offtimes_df.to_csv(f"{out_dir}/offtimes/buffered_offtimes_A{A}_cutoff_on{cutoff_on}_cutoff_off{cutoff_off}_r0_0.25_spacing_{spacing}.csv", index=False)
        print(f"Processed offtimes for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}, total offtimes recorded: {len(offtimes)} in {time() - start:.2f} seconds")



if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument('--in_dir', type=str, default=in_dir, help='Input directory')
    args.add_argument('--out_dir', type=str, default=out_dir, help='Output directory')
    args.add_argument('--cutoffs', type=float, nargs='*', default=cutoffs, help='Cutoff distances to process')
    args.add_argument('--As', type=int, nargs='*', default=As, help='A values to process')
    args.add_argument('--spacing', type=int, default=spacing, help='Patch to patch spacing')
    args = args.parse_args()

    print(args.cutoffs)

    in_dir = args.in_dir
    out_dir = args.out_dir
    cutoffs = args.cutoffs
    As = args.As
    spacing = args.spacing
    # print(cutoffs)

    for A in As:
        print(f"Processing lifetimes for A={A}")
        lifetimes = process_and_save_lifetimes(args.cutoffs, A)
        print(f"Processing offtimes for A={A}")
        offtimes = process_and_save_offtimes(args.cutoffs, A)
