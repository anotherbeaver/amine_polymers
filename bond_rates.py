import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

cutoffs = [0.4, 0.5, 0.6, 0.7]
As = [60, 70, 80, 90, 100, 110, 120, 130, 140]
spacing = 7

# in_dir = "data/patches"
in_dir = "/mnt/hdd1/patchy_polymers_files"
out_dir = "data"

def biexponential(x, a1, b1, a2, b2):
    return a1 * np.exp(-b1 * x) + a2 * np.exp(-b2 * x)

def single_exponential(x, a, b):
    return a * np.exp(-b * x)


if __name__ == "__main__":
    import argparse

    args = argparse.ArgumentParser()
    args.add_argument('--in_dir', type=str, default=in_dir, help='Input directory')
    args.add_argument('--out_dir', type=str, default=out_dir, help='Output directory')
    args.add_argument('--cutoffs', type=float, nargs='*', default=cutoffs, help='Cutoff distances to process')
    args.add_argument('--As', type=int, nargs='*', default=As, help='A values to process')
    args.add_argument('--spacing', type=int, default=spacing, help='Patch to patch spacing')
    args = args.parse_args()

    in_dir = args.in_dir
    out_dir = args.out_dir
    cutoffs = args.cutoffs
    As = args.As
    spacing = args.spacing

    # tau_list = []

    # cutoff_on_list = [0.4, 0.5, 0.6, 0.7]
    # cutoff_off_list = [0.5, 0.6, 0.7]

    cutoff_on_list = cutoffs
    cutoff_off_list = cutoffs

    for A in As:
        for cutoff_off in cutoff_off_list:
            # fig = plt.figure(figsize=(8, 6))
            for cutoff_on in cutoff_on_list:
                if cutoff_on <= cutoff_off:
                
                    df = pd.read_csv(f"{out_dir}/lifetimes/buffered_lifetimes_A{A}_cutoff_on{cutoff_on}_cutoff_off{cutoff_off}_r0_0.25_spacing_{spacing}.csv")
                    lifetimes = np.array(df["lifetime"].values)

                    # values, counts = np.unique(lifetimes, return_counts=True)
                    hist, bins = np.histogram(lifetimes, bins=50)
                    # print(f"Processed lifetimes for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}, total lifetimes recorded: {len(lifetimes)}")
                    hist = hist[1:] # remove first bin
                    bins = bins[1:] # remove first bin
                    # print(f"Removing, remaining lifetimes recorded: {np.sum(hist)}")
                    try:
                        coeffs, cov = curve_fit(biexponential, bins[:-1], hist, p0=[1, 0.01, 1, 0.01], maxfev=100000000)
                        # print(f"Fitted coefficients for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}: {coeffs}")
                        tau_1, tau_2 = sorted((1/coeffs[1], 1/coeffs[3]))
                        # print(f"Fitted lifetimes for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}: tau_1={tau_1:.2f}, tau_2={tau_2:.2f}")

                    except:
                        print(f"Failed to fit biexponential for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}")
                        tau_1, tau_2 = np.nan, np.nan

                    tau_list = [A, "life", cutoff_on, cutoff_off, tau_1, tau_2]
                    df = pd.DataFrame([tau_list], columns=["A", "off/life", "cutoff_on", "cutoff_off", "tau_1", "tau_2"])
                    df.to_csv(f"{out_dir}/lifetimes/tau_lifetimes_A{A}_cutoff_on{cutoff_on}_cutoff_off{cutoff_off}_r0_0.25_spacing_{spacing}.csv", index=False)



                    # plt.plot(bins[:-1], hist, alpha=0.5, label=f"A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}")
                    # plt.plot(bins[:-1], biexponential(bins[:-1], *coeffs), '--', label=r"Fit $\tau_1=$" + f"{tau_1:.2f}" + r", $\tau_2=$" + f"{tau_2:.2f}, lifetimes")
                    
            # plt.xlabel("Lifetime (LJ units)")
            # plt.ylabel("Density")
            # plt.yscale("log")
            # # plt.ylim(0, 0.2)
            # plt.grid()
            # plt.title(f"Distribution of Bond Lifetimes Assorted r_on, r_off={cutoff_off} for A={A}")
            # plt.legend()
            # plt.show()


    
    # tau_list = []


    for A in As:
        for cutoff_on in cutoff_on_list:
            # fig = plt.figure(figsize=(8, 6))
            for cutoff_off in cutoff_off_list:
                if cutoff_on <= cutoff_off:
                
                    df = pd.read_csv(f"{out_dir}/offtimes/buffered_offtimes_A{A}_cutoff_on{cutoff_on}_cutoff_off{cutoff_off}_r0_0.25_spacing_{spacing}.csv")
                    offtimes = np.array(df["offtime"].values)

                    # values, counts = np.unique(lifetimes, return_counts=True)
                    hist, bins = np.histogram(offtimes, bins=50)
                    # print(f"Processed offtimes for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}, total lifetimes recorded: {len(offtimes)}")
                    hist = hist[1:] # remove first bin
                    bins = bins[1:] # remove first bin
                    # print(f"Removing, remaining offtimes recorded: {np.sum(hist)}")

                    try:
                        coeffs, cov = curve_fit(biexponential, bins[:-1], hist, p0=[1, 0.01, 1, 0.01], maxfev=100000000)
                        # print(f"Fitted coefficients for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}: {coeffs}")
                        tau_1, tau_2 = sorted((1/coeffs[1], 1/coeffs[3]))
                        # print(f"Fitted offtimes for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}: tau_1={tau_1:.2f}, tau_2={tau_2:.2f}")

                    except:
                        print(f"Failed to fit biexponential for A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}")
                        tau_1, tau_2 = np.nan, np.nan

                    tau_list = [A, "off", cutoff_on, cutoff_off, tau_1, tau_2]
                    df = pd.DataFrame([tau_list], columns=["A", "off/life", "cutoff_on", "cutoff_off", "tau_1", "tau_2"])
                    df.to_csv(f"{out_dir}/offtimes/tau_offtimes_A{A}_cutoff_on{cutoff_on}_cutoff_off{cutoff_off}_r0_0.25_spacing_{spacing}.csv", index=False)
                    

                    # plt.plot(bins[:-1], hist, alpha=0.5, label=f"A={A}, cutoff_on={cutoff_on}, cutoff_off={cutoff_off}")
                    # plt.plot(bins[:-1], biexponential(bins[:-1], *coeffs), '--', label=r"Fit $\tau_1=$" + f"{tau_1:.2f}" + r", $\tau_2=$" + f"{tau_2:.2f}, offtimes")

            # plt.xlabel("Lifetime (LJ units)")
            # plt.ylabel("Density")
            # plt.yscale("log")
            # # plt.ylim(0, 0.2)
            # plt.grid()
            # plt.title(f"Distribution of Bond Offtimes Assorted r_on, r_off={cutoff_off} for A={A}")
            # plt.legend()
            # plt.show()




    # life_tau_df = pd.read_csv(f"{out_dir}/lifetimes/tau_lifetimes.csv")
    # off_tau_df = pd.read_csv(f"{out_dir}/offtimes/tau_offtimes.csv")

    # # plot tau_2 for life over tau_2 for off for each A, cutoff_on, cutoff_off
    # plt.figure(figsize=(8, 6))
    # As = [80, 110]

    # for A in As:
    #     for cutoff_off in cutoff_off_list:
    #         tau_life_over_off_list = []
    #         saved_cutoff_on_list = []
    #         for cutoff_on in cutoff_on_list:
    #             if cutoff_on <= cutoff_off:
    #                 tau_life = life_tau_df.loc[(life_tau_df["A"] == A) & (life_tau_df["off/life"] == "life") & (life_tau_df["cutoff_on"] == cutoff_on) & (life_tau_df["cutoff_off"] == cutoff_off), "tau_2"].values[0]
    #                 tau_off = off_tau_df.loc[(off_tau_df["A"] == A) & (off_tau_df["off/life"] == "off") & (off_tau_df["cutoff_on"] == cutoff_on) & (off_tau_df["cutoff_off"] == cutoff_off), "tau_2"].values[0]
    #                 tau_life_over_off_list.append(tau_life / tau_off)
    #                 print(tau_life_over_off_list)
    #                 saved_cutoff_on_list.append(cutoff_on)

    #         plt.plot(saved_cutoff_on_list, np.array(tau_life_over_off_list), label=f"A={A}, cutoff_on={cutoff_on}", alpha=0.7)
    #         plt.xlabel(r"$r_{\text{on}}$ for (LJ distance units)")
    #         plt.ylabel(r"$\tau_{life}/\tau_{off}$ (LJ time units)")
    #         plt.title(f"Comparison of \\alpha/\\beta for Bond Lifetimes and Offtimes, cutoff_off={cutoff_off}")
    #         # plt.xscale("log")
    #         # plt.yscale("log")
    #         plt.grid()
    #         plt.legend()
    #         plt.show()

