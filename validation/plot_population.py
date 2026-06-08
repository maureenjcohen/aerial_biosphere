"""
plot_population.py
Plot population time series from ensemble pop files to verify stationarity
(ADF test analog, Yates Section 2.3).

Reads model/output/ensemble_NNN_pop.dat and plots:
  - n_orgs vs time for all ensemble members
  - Mean and std-dev band across ensemble
  - Flag members that went to zero (extinction)

Usage:
    python plot_population.py [--outdir model/output] [--n_ens 20]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_pop_series(outdir, n_ens):
    series = []
    for iens in range(1, n_ens + 1):
        fname = os.path.join(outdir, f"ensemble_{iens:03d}_pop.dat")
        if not os.path.isfile(fname):
            continue
        data = np.loadtxt(fname, comments="#")
        if data.ndim == 1 or data.shape[0] == 0:
            continue
        series.append((iens, data[:, 0], data[:, 1].astype(int)))
    return series


def main():
    parser = argparse.ArgumentParser(description="Plot population time series")
    parser.add_argument("--outdir", default="../model/output")
    parser.add_argument("--n_ens", type=int, default=20)
    parser.add_argument("--savefig", default="population_timeseries.png")
    args = parser.parse_args()

    series = load_pop_series(args.outdir, args.n_ens)
    if not series:
        sys.exit("No pop files found")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    t_ref = series[0][1]
    all_n = []
    extinct = []
    for iens, t, n in series:
        color = "steelblue" if np.any(n == 0) else "steelblue"
        alpha = 0.3
        ax1.plot(t / 365.25, n, lw=0.7, alpha=alpha, color="steelblue")
        all_n.append(np.interp(t_ref, t, n))
        if n[-1] == 0:
            extinct.append(iens)

    arr = np.array(all_n)
    mu = arr.mean(axis=0)
    sigma = arr.std(axis=0)
    ax1.plot(t_ref / 365.25, mu, "k-", lw=2, label="ensemble mean")
    ax1.fill_between(t_ref / 365.25, mu - sigma, mu + sigma, alpha=0.25, color="k",
                     label="±1 σ")
    ax1.set_ylabel("n_orgs")
    ax1.set_title(f"Population time series ({len(series)} members)")
    ax1.legend(fontsize=9)
    if extinct:
        ax1.set_title(ax1.get_title() + f" — EXTINCTION in members: {extinct}")

    # Last 75 years: check stationarity (Yates ADF approach)
    t_max = t_ref.max() / 365.25
    t75_start = t_max - 75
    mask = t_ref / 365.25 >= t75_start
    ax2.plot(t_ref[mask] / 365.25, mu[mask], "k-", lw=2)
    ax2.fill_between(t_ref[mask] / 365.25,
                     (mu - sigma)[mask], (mu + sigma)[mask],
                     alpha=0.25, color="k")
    ax2.set_xlabel("Time (Earth years)")
    ax2.set_ylabel("n_orgs")
    ax2.set_title("Last 75 years (stationarity window, Yates Sect. 2.3)")

    plt.tight_layout()
    fig.savefig(args.savefig, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.savefig}")
    print(f"Extinct members: {extinct if extinct else 'none'}")
    print(f"Final mean population: {mu[-1]:.0f} ± {sigma[-1]:.0f}")


if __name__ == "__main__":
    main()
