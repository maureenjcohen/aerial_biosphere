"""
plot_distributions.py
Reproduce Yates et al. (2017) Figures 3, 4, 5 from IBM ensemble output.

Reads state files from model/output/ensemble_NNN_state.dat and plots:
  - Fig 3 analog: organism mass distribution (log scale)
  - Fig 4 analog: organism radius distribution
  - Fig 5 analog: organism altitude distribution

Usage:
    python plot_distributions.py [--outdir model/output] [--n_ens 20]
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Column indices in state files
# # z[m]  mass[kg]  radius[m]  G[-]  rho_org[kg/m3]  age[days]  skin_width[m]
# ---------------------------------------------------------------------------
COL_Z        = 0
COL_MASS     = 1
COL_RADIUS   = 2
COL_G        = 3
COL_RHO      = 4
COL_AGE      = 5
COL_SKIN     = 6


def load_ensemble(outdir, n_ens):
    """Load all state files; return concatenated arrays."""
    all_z, all_mass, all_radius, all_G, all_rho = [], [], [], [], []
    loaded = 0
    for iens in range(1, n_ens + 1):
        fname = os.path.join(outdir, f"ensemble_{iens:03d}_state.dat")
        if not os.path.isfile(fname):
            print(f"  [skip] {fname} not found")
            continue
        data = np.loadtxt(fname, comments="#")
        if data.ndim == 1:
            data = data[np.newaxis, :]
        if data.shape[0] == 0:
            print(f"  [skip] {fname} empty")
            continue
        all_z.append(data[:, COL_Z])
        all_mass.append(data[:, COL_MASS])
        all_radius.append(data[:, COL_RADIUS])
        all_G.append(data[:, COL_G])
        all_rho.append(data[:, COL_RHO])
        loaded += 1
        print(f"  loaded {fname}: {data.shape[0]} organisms")

    if loaded == 0:
        sys.exit("No state files found — has the model run completed?")

    return (
        np.concatenate(all_z),
        np.concatenate(all_mass),
        np.concatenate(all_radius),
        np.concatenate(all_G),
        np.concatenate(all_rho),
        loaded,
    )


def plot_mass_distribution(mass, ax):
    """Yates Fig 3: mass distribution on log scale."""
    log_m = np.log10(mass)
    bins = np.linspace(log_m.min() - 0.1, log_m.max() + 0.1, 50)
    ax.hist(log_m, bins=bins, density=True, color="steelblue", edgecolor="none", alpha=0.8)
    ax.axvline(np.log10(np.median(mass)), color="red", ls="--", lw=1.5,
               label=f"median = {np.median(mass):.2e} kg")
    ax.set_xlabel("log$_{10}$(mass / kg)")
    ax.set_ylabel("Probability density")
    ax.set_title("Mass distribution (Yates Fig. 3 analog)")
    ax.legend(fontsize=9)


def plot_radius_distribution(radius, ax):
    """Yates Fig 4: radius distribution."""
    r_um = radius * 1e6   # convert m -> μm
    bins = np.linspace(0, r_um.max() * 1.05, 60)
    ax.hist(r_um, bins=bins, density=True, color="darkorange", edgecolor="none", alpha=0.8)
    ax.axvline(np.median(r_um), color="red", ls="--", lw=1.5,
               label=f"median = {np.median(r_um):.1f} μm")
    ax.set_xlabel("Radius (μm)")
    ax.set_ylabel("Probability density")
    ax.set_title("Radius distribution (Yates Fig. 4 analog)")
    ax.legend(fontsize=9)


def plot_altitude_distribution(z, ax):
    """Yates Fig 5: altitude distribution."""
    z_km = z * 1e-3   # convert m -> km
    bins = np.linspace(0, 105, 53)
    ax.hist(z_km, bins=bins, density=True, color="seagreen", edgecolor="none", alpha=0.8,
            orientation="horizontal")
    ax.axhline(np.median(z_km), color="red", ls="--", lw=1.5,
               label=f"median = {np.median(z_km):.1f} km")
    ax.set_ylabel("Altitude (km)")
    ax.set_xlabel("Probability density")
    ax.set_title("Altitude distribution (Yates Fig. 5 analog)")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9)


def plot_G_distribution(G, ax):
    """G distribution: should skew toward low G (Yates Sect. 3.2)."""
    bins = np.linspace(0, 1, 51)
    ax.hist(G, bins=bins, density=True, color="purple", edgecolor="none", alpha=0.8)
    ax.axvline(np.median(G), color="red", ls="--", lw=1.5,
               label=f"median G = {np.median(G):.3f}")
    ax.set_xlabel("G (growth strategy, 0=solid, 1=balloon)")
    ax.set_ylabel("Probability density")
    ax.set_title("G distribution (Yates Sect. 3.2: skew toward G→0)")
    ax.legend(fontsize=9)


def print_summary(z, mass, radius, G, rho, n_ens):
    print("\n=== Summary statistics ===")
    print(f"  Ensemble members with data: {n_ens}")
    print(f"  Total organisms (snapshot): {len(mass)}")
    print(f"  Mass:    median={np.median(mass):.3e} kg,  mean={np.mean(mass):.3e} kg")
    print(f"           10/90th pct = [{np.percentile(mass,10):.2e}, {np.percentile(mass,90):.2e}] kg")
    print(f"  Radius:  median={np.median(radius)*1e6:.2f} μm,  mean={np.mean(radius)*1e6:.2f} μm")
    print(f"  G:       median={np.median(G):.3f},  mean={np.mean(G):.3f}  (Yates: G→0)")
    print(f"  rho_org: median={np.median(rho):.0f} kg/m3")
    print(f"  Altitude: median={np.median(z)*1e-3:.1f} km,  "
          f"10/90th=[{np.percentile(z,10)*1e-3:.1f}, {np.percentile(z,90)*1e-3:.1f}] km")
    print()


def main():
    parser = argparse.ArgumentParser(description="Plot Yates 2017 IBM validation figures")
    parser.add_argument("--outdir", default="../model/output",
                        help="Path to model output directory (default: ../model/output)")
    parser.add_argument("--n_ens", type=int, default=20,
                        help="Number of ensemble members (default: 20)")
    parser.add_argument("--savefig", default="yates_validation.png",
                        help="Output figure filename")
    args = parser.parse_args()

    print(f"Loading ensemble output from: {args.outdir}")
    z, mass, radius, G, rho, n_loaded = load_ensemble(args.outdir, args.n_ens)
    print_summary(z, mass, radius, G, rho, n_loaded)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.suptitle(f"Yates 2017 IBM validation — {n_loaded} ensemble members, "
                 f"{len(mass):,} organisms", fontsize=12)

    plot_mass_distribution(mass, axes[0, 0])
    plot_radius_distribution(radius, axes[0, 1])
    plot_altitude_distribution(z, axes[1, 0])
    plot_G_distribution(G, axes[1, 1])

    plt.tight_layout()
    fig.savefig(args.savefig, dpi=150, bbox_inches="tight")
    print(f"Figure saved to: {args.savefig}")


if __name__ == "__main__":
    main()
