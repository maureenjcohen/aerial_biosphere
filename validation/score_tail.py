#!/usr/bin/env python3
"""
score_tail.py
Score IBM runs against the Yates (2017) Figure 3 criteria, with emphasis on the
low-mass tail (the metric the grow-and-fission calibration is chasing).

For each run directory containing a member subdir with
ensemble_001_{state,pop}.dat, prints establishment, steady population, mass peak,
median radius, mean G, and tail metrics (10th-percentile mass, fraction of the
population below 1e-8 g and 1e-9 g).

Usage:
    python score_tail.py ../sweep_output/grow_*        # glob of run dirs
    python score_tail.py --label-from grow_ ../sweep_output/grow_*
"""
import argparse
import glob
import os
import sys

import numpy as np

YATES = dict(pop=50000, mass_peak_g=2e-8, radius_um=15.6)


def find_files(run_dir):
    """Return (state, pop) paths under run_dir, supporting member_NN/ or m/ layouts."""
    for sub in ("m", "member_01", "."):
        s = os.path.join(run_dir, sub, "ensemble_001_state.dat")
        p = os.path.join(run_dir, sub, "ensemble_001_pop.dat")
        if os.path.isfile(s):
            return s, p
    return None, None


def load_state(path):
    if os.path.getsize(path) == 0:
        return np.empty((0, 7))
    d = np.genfromtxt(path, comments="#", invalid_raise=False)
    if d.ndim == 1:
        d = d[np.newaxis, :]
    if d.size == 0 or d.shape[1] < 7:
        return np.empty((0, 7))
    good = (np.all(np.isfinite(d), axis=1)
            & (d[:, 1] > 0) & (d[:, 1] < 1e-8))   # finite, physical mass [kg]
    return d[good]


def steady_pop(pop_path, frac=0.5):
    """Mean population over the last `frac` of the run (steady-state estimate)."""
    p = np.genfromtxt(pop_path, comments="#", invalid_raise=False)
    if p.ndim == 1 or p.shape[0] < 2:
        return 0.0
    n = p.shape[0]
    return p[int(n * (1 - frac)):, 1].mean()


def score(run_dir, label):
    s_path, p_path = find_files(run_dir)
    if s_path is None:
        print(f"  {label:>8}: no state file (run incomplete?)")
        return
    state = load_state(s_path)
    if state.shape[0] == 0:
        print(f"  {label:>8}: EXTINCT / empty")
        return
    mass_g = state[:, 1] * 1e3
    rad_um = state[:, 2] * 1e6
    G = state[:, 3]
    pop = steady_pop(p_path)

    median_g = np.median(mass_g)
    p10_g = np.percentile(mass_g, 10)
    frac_8 = 100 * np.mean(mass_g < 1e-8)
    frac_9 = 100 * np.mean(mass_g < 1e-9)

    print(f"  {label:>8}: pop={pop:6.0f}  mass_med={median_g:.2e}g  "
          f"R_med={np.median(rad_um):5.2f}um  G={np.mean(G):.3f}  "
          f"| tail: p10={p10_g:.2e}g  <1e-8={frac_8:4.1f}%  <1e-9={frac_9:4.2f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+", help="run directories (globs ok)")
    ap.add_argument("--label-from", default=None,
                    help="strip this prefix from the basename for the label")
    args = ap.parse_args()

    dirs = []
    for pattern in args.run_dirs:
        dirs.extend(sorted(glob.glob(pattern)) or [pattern])

    print(f"Yates targets: pop~{YATES['pop']}, mass peak ~{YATES['mass_peak_g']:.0e} g, "
          f"R ~{YATES['radius_um']} um, substantial low-mass tail\n")
    for d in dirs:
        if not os.path.isdir(d):
            continue
        base = os.path.basename(d.rstrip("/"))
        label = base.replace(args.label_from, "") if args.label_from else base
        score(d, label)


if __name__ == "__main__":
    main()
