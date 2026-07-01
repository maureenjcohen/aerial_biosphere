#!/usr/bin/env python3
"""Venus cloud biosphere — population oscillation + heritable-trait distributions.

Reads the per-snapshot histogram CSVs written by evolve_test when the driver
namelist has trait_hist=.true. (traits_r.csv, traits_X.csv, traits_summary.csv),
and renders a 3-panel figure:
  (1) population time series (ACTIVE / DORMANT / DEPOT) — the boom-bust cycle
  (2) cell-radius r distribution over time (column-normalised) + mean
  (3) reproduction half-life X distribution over time (column-normalised) + mean

Usage:  python3 plot_venus_traits.py [csv_dir] [out_png]
        csv_dir defaults to ../model, out_png to venus_trait_evolution.png (here).
"""
import os, sys, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
csv_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'model')
out_png = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'venus_trait_evolution.png')

# Trait bounds MUST match bio_venus.F90 (r_min..r_max, X_min_s..X_max_s) and NTBIN.
NB = 24
r_edges = np.exp(np.linspace(np.log(0.2), np.log(20.0), NB + 1))          # um
X_edges = np.exp(np.linspace(np.log(3600 / 86400.), np.log(100.0), NB + 1))  # days


def load_hist(path):
    raw = np.genfromtxt(path, delimiter=',', names=True)
    t = raw['t_day']
    H = np.vstack([raw[n_] for n_ in raw.dtype.names[1:]]).T   # (ntime, nbin)
    return t, H


def colnorm(H):
    tot = H.sum(axis=1, keepdims=True)
    frac = np.divide(H, tot, out=np.zeros_like(H), where=tot > 0)
    return np.ma.masked_where(frac <= 0, frac)


s = np.genfromtxt(os.path.join(csv_dir, 'traits_summary.csv'), delimiter=',', names=True)
tr, R = load_hist(os.path.join(csv_dir, 'traits_r.csv'))
tx, X = load_hist(os.path.join(csv_dir, 'traits_X.csv'))


def t_edges(t):
    return np.append(t, t[-1] + (t[-1] - t[-2]))


fig, ax = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

ax[0].plot(s['t_day'], s['nACT'], lw=1.2, label='ACTIVE (in colonies)')
ax[0].plot(s['t_day'], s['nDOR'], lw=1.0, label='DORMANT (free spores)')
ax[0].plot(s['t_day'], s['nDEP'], lw=1.0, label='DEPOT (haze bank)')
ax[0].set_ylabel('population'); ax[0].set_yscale('log')
ax[0].legend(loc='upper right', fontsize=8)
ax[0].set_title('Population time series (boom-bust oscillation)')

pm = ax[1].pcolormesh(t_edges(tr), r_edges, colnorm(R).T, cmap='viridis', shading='flat')
ax[1].plot(s['t_day'], s['r_mean_um'], color='w', lw=1.3, label='mean')
ax[1].axhline(0.2, color='r', ls=':', lw=1, label='r_min (Seager floor)')
ax[1].set_yscale('log'); ax[1].set_ylabel('cell radius r  [um]')
ax[1].set_title('Cell-size distribution (column-normalised fraction of ACTIVE pop)')
ax[1].legend(loc='upper right', fontsize=8); fig.colorbar(pm, ax=ax[1], label='fraction')

pm = ax[2].pcolormesh(t_edges(tx), X_edges, colnorm(X).T, cmap='magma', shading='flat')
ax[2].plot(s['t_day'], s['X_mean_d'], color='c', lw=1.3, label='mean')
ax[2].set_yscale('log'); ax[2].set_ylabel('repro half-life X  [days]')
ax[2].set_xlabel('time [Earth-days]')
ax[2].set_title('Reproduction-rate-gene distribution (column-normalised)')
ax[2].legend(loc='upper right', fontsize=8); fig.colorbar(pm, ax=ax[2], label='fraction')

fig.tight_layout()
fig.savefig(out_png, dpi=110)
print('wrote', out_png)
