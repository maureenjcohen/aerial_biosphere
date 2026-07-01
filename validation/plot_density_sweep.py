#!/usr/bin/env python3
"""Density sweep: how dense must in-droplet cells be for colonies to sediment out
and close the Seager return leg (reseed the sub-cloud haze depot)?

Reads density_sweep.csv (columns: config, rho_cell, n_rainout, ..., cells_depot,
rmin_settle_um, final_nDEP) and plots, per cloud config:
  (top)  cells delivered to the depot over the run   vs cell bulk density
  (bot)  final depot standing population (loop closed?) vs cell bulk density
Reference densities (water, hydrated cell, H2SO4, silica, magnetite, iron, osmium)
are marked to place the required density in physical context.
"""
import sys, os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else 'density_sweep.csv'
out = sys.argv[2] if len(sys.argv) > 2 else 'venus_density_sweep.png'

rows = {}
with open(src) as f:
    for r in csv.DictReader(f):
        rows.setdefault(r['config'], []).append(r)

REFS = [('water', 1000), ('hydrated cell', 1100), ('H2SO4 drop', 1800),
        ('silica', 2200), ('magnetite', 5200), ('iron', 7870), ('osmium', 22600)]

fig, ax = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
colors = {'VPCM12': 'tab:blue', 'VPCM12+Haus3': 'tab:red'}
for cfg, rr in rows.items():
    rho = np.array([float(x['rho_cell']) for x in rr])
    cells = np.array([float(x['cells_depot'] or 0) for x in rr])
    ndep = np.array([float(x['final_nDEP'] or 0) for x in rr])
    c = colors.get(cfg, 'k')
    ax[0].plot(rho, np.maximum(cells, 0.5), 'o-', color=c, label=cfg)
    ax[1].plot(rho, np.maximum(ndep, 0.5), 'o-', color=c, label=cfg)

for a in ax:
    a.set_xscale('log'); a.set_yscale('log')
    for name, d in REFS:
        a.axvline(d, color='gray', ls=':', lw=0.8)
    a.grid(True, which='both', alpha=0.15)
# label refs along the top axis
ymax = ax[0].get_ylim()[1]
for name, d in REFS:
    ax[0].text(d, ymax, name, rotation=90, va='top', ha='right', fontsize=7, color='gray')

ax[1].axhline(5000, color='green', ls='--', lw=1, label='seed bank (5000)')
ax[0].set_ylabel('cells delivered to depot (cumulative)')
ax[0].set_title('Seager return-leg closure vs in-droplet cell bulk density')
ax[0].legend(fontsize=8)
ax[1].set_ylabel('final depot population')
ax[1].set_xlabel(r'cell bulk density $\rho_{cell}$  [kg m$^{-3}$]')
ax[1].legend(fontsize=8, loc='upper left')

fig.tight_layout()
here = os.path.dirname(os.path.abspath(__file__))
path = out if os.path.isabs(out) else os.path.join(here, out)
fig.savefig(path, dpi=120)
print('wrote', path)
