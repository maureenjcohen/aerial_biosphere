# Aerial Biosphere

A 1-D individual-based model (IBM) of an aerial microbial ecosystem living in a
planetary atmosphere, following Yates et al. (2017). Organisms are modelled as
buoyant spherical shells that grow by consuming biomass, are carried by
convection and gravitational settling, reproduce, and die — allowing a population
to evolve a survival strategy suited to its atmospheric environment.

The reference application is the Atmospheric Habitable Zone (AHZ) of a cool Y
dwarf, but the framework is intended to be adapted to other atmospheres.

## Repository layout

```
aerial_biosphere/
├── model/          Standalone Fortran IBM
│   ├── bio_params.F90       Physical constants and default parameters (SI units)
│   ├── bio_atmosphere.F90   Atmosphere profiles (T, density, viscosity) over the AHZ
│   ├── bio_model.F90        Organism lifecycle: growth, transport, reproduction, death
│   ├── bio_driver.F90       Ensemble runner, namelist I/O, output files
│   ├── bio_run.nml          Control-run configuration (editable)
│   └── Makefile             Build the `bio_1d` executable
│
├── validation/     Python tools to run ensembles and compare against Yates (2017)
│   ├── ensemble_launcher.py        Run a multi-member ensemble (one member per core)
│   ├── reproduce_yates_fig3.ipynb  Notebook reproducing Yates (2017) Figure 3
│   ├── plot_distributions.py       Organism property distributions (cf. Yates Figs 3–5)
│   ├── plot_population.py           Population time series / steady-state check
│   └── sensitivity_runs.py          Parameter sensitivity experiments
│
└── *.pdf           Reference papers (Yates 2017, Morley 2014, Sagan & Salpeter 1976,
                    Seager 2021)
```

## Building

```
cd model
make        # debug build (-O2, bounds checking) — good for development
make fast   # optimized build (-O3 -march=native) — use for production ensembles
```

Either produces the `bio_1d` executable. Use `make clean` before rebuilding on a
different machine or compiler.

## Running

A single run reads its configuration from a namelist:

```
cd model
./bio_1d bio_run.nml
```

To run a full ensemble in parallel, use the launcher:

```
cd validation
python3 ensemble_launcher.py
```

Output (organism state snapshots and population time series) is written per
member; the plotting scripts and `reproduce_yates_fig3.ipynb` in `validation/`
turn it into figures.

The control run reproduces the Yates (2017) Figure 3 distributions: a steady
population of ~50,000 organisms with a mass distribution peaking near 2×10⁻⁸ g
(radius ~15.6 μm) and growth strategy *G* evolving toward the solid end.

## Differences from Yates (2017)

This model follows Yates closely but changes two things — the initial conditions
and how biomass is supplied — both controlled from `model/bio_run.nml`.

### Initial conditions

Like Yates, the population starts cold: 100 founder organisms at an initial mass
`m_init = 1×10⁻¹² kg` (Yates' "approximate mass of 10⁻⁹ g").

Organisms reproduce at a fixed reproduction mass `m_repr`. If every founder is
pinned at `m_init`, none can grow to the neutral-buoyancy mass before being swept
out of the AHZ, and the population washes out before it can establish. To let the
population establish and self-organize to the Yates steady state, the founders'
reproduction mass is seeded from a log-uniform spread up to `mrepr_seed_max`
(default 3×10⁻¹¹ kg).

### Biomass supply and boundaries

Yates uses a fixed, conserved biomass pool. Here the food supply is an open
conveyor instead:

- **Source.** Biomass enters at the bottom of the AHZ at a tunable rate
  `biomass_flux` [kg/day], is advected upward by convection, and is consumed by
  organisms along the way.
- **Sinks.** Both AHZ boundaries are absorbing: organisms that leave through the
  top or bottom are removed and their biomass is lost. The top in particular is a
  pure sink — microbes carried aloft by updrafts have no mechanism to rain back
  down.

This makes `biomass_flux` the control parameter that sets the sustainable
population size. The population scales approximately linearly with the flux: a
flux of ~2×10⁻¹¹ kg/day is the minimum that sustains a population, while
~2.5×10⁻⁸ kg/day reproduces the Yates control carrying capacity (~50,000
organisms). The per-organism mass distribution is set by the atmospheric physics
and is independent of the flux — the flux changes *how many* organisms persist,
not *how large* they are.

## Reference

Yates, J. S., Palmer, P. I., Biller, B., & Cockell, C. S. (2017).
*Atmospheric Habitable Zones in Y Dwarf Atmospheres.*
The Astrophysical Journal, 836, 184. https://doi.org/10.3847/1538-4357/836/2/184
