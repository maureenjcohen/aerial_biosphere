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

The control run reproduces the Yates (2017) Figure 3 *distributions*: a steady
population whose mass distribution peaks near 2×10⁻⁸ g (radius ~15.6 μm), with
growth strategy *G* evolving toward the solid end and organisms concentrated in
the upper AHZ. Validation targets the per-organism distributions and steady-state
stationarity — not an absolute population count, which Yates do not report and
which here is simply set by the chosen food supply.

## Differences from Yates (2017)

This model follows Yates closely but makes a few implementation choices the paper
leaves open. All are controlled from `model/bio_run.nml`.

### Initial conditions — genuine cold start

The population starts cold, as in Yates: 100 founder organisms at an initial mass
`m_init = 1×10⁻¹² kg` (Yates' "approximate mass of 10⁻⁹ g"), with random growth
strategy, density, and altitude. Each founder's reproduction mass equals its birth
mass, so the population must grow and self-organize to a stable mass/size strategy
on its own. A high maximum specific growth rate (`growth_rate_day`, the μ_max
ceiling) lets founders grow toward the neutral-buoyancy mass before convection
sweeps them out of the AHZ, so the population establishes from the cold start with
**no founder seeding**.

### Biomass supply and boundaries

Yates uses a fixed, conserved biomass pool. Here the food supply is an open
conveyor instead:

- **Source.** Biomass enters at the bottom of the AHZ at a tunable rate
  `biomass_flux` [kg/day] and is advected upward by convection.
- **Recycling.** When an organism dies of old age inside the domain, its mass is
  returned to the biomass of the layer where it died.
- **Sinks.** Both AHZ boundaries are absorbing: organisms (and any biomass) that
  leave through the top or bottom are removed. The top in particular is a pure
  sink — microbes carried aloft by updrafts have no mechanism to rain back down.

`biomass_flux` therefore sets the carrying capacity: more food supports more
organisms. The absolute population is a free parameter of the model — Yates report
no population count, only that the total is steady — so the per-organism mass
distribution (which the atmospheric physics fixes) is what we validate, not the
headcount.

### Biomass consumption and growth

Each timestep an organism grows by consuming biomass from its atmospheric layer.
The layer's available biomass `B_layer` is **shared** among the organisms there in
proportion to `mass^p`, where `p = uptake_exp`:

```
dBᵢ = ( B_layer · mᵢ^p / Σⱼ mⱼ^p ),   capped at   μ_max · Δt · mᵢ
```

Consumed biomass becomes organism mass one-to-one (`mᵢ → mᵢ + dBᵢ`), and the
layer's biomass is debited so it cannot go negative.

- With `uptake_exp = 1` (default) the share is proportional to mass, so every
  organism grows at the same specific rate `dBᵢ/mᵢ` regardless of size — consistent
  with Yates, who redistribute biomass "as a function of organism weight."
- With `uptake_exp < 1` (e.g. 2/3, surface-area-limited uptake) smaller organisms
  receive a higher specific growth rate. Yates did not include this; it is a
  physically-motivated lever kept for future experiments (e.g. the Venus AHZ).

The maximum specific growth rate `growth_rate_day` (μ_max) caps how fast an
organism can grow in one step. It only binds transiently — when food is locally
abundant, such as during cold-start establishment — because in a food-limited
steady state the shared portion is the smaller term and actual growth runs well
below the cap.

(An experimental alternative consumption model, `eat_mode = 1`, lets organisms
deplete the layer biomass at μ_max in turn rather than share it; the default
`eat_mode = 0` is the shared model described above.)

## Reference

Yates, J. S., Palmer, P. I., Biller, B., & Cockell, C. S. (2017).
*Atmospheric Habitable Zones in Y Dwarf Atmospheres.*
The Astrophysical Journal, 836, 184. https://doi.org/10.3847/1538-4357/836/2/184
