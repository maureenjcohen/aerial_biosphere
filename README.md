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
│   ├── ensemble_launcher.py   Run a multi-member ensemble (one member per core)
│   ├── plot_distributions.py  Organism property distributions (cf. Yates Figs 3–5)
│   ├── plot_population.py      Population time series / steady-state check
│   └── sensitivity_runs.py     Parameter sensitivity experiments
│
└── *.pdf           Reference papers (Yates 2017, Morley 2014, Sagan & Salpeter 1976,
                    Seager 2021)
```

## Building

```
cd model
make
```

This produces the `bio_1d` executable.

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
member; the plotting scripts in `validation/` turn it into figures.

## Reference

Yates, J. S., Palmer, P. I., Biller, B., & Cockell, C. S. (2017).
*Atmospheric Habitable Zones in Y Dwarf Atmospheres.*
The Astrophysical Journal, 836, 184. https://doi.org/10.3847/1538-4357/836/2/184
