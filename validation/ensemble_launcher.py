#!/usr/bin/env python3
"""
ensemble_launcher.py
Run the Yates IBM ensemble with one member per processor.

Usage
-----
  # Run all 20 members in parallel (one per core):
  python ensemble_launcher.py

  # Limit to N simultaneous jobs (e.g. if cluster has fewer free cores):
  python ensemble_launcher.py --jobs 8

  # Run a single member (useful for testing):
  python ensemble_launcher.py --member 3

  # Override key parameters at the command line:
  python ensemble_launcher.py --n_members 20 --t_sim_years 100 --b_factor 3

Directory layout produced
-------------------------
  ensemble_output/
    namelists/   member_01.nml  member_02.nml  ...
    member_01/   ensemble_001_pop.dat   ensemble_001_state.dat   run.log
    member_02/   ...
    run_summary.txt
"""

import argparse
import os
import subprocess
import sys
import textwrap
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration — edit these to match your local paths / run settings
# ---------------------------------------------------------------------------

# Path to compiled model binary, relative to this script
EXE_REL = "../model/bio_1d"

# Root output directory, relative to this script
OUTDIR_REL = "../ensemble_output"

# Default ensemble settings
DEFAULT_N_MEMBERS   = 20
DEFAULT_SEED_BASE   = 42      # seed for member k = SEED_BASE + (k-1)*1000

# Default namelist parameters (Yates 2017 Table 1 control run)
DEFAULT_PARAMS = dict(
    n_levels        = 106,
    n_init_orgs     = 100,
    m_init          = 1.857e-11,    # warm start near neutral buoyancy [kg]
    halflife        = 30.0,         # organism half-life [days]
    growth_rate_day = 0.0,          # optional max growth-rate cap [/day]; 0 = biomass-limited (Yates)
    b_total_kg      = 1.0e-6,       # conserved biomass pool [kg]  (B=1 control)
    b_factor        = 1.0,
    v_conv          = 1.0,          # convective velocity [m/s]
    dt_hrs          = 6.0,
    t_sim_years     = 100.0,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent


def resolve_exe():
    p = (HERE / EXE_REL).resolve()
    if not p.is_file():
        sys.exit(
            f"ERROR: model executable not found at {p}\n"
            f"Compile first (from aerial_biosphere/model/):\n"
            f"  gfortran -O3 -march=native -c bio_params.F90 bio_atmosphere.F90 "
            f"bio_model.F90 bio_driver.F90 && "
            f"gfortran -O3 -march=native -o bio_1d bio_params.o bio_atmosphere.o "
            f"bio_model.o bio_driver.o"
        )
    return p


def outdir_root():
    return (HERE / OUTDIR_REL).resolve()


def member_dir(k: int) -> Path:
    return outdir_root() / f"member_{k:02d}"


def nml_file(k: int) -> Path:
    return outdir_root() / "namelists" / f"member_{k:02d}.nml"


def write_namelist(k: int, params: dict, seed_base: int) -> Path:
    """Write a 1-member namelist for member k (1-indexed)."""
    member_dir(k).mkdir(parents=True, exist_ok=True)
    nml_dir = outdir_root() / "namelists"
    nml_dir.mkdir(parents=True, exist_ok=True)

    # seed_base chosen so parallel member k reproduces the same trajectory as
    # the sequential driver's member k:
    #   sequential: seed = seed_base + iens*1000  (iens=k)
    #   parallel:   seed = (seed_base + (k-1)*1000) + 1*1000  — same result
    seed = seed_base + (k - 1) * 1000
    p = params

    nml_text = textwrap.dedent(f"""\
        &bio_run
          n_ensemble      = 1
          seed_base       = {seed}
          n_levels        = {p['n_levels']}
          n_init_orgs     = {p['n_init_orgs']}
          m_init          = {p['m_init']:.3e}
          halflife        = {p['halflife']}
          growth_rate_day = {p['growth_rate_day']}
          b_total_kg      = {p['b_total_kg']:.3e}
          b_factor        = {p['b_factor']}
          v_conv          = {p['v_conv']}
          dt_hrs          = {p['dt_hrs']}
          t_sim_years     = {p['t_sim_years']}
          outdir          = '{member_dir(k)}'
        /
    """)

    path = nml_file(k)
    path.write_text(nml_text)
    return path


def run_member(k: int):
    """Run member k; returns (k, exit_code)."""
    exe = resolve_exe()
    nml = nml_file(k)
    log = member_dir(k) / "run.log"

    with open(log, "w") as f:
        ret = subprocess.run([str(exe), str(nml)], stdout=f, stderr=subprocess.STDOUT)

    return k, ret.returncode


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run the Yates IBM ensemble in parallel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--member", type=int, metavar="N",
                        help="Run only member N (1-indexed); skip the rest.")
    parser.add_argument("--jobs", type=int, default=DEFAULT_N_MEMBERS,
                        help=f"Max simultaneous members (default: {DEFAULT_N_MEMBERS}).")
    parser.add_argument("--n_members", type=int, default=DEFAULT_N_MEMBERS)
    parser.add_argument("--t_sim_years", type=float, default=DEFAULT_PARAMS["t_sim_years"])
    parser.add_argument("--b_factor",    type=float, default=DEFAULT_PARAMS["b_factor"])
    parser.add_argument("--v_conv",      type=float, default=DEFAULT_PARAMS["v_conv"])
    parser.add_argument("--seed_base",   type=int,   default=DEFAULT_SEED_BASE)
    args = parser.parse_args()

    params = dict(DEFAULT_PARAMS)
    params["t_sim_years"] = args.t_sim_years
    params["b_factor"]    = args.b_factor
    params["v_conv"]      = args.v_conv

    n = args.n_members
    members = [args.member] if args.member else list(range(1, n + 1))

    # Write all namelists up front so workers don't race on directory creation
    for k in members:
        write_namelist(k, params, args.seed_base)

    print(f"Output directory : {outdir_root()}")
    print(f"Members to run   : {members[0]}–{members[-1]}  ({len(members)} total)")
    print(f"Parallel jobs    : {min(args.jobs, len(members))}")
    print()

    results = {}
    with ProcessPoolExecutor(max_workers=min(args.jobs, len(members))) as pool:
        futures = {pool.submit(run_member, k): k for k in members}
        for fut in as_completed(futures):
            k, rc = fut.result()
            status = "OK" if rc == 0 else f"FAILED (exit {rc})"
            print(f"  member {k:02d}: {status}", flush=True)
            results[k] = rc

    # Write a summary
    summary = outdir_root() / "run_summary.txt"
    with open(summary, "w") as f:
        f.write("# member  exit_code  status\n")
        for k in sorted(results):
            rc = results[k]
            f.write(f"  {k:6d}  {rc:9d}  {'ok' if rc == 0 else 'FAILED'}\n")

    failed = [k for k, rc in results.items() if rc != 0]
    print()
    if failed:
        print(f"FAILED members: {failed}")
        sys.exit(1)
    print(f"All {len(members)} members completed successfully.")
    print(f"Results in: {outdir_root()}")


if __name__ == "__main__":
    main()
