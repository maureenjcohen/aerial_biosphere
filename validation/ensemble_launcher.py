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

# Control namelist, relative to this script.  Run values are read from here;
# DEFAULT_PARAMS below are only fallbacks for keys absent from the namelist.
NML_REL = "../model/bio_run.nml"

# Root output directory, relative to this script (override with --outdir)
OUTDIR_REL = "../ensemble_output"
_OUTDIR_OVERRIDE = None   # absolute path set from --outdir; see outdir_root()

# Default ensemble settings (used only if absent from the namelist)
DEFAULT_N_MEMBERS   = 20
DEFAULT_SEED_BASE   = 42      # seed for member k = SEED_BASE + (k-1)*1000

# Fallback namelist parameters (Yates 2017 Table 1 control run).  Any key also
# present in bio_run.nml is overridden by the namelist value; CLI flags override
# both.
DEFAULT_PARAMS = dict(
    n_levels        = 106,
    n_init_orgs     = 100,
    m_init          = 1.0e-12,      # cold start (Yates Table 1 control: 1e-9 g) [kg]
    halflife        = 30.0,         # organism half-life [days]
    growth_rate_day = 2.5,          # max specific growth rate (NPZ mu_max) [/day]
    mrepr_seed_max  = 0.0,          # founder m_repr spread [kg]; <=m_init = no spread
    b_total_kg      = 1.0e-6,       # conserved biomass pool [kg]  (B=1 control)
    b_factor        = 1.0,
    v_conv          = 1.0,          # convective velocity [m/s]
    dt_hrs          = 6.0,
    t_sim_years     = 100.0,
)

# Keys read from the namelist as ensemble-level controls (not per-member params)
ENSEMBLE_KEYS = {"n_ensemble", "seed_base", "outdir"}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent


def _coerce(val: str):
    """Convert a namelist scalar string to int / float / bool / str."""
    if (val.startswith("'") and val.endswith("'")) or \
       (val.startswith('"') and val.endswith('"')):
        return val[1:-1]
    low = val.lower()
    if low in (".true.", "t", ".t."):
        return True
    if low in (".false.", "f", ".f."):
        return False
    try:
        return int(val)
    except ValueError:
        pass
    try:
        # Fortran allows 'd'/'D' exponent markers (e.g. 1.0d-6)
        return float(val.replace("d", "e").replace("D", "e"))
    except ValueError:
        return val


def parse_namelist(path: Path) -> dict:
    """Parse a simple Fortran namelist (&group ... /) into a typed dict.

    Reads ``key = value`` lines from the first namelist group, stripping ``!``
    comments and inferring types.  Only scalar values are handled, which is all
    bio_run.nml contains.  Returns {} if the file is missing.
    """
    params: dict = {}
    if not path.is_file():
        return params
    in_group = False
    for raw in path.read_text().splitlines():
        line = raw.split("!", 1)[0].strip()      # drop trailing comment
        if not line:
            continue
        if line.startswith("&"):
            in_group = True
            continue
        if line.startswith("/"):
            break                                 # end of group
        if not in_group or "=" not in line:
            continue
        key, val = line.split("=", 1)
        params[key.strip().lower()] = _coerce(val.strip().rstrip(",").strip())
    return params


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
    if _OUTDIR_OVERRIDE is not None:
        return Path(_OUTDIR_OVERRIDE).resolve()
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
          mrepr_seed_max  = {p['mrepr_seed_max']:.3e}
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
    global _OUTDIR_OVERRIDE

    parser = argparse.ArgumentParser(
        description="Run the Yates IBM ensemble in parallel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # Precedence for every run value: DEFAULT_PARAMS < bio_run.nml < CLI flag.
    # Param-overriding flags default to None so we can tell when they were given.
    parser.add_argument("--nml", default=NML_REL,
                        help=f"Control namelist to read (default: {NML_REL}).")
    parser.add_argument("--member", type=int, metavar="N",
                        help="Run only member N (1-indexed); skip the rest.")
    parser.add_argument("--jobs", type=int, default=None,
                        help="Max simultaneous members (default: all of them).")
    parser.add_argument("--n_members",   type=int,   default=None)
    parser.add_argument("--t_sim_years", type=float, default=None)
    parser.add_argument("--b_factor",    type=float, default=None)
    parser.add_argument("--v_conv",      type=float, default=None)
    parser.add_argument("--growth_rate_day", type=float, default=None,
                        help="Max specific growth rate (NPZ mu_max) [/day].")
    parser.add_argument("--m_init",      type=float, default=None,
                        help="Initial organism mass [kg].")
    parser.add_argument("--mrepr_seed_max", type=float, default=None,
                        help="Founder m_repr log-uniform spread up to this mass [kg]; "
                             "<= m_init disables (all founders at m_init).")
    parser.add_argument("--seed_base",   type=int,   default=None)
    parser.add_argument("--outdir",      default=None,
                        help="Output root (default: ../ensemble_output).")
    args = parser.parse_args()

    # 1. fallbacks, 2. namelist, 3. CLI
    nml_path = Path(args.nml)
    if not nml_path.is_absolute():
        nml_path = (HERE / nml_path).resolve()
    nml = parse_namelist(nml_path)

    params = dict(DEFAULT_PARAMS)
    for key in params:
        if key in nml:
            params[key] = nml[key]
    for key in ("t_sim_years", "b_factor", "v_conv", "growth_rate_day",
                "m_init", "mrepr_seed_max"):
        if getattr(args, key) is not None:
            params[key] = getattr(args, key)

    n_members = args.n_members if args.n_members is not None \
        else nml.get("n_ensemble", DEFAULT_N_MEMBERS)
    seed_base = args.seed_base if args.seed_base is not None \
        else nml.get("seed_base", DEFAULT_SEED_BASE)
    if args.outdir is not None:
        _OUTDIR_OVERRIDE = args.outdir

    members = [args.member] if args.member else list(range(1, n_members + 1))
    jobs = args.jobs if args.jobs is not None else len(members)

    # Write all namelists up front so workers don't race on directory creation
    for k in members:
        write_namelist(k, params, seed_base)

    print(f"Namelist read    : {nml_path}{'' if nml else '  (missing — using defaults)'}")
    print(f"Output directory : {outdir_root()}")
    print(f"m_init / growth  : {params['m_init']:.3e} kg  /  {params['growth_rate_day']} /day")
    print(f"Members to run   : {members[0]}–{members[-1]}  ({len(members)} total)")
    print(f"Parallel jobs    : {min(jobs, len(members))}")
    print()

    results = {}
    with ProcessPoolExecutor(max_workers=min(jobs, len(members))) as pool:
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
