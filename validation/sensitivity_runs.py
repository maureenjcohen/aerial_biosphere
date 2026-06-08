"""
sensitivity_runs.py
Automate the Yates Table 1 sensitivity runs by writing namelist files
and launching the model for each parameter set.

Yates Table 1 (control run + variations):
  Control:  v_conv=1.0 m/s, B=1  (b_factor=1)
  B runs:   b_factor = 1, 3, 6    (v_conv fixed at 1.0)
  v runs:   v_conv = 0.01, 1.0, 10.0 m/s  (b_factor=1)
            m_init must scale with analytical equilibrium mass for each v_conv

m_init scaling: m_eq ∝ v_conv (from Stokes, R_eq ∝ sqrt(v), m_eq ∝ R_eq^3 ∝ v^(3/2))
  v_conv=0.01 m/s: m_init ~ 1.857e-11 × (0.01/1.0)^(3/2) = 1.857e-11 × 0.001 = 1.857e-14 kg
  v_conv=1.0  m/s: m_init = 1.857e-11 kg  (control)
  v_conv=10.0 m/s: m_init ~ 1.857e-11 × (10/1.0)^(3/2) = 1.857e-11 × 31.62 = 5.87e-10 kg

Usage:
    python sensitivity_runs.py --exe ../model/bio_1d --outbase /tmp/sensitivity
"""

import argparse
import os
import subprocess
import sys

BASE_NML = """
&bio_run
  n_ensemble      = 20
  seed_base       = 42
  n_levels        = 106
  n_init_orgs     = 100
  halflife        = 30.0
  growth_rate_day = 0.70
  dt_hrs          = 6.0
  t_sim_years     = 100.0
  m_init          = {m_init}
  v_conv          = {v_conv}
  b_total_kg      = 1.0e-6
  b_factor        = {b_factor}
  outdir          = '{outdir}'
/
"""

# Warm-start mass for control v_conv=1 m/s
M_INIT_CONTROL = 1.857e-11

def m_init_for_vconv(v):
    """Scale warm-start mass with v_conv^(3/2) (Stokes equilibrium scaling)."""
    return M_INIT_CONTROL * (v / 1.0) ** 1.5


RUNS = [
    # name           v_conv   b_factor
    ("control",      1.0,     1.0),
    ("B3",           1.0,     3.0),
    ("B6",           1.0,     6.0),
    ("v001",         0.01,    1.0),
    ("v10",          10.0,    1.0),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe",     default="../model/bio_1d")
    parser.add_argument("--outbase", default="/tmp/sensitivity")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print nmls without running")
    args = parser.parse_args()

    exe = os.path.abspath(args.exe)
    if not os.path.isfile(exe):
        sys.exit(f"Executable not found: {exe}")

    os.makedirs(args.outbase, exist_ok=True)

    for name, v_conv, b_factor in RUNS:
        rundir = os.path.join(args.outbase, name)
        os.makedirs(rundir, exist_ok=True)
        m_init = m_init_for_vconv(v_conv)
        nml = BASE_NML.format(
            m_init=f"{m_init:.3e}",
            v_conv=v_conv,
            b_factor=b_factor,
            outdir=rundir,
        )
        nml_path = os.path.join(rundir, "bio_run.nml")
        with open(nml_path, "w") as f:
            f.write(nml)
        print(f"[{name}] v_conv={v_conv} m/s, B={b_factor}, m_init={m_init:.2e} kg")
        print(f"        outdir: {rundir}")
        if not args.dry_run:
            print(f"        running {exe} ...")
            ret = subprocess.run([exe, nml_path], cwd=rundir, capture_output=False)
            if ret.returncode != 0:
                print(f"        *** FAILED (return code {ret.returncode}) ***")
            else:
                print(f"        done.")
        print()


if __name__ == "__main__":
    main()
