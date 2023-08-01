# KIMMDY

Reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

## Quick start

* clone repository, e.g. `git clone https://github.com/hits-mbm-dev/kimmdy.git`
* `cd kimmdy`
* `python -m venv .venv`
* `source ./venv/bin/activate`
* `python -m pip install -r requirments.txt`
* Some rections need a GROMACS version patched with PLUMED

## First simulation

* change directory to `example_triala`
* `ln -s ../../tests/test_files/assets/amber99sb-star-ildnp.ff ./amber99sb-star-ildnp.ff`
* run kimmdy: `kimmdy -l INFO`
* check output: `kimmdy.log`, `test_out_00X/`


