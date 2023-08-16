# KIMMDY

[![tests on latest release](https://github.com/hits-mbm-dev/kimmdy/actions/workflows/test-release.yml/badge.svg)](https://github.com/hits-mbm-dev/kimmdy/actions/workflows/test-release.yml)

Reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

## Quick start

* clone repository, e.g. `git clone https://github.com/hits-mbm-dev/kimmdy.git`
* `cd kimmdy`
* `python -m venv .venv`
* `source ./venv/bin/activate`
* `python -m pip install -e ./`
* Some reactions need a GROMACS version patched with PLUMED, gromacs name should then contain `MODIFIED` or `plumed`

## Full installation

* install conda and activate it
* `conda create -n "kimmdy" python=3.10`, i.e. create an environment for python 3.10 
* `conda activate kimmdy`
* `conda install -c conda-forge tensorflow==2.10` for the HAT plugin
* `conda install -c conda-forge openmm` for grappa parameterization
* `git clone https://github.com/hits-mbm-dev/HAT_reaction_plugin.git`
* `cd HAT_reaction_plugin/`
* `pip install -r requirements.txt`
* `cd ..`
* `git clone https://github.com/hits-mbm-dev/kimmdy.git`
* `cd kimmdy`
* `pip install -r requirements.txt`
* Some reactions need a GROMACS version patched with PLUMED, gromacs name should then contain `MODIFIED` or `plumed`

Other ways to install kimmdy are currently discouraged because of the high number of dependencies.


## Development setup

* `python -m pip install -r requirements.txt`
* code style: black
* docstrings: numpy
* [Conventional commit](https://www.conventionalcommits.org/en/v1.0.0/) messages when possible for pretty release notes.


## First simulation

* change directory to `example_triala`
* `ln -s ../../tests/test_files/assets/amber99sb-star-ildnp.ff ./amber99sb-star-ildnp.ff`
* run kimmdy: `kimmdy -l INFO`
* check output: `kimmdy.log`, `test_out_00X/`


## Local testing

For developoment, we provide a docker image containing gromacs and multiple python versions to test against.  
To run the test locally, you must:
- install docker
- install [act](https://github.com/nektos/act), easiest option is with github cli
    - install github cli (`gh`)
    - `gh extension install https://github.com/nektos/gh-act`
- run tests with `gh extension exec act -j test --artifact-server-path ./artifacts`
    - customize which python versions to test in `tox.ini` 
    - html coverage report is exported into `artifacts`
