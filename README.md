# KIMMDY

[![tests on latest release](https://github.com/hits-mbm-dev/kimmdy/actions/workflows/test-release.yml/badge.svg)](https://github.com/hits-mbm-dev/kimmdy/actions/workflows/test-release.yml)

Reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

## Installation

**Note**: Some reactions need a GROMACS version patched with PLUMED, gromacs name should then contain `MODIFIED` or `plumed`.

### Bare Installation
```bash
git clone https://github.com/hits-mbm-dev/kimmdy.git
cd kimmdy
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ./
```

This installation includes only the most basic functionality as no plugins and analysis tool dependencies are installed. Plugins can be installed using `python -m pip install -e ./` in their module directories (e.g. `kimmdy/plugins/default_reactions`). The analysis tool dependencies can be installed with `python -m pip install -e ./[analysis]` in `kimmdy/`.

### Full installation
```bash
conda create -n kimmdy_full python=3.10 tensorflow==2.10 openmm
conda activate kimmdy_full
git clone https://github.com/hits-mbm-dev/HAT_reaction_plugin.git
cd HAT_reaction_plugin/
pip install -r requirements.txt
cd ..
git clone https://github.com/hits-mbm-dev/grappa.git
cd grappa
pip install -e .
cd ..
git clone https://github.com/hits-mbm-dev/kimmdy.git
cd kimmdy
pip install -r requirements.txt
pip install -e ./[parameterization_plugins]
```

Other ways to install kimmdy with all plugins are currently discouraged because of the high number of dependencies.

## Documentation

### GitHub pages

[GitHub Pages documentation](https://hits-mbm-dev.github.io/kimmdy/)

### Local documentation
```bash
cd <kimmdy root directory>
pip install -r requirements.txt # Documentation depends on packages in requirements.txt
make docs
quarto preview # Should open the documentation in your standard browser
```

## Development setup

* install packages necessary for developmentmake `python -m pip install -r requirements.txt`
* code style: black
* docstrings: numpy
* [Conventional commit](https://www.conventionalcommits.org/en/v1.0.0/) messages when possible for pretty release notes.


## First simulation

* change directory to `example_triala`
* run kimmdy: `kimmdy -l INFO`
* check output: `kimmdy.log`, `test_out_00X/`


## Local testing

For developoment, we provide a docker image containing gromacs and multiple python versions to test against.  
To run the test locally, you must:

- install docker
- install [act](https://github.com/nektos/act), the easiest option is with github cli via `gh extension install https://github.com/nektos/gh-act`
- run tests with `gh extension exec act -j test --artifact-server-path ./artifacts`
- html coverage report is exported into `artifacts`

