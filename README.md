# KIMMDY

[![test latest release](https://github.com/graeter-group/kimmdy/actions/workflows/tests.yml/badge.svg?branch=release-please--branches--main)](https://github.com/graeter-group/kimmdy/actions/workflows/tests.yml/?branch=release-please--branches--main)

Reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

## Installation

**Note**: KIMMDY requires [GROMACS](https://www.gromacs.org/) to be installed.
Some reactions need a GROMACS version patched with [PLUMED](https://www.plumed.org/).
The gromacs version name should then contain `MODIFIED` or `plumed`.

```bash
pip install kimmdy
```

This installation includes only the most basic functionality as no plugins and analysis tools are installed.

To install the builtin reaction plugins, use

```bash
pip install kimmdy[reactions]
```

To install the builtin reactions and analysis tools use

```bash
pip install kimmdy[reactions,analysis]
```

However, this is only half the fun!

KIMMDY has two exciting plugins in the making, which properly parameterize your molecules
for radicals using GrAPPa (Graph Attentional Protein Parametrization) and predict
Hydrogen Atom Transfer (HAT) rates.

Full installation instructions are available [here](https://grater-group.github.io/kimmdy/guide/how-to/install-ml-plugins.html)

## Documentation

The documentation is available [here](https://grater-group.github.io/kimmdy/).

### Getting started

Head over to the [getting started](https://graeter-group.github.io/kimmdy/guide/tutorials/getting-started.html) tutorial.

## Development

### Development setup

Clone kimmdy and the default reaction and parameterization plugins and install requirements and kimmdy as editable via

```bash
git clone git@github.com:graeter-group/kimmdy.git
git clone git@github.com:graeter-group/kimmdy-reactions.git
git clone git@github.com:graeter-group/kimmdy-grappa.git
cd kimmdy
python -m venv .venv
source ./venv/bin/activate
pip install -r requirements.txt
```

Conventions:

* code style: black
* docstrings: numpy
* [Conventional commit](https://www.conventionalcommits.org/en/v1.0.0/) messages when possible for pretty release notes.

### Local testing

For developoment, we provide a docker image containing gromacs and multiple python versions to test against.  
To run the test locally, you must:

- install docker
- install [act](https://github.com/nektos/act), the easiest option is with github cli via `gh extension install https://github.com/nektos/gh-act`
- run tests with `gh extension exec act -j test --artifact-server-path ./artifacts`
- html coverage report is exported into `artifacts`

