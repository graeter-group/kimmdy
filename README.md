# KIMMDY

[![](https://github.com/graeter-group/kimmdy/actions/workflows/tests.yml/badge.svg?branch=release-please--branches--main)](https://github.com/graeter-group/kimmdy/actions/workflows/tests.yml/?branch=release-please--branches--main)

Reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

## Installation

**Note**: KIMMDY requires [GROMACS](https://www.gromacs.org/) to be installed.
Some reactions need a GROMACS version patched with [PLUMED](https://www.plumed.org/).
The gromacs version name should then contain `MODIFIED` or `plumed`.

While it is possible to install KIMMDY with just `pip install kimmdy`,
this can take a while due to dependency resolution.
We recommend installing KIMMDY with [uv](https://docs.astral.sh/uv/) instead:

```bash
uv tool install -p 3.11 kimmdy
```

This installation includes only the most basic functionality as no plugins and analysis tools are installed.

To install the builtin reaction plugins, use

```bash
uv tool install -p 3.11 kimmdy[plugins]
```

To install the builtin example reactions and analysis tools use

```bash
uv tool install -p 3.11 kimmdy[reactions,analysis]
```

However, this is only half the fun!

To install KIMMDY with all currently available official plugins, like kimmdy-grappa, which properly parameterizes
your molecules for radicals using GrAPPa (Graph Attentional Protein
Parametrization) and reaction plugins like kimmdy-hat (for Hydrogen Atom Transfer) or kimmdy-hydrolysis use

```bash
uv tool install -p 3.11 kimmdy[plugins]
```

To uninstall KIMMDY again, use

```bash
uv tool uninstall kimmdy
```

## Documentation

The documentation is available [here](https://graeter-group.github.io/kimmdy/).

### Getting started

Head over to the [getting started](https://graeter-group.github.io/kimmdy/guide/tutorials/getting-started.html) tutorial.

## Development

### Development setup

Clone kimmdy and the default reaction and parameterization plugins and install requirements and kimmdy as editable via

```bash
git clone git@github.com:graeter-group/kimmdy.git --recurse-submodules
cd kimmdy
uv sync --extra plugins
```

Conventions:

* code style: black
* docstrings: numpy
* [Conventional commit](https://www.conventionalcommits.org/en/v1.0.0/) messages when possible for pretty release notes.

### Local testing

For developoment, we provide a docker image containing gromacs and multiple python versions to test against.  
To run the test locally, you must:

* install docker
* install [act](https://github.com/nektos/act), the easiest option is with github cli via `gh extension install https://github.com/nektos/gh-act`
* run tests with `gh extension exec act -j test --artifact-server-path ./artifacts`
* html coverage report is exported into `artifacts`


