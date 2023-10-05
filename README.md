
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

| 🚧    | GitHub pages documentation will be available once kimmdy is public |
|---------------|:------------------------|

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
- install [act](https://github.com/nektos/act), easiest option is with github cli
    - install github cli (`gh`)
    - `gh extension install https://github.com/nektos/gh-act`
- run tests with `gh extension exec act -j test --artifact-server-path ./artifacts`
    - customize which python versions to test in `tox.ini` 
    - html coverage report is exported into `artifacts`
