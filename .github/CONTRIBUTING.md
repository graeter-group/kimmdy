# Contributing to KIMMDY

We welcome contributions to KIMMDY!

You can contribute in many ways:

- By opening issues to provide feedback and share ideas.
- By fixing typos in documentation
- By submitting Pull Request (PR) to fix opened issues
- By submitting Pull Request (PR) to suggest new features (it is considered good practice to open an issue for discussion before working on a pull request for a new feature).

## To submit a contribution using a Pull Request

1.  [Fork](https://github.com/hits-mbm-dev/kimmdy/fork) the repository, clone it locally, and make your changes in a new branch specific to the PR. For example:

    ```bash
    # clone your fork
    git clone git@github.com:<username>/kimmdy

    # checkout a new branch
    git switch -c fix/myfix
    ```

3.  Submit the [pull request](https://help.github.com/articles/using-pull-requests). It is ok to submit as a draft if you are still working on it but would like some feedback from us. It is always good to share in the open that you are working on it.

We'll try to be as responsive as possible in reviewing and accepting pull requests. Very much appreciate your contributions!

## Development

Quickstart after cloning the repository:

```bash
python -m venv .venv
source ./venv/bin/activate
python -m pip install -r requirments.txt
```

* Docstrings should be in [numpy style](https://numpydoc.readthedocs.io/en/latest/format.html#documenting-classes)
* Code should be formatted with [black](https://github.com/psf/black)
* The main code is in `src/kimmdy/`
* Reaction plugins are in `plugins/`
* Releases with semantic versioning are created based on [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/#summary)
* Python scripts with command line interfaces are in src/kimmdy/cmd and must be registered in setup.cfg
* requirements.txt containes packages necessary for development, like for testing and linting
* pytest is used for tests, these are located in `tests/`
* tox is used to test automated against multiple python versions
* [quartodoc](https://github.com/machow/quartodoc) is used to generate the documentation website.
  Run `make setup-docs` (once) and/or `make docs` to build and render or `make preview` to build and preview.

