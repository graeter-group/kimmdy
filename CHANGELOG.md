# Changelog

## [2.3.0](https://github.com/hits-mbm-dev/kimmdy/compare/v2.2.0...v2.3.0) (2023-05-11)


### Features

* add MOVE conversion for HAT ([0dbe4b1](https://github.com/hits-mbm-dev/kimmdy/commit/0dbe4b1b53272bd9c813550cae875134bc679f55))
* make periodicity part of the id for dihedrals ([0710086](https://github.com/hits-mbm-dev/kimmdy/commit/0710086473a81878d86ed8ca4daa2b31b921e868))
* parse dihedral and position restraints ([1209714](https://github.com/hits-mbm-dev/kimmdy/commit/12097148ff8fe64b4f40eb0d19ad30def6070c7f))
* test for radicals based on bondorder in top initialization ([a446a70](https://github.com/hits-mbm-dev/kimmdy/commit/a446a70348bc575b5e5212a23301c51cd9c41233))


### Bug Fixes

* more robust fix for [#93](https://github.com/hits-mbm-dev/kimmdy/issues/93) ([3cea313](https://github.com/hits-mbm-dev/kimmdy/commit/3cea31394cb7c4abb02439796f724642e4b42e2b))
* rtp parsing for aminoacids.rtp ([025e6c1](https://github.com/hits-mbm-dev/kimmdy/commit/025e6c14fc50982db7809079b8418644e3d474df))


### Documentation

* Updates reaction api documentation ([87007b4](https://github.com/hits-mbm-dev/kimmdy/commit/87007b485c23262451fe86cb1fb8a9ac0a90fce6)), closes [#114](https://github.com/hits-mbm-dev/kimmdy/issues/114)

## [2.2.0](https://github.com/hits-mbm-dev/kimmdy/compare/v2.1.0...v2.2.0) (2023-03-07)


### Features

* add --version command ([bf3e183](https://github.com/hits-mbm-dev/kimmdy/commit/bf3e183ee57c7b78b3913070bde861a198adc441))


### Bug Fixes

* update importlib syntax and and importlib dependency ([abc89e0](https://github.com/hits-mbm-dev/kimmdy/commit/abc89e00cbdc8e4bc5a235bd4a5f8e387deb866e)), closes [#93](https://github.com/hits-mbm-dev/kimmdy/issues/93)

## [2.1.0](https://github.com/hits-mbm-dev/kimmdy/compare/v2.0.0...v2.1.0) (2023-02-27)


### Features

* add kimmdy checkpoint files ([b468378](https://github.com/hits-mbm-dev/kimmdy/commit/b468378bf57219e9b55e7a15b164b82d5fe4a7d3))
* test CI ([d3dff98](https://github.com/hits-mbm-dev/kimmdy/commit/d3dff988f0810e3035a4b06c20f842f98135a390))
* write task file history to separate file ([bdefba1](https://github.com/hits-mbm-dev/kimmdy/commit/bdefba17dc21a6ad20d79031f70f6c120da1ce9e))


### Bug Fixes

* [#84](https://github.com/hits-mbm-dev/kimmdy/issues/84) ([1bef009](https://github.com/hits-mbm-dev/kimmdy/commit/1bef00994e841c77e5510ad5a5d807ec8b61a76a))
* 23 ([a6eb2b7](https://github.com/hits-mbm-dev/kimmdy/commit/a6eb2b7c4c57ad51fd56de0bc95d1c4f4219144c))
* abort kimmdy if gromacs fails ([080ecc0](https://github.com/hits-mbm-dev/kimmdy/commit/080ecc098b09feeac4b166e5639e6ccf26982650)), closes [#83](https://github.com/hits-mbm-dev/kimmdy/issues/83)
* allow existing directories of starting from checkpoint ([8370713](https://github.com/hits-mbm-dev/kimmdy/commit/83707135c42534b31824bac5393d1d24fcb1642d))
* correctly apply atomtype and residue to jumping H ([dbdd21e](https://github.com/hits-mbm-dev/kimmdy/commit/dbdd21e629a9cadac40e4573bd9dde5fc453f97e)), closes [#84](https://github.com/hits-mbm-dev/kimmdy/issues/84)
* don't create empty task files if no files are changed ([615d11c](https://github.com/hits-mbm-dev/kimmdy/commit/615d11cb582da0604d005433ee340eaa208c98e6))
* don't create logfile christmas trees. ([fca231f](https://github.com/hits-mbm-dev/kimmdy/commit/fca231f5805ffbfbdf34a043a3cf5b84fc9c7bbf)), closes [#82](https://github.com/hits-mbm-dev/kimmdy/issues/82)
* make plumed really optional ([beb09cd](https://github.com/hits-mbm-dev/kimmdy/commit/beb09cde189ce3f110537def1bc3450c75e4500f))
* minimization path to mdp ([971453b](https://github.com/hits-mbm-dev/kimmdy/commit/971453b8f2fa34f15cbad1d4256a6bd0cf8408ee))
* use output file dir for histfile name ([a069c3c](https://github.com/hits-mbm-dev/kimmdy/commit/a069c3cce53987eeb3966373361d3d3054753239))
