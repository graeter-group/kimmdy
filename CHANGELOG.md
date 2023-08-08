# Changelog

## [3.5.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.4.0...v3.5.0) (2023-08-08)


### Features

* **ci:** first try! ([f2c686f](https://github.com/hits-mbm-dev/kimmdy/commit/f2c686f5336f5e3e5e0c819990c0fb89bec9f647))


### Bug Fixes

* ci ([3bcc81a](https://github.com/hits-mbm-dev/kimmdy/commit/3bcc81a3231ddc2e4a6b7b3d8f112228e2c10226))
* **ci:** almost did it first try... ([2850379](https://github.com/hits-mbm-dev/kimmdy/commit/2850379ca5c0d3d0e1e755a1e2170caa4bef38b8))

## [3.4.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.3.1...v3.4.0) (2023-08-08)


### Features

* add gmx_mdrun_flags option. fixes [#145](https://github.com/hits-mbm-dev/kimmdy/issues/145) ([#161](https://github.com/hits-mbm-dev/kimmdy/issues/161)) ([e780402](https://github.com/hits-mbm-dev/kimmdy/commit/e780402e33e796d2672fc286a81031154ecd1569))
* **config:** validation, completion and defaults for kimmdy.yml ([#153](https://github.com/hits-mbm-dev/kimmdy/issues/153)) ([de7c78f](https://github.com/hits-mbm-dev/kimmdy/commit/de7c78f30de9fb5d5075bff239e82888abc0e35d))
* example from tests ([#175](https://github.com/hits-mbm-dev/kimmdy/issues/175)) ([a096a45](https://github.com/hits-mbm-dev/kimmdy/commit/a096a45af9a3cd0f7974301a307265404b105424))


### Bug Fixes

* consistent use of `edis` in homolysis config ([#170](https://github.com/hits-mbm-dev/kimmdy/issues/170)) ([e847707](https://github.com/hits-mbm-dev/kimmdy/commit/e8477074810a523f2cde79f49bc3cc103fa869a5))
* don't allow `Move` w/o ix/id_to_move ([10c9231](https://github.com/hits-mbm-dev/kimmdy/commit/10c9231f6dd98e2d51613ad94a2c690058bb84b2))
* inconsistent inheritance ([#154](https://github.com/hits-mbm-dev/kimmdy/issues/154)) ([2bec671](https://github.com/hits-mbm-dev/kimmdy/commit/2bec67167207bf3c38e2985e5bef7ed1f3a60870))
* load plugins only once in __init__ ([#171](https://github.com/hits-mbm-dev/kimmdy/issues/171)) ([88b55ee](https://github.com/hits-mbm-dev/kimmdy/commit/88b55ee281f0e9c39b5f527aa21eed88c57e5c8a))
* merge dihedrals ([#174](https://github.com/hits-mbm-dev/kimmdy/issues/174)) ([37f6be3](https://github.com/hits-mbm-dev/kimmdy/commit/37f6be3b322cdba4175817d40fcf489be992a1f8))
* scheme keys not plugin names but entrypoints ([b878b69](https://github.com/hits-mbm-dev/kimmdy/commit/b878b6966ead7dee83a21d274d3de0a8e6aa5291))
* test input files now compatible ([eda027a](https://github.com/hits-mbm-dev/kimmdy/commit/eda027a297d0d33bef404907ac7090f5b0b29204))
* test_merge_prm_top ([9d333bc](https://github.com/hits-mbm-dev/kimmdy/commit/9d333bc37423ae007019bc742b5eba36a70192e3))
* update scheme paths after building examples ([5076ca4](https://github.com/hits-mbm-dev/kimmdy/commit/5076ca4c5196a0d1f12e8f8245f6110ed24417cb))


### Documentation

* improve config option descriptions, names and error reporting ([#163](https://github.com/hits-mbm-dev/kimmdy/issues/163)) ([3318bdc](https://github.com/hits-mbm-dev/kimmdy/commit/3318bdc0f16fff99efb7493e03d6f9a086814325))
* **types:** add custom renderer ([3318bdc](https://github.com/hits-mbm-dev/kimmdy/commit/3318bdc0f16fff99efb7493e03d6f9a086814325))

## [3.3.1](https://github.com/hits-mbm-dev/kimmdy/compare/v3.3.0...v3.3.1) (2023-07-12)


### Documentation

* add a lot more docstrings to functions, classes, modules ([#149](https://github.com/hits-mbm-dev/kimmdy/issues/149)) ([a0d2a43](https://github.com/hits-mbm-dev/kimmdy/commit/a0d2a43ca257934222227897e586ba81d7507497))

## [3.3.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.2.0...v3.3.0) (2023-07-05)


### Features

* linked ixs and ids for RecipSteps ([#137](https://github.com/hits-mbm-dev/kimmdy/issues/137)) ([8238df9](https://github.com/hits-mbm-dev/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))
* **topology:** enable multiple dihedrals (unique by their periodicity) for proper dihedrals ([8238df9](https://github.com/hits-mbm-dev/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))


### Bug Fixes

* add default funct for Bonds and BondTypes ([8238df9](https://github.com/hits-mbm-dev/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))
* add more nestable sections in topology parsing ([8238df9](https://github.com/hits-mbm-dev/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))
* match X as wildcards for atomic types ([8238df9](https://github.com/hits-mbm-dev/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))

## [3.2.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.1.0...v3.2.0) (2023-07-03)


### Features

* quartodoc ([f10e531](https://github.com/hits-mbm-dev/kimmdy/commit/f10e53152d7a9ddc2655bb3744972d06d64eb2b2))
* quartodoc ([#142](https://github.com/hits-mbm-dev/kimmdy/issues/142)) ([f10e531](https://github.com/hits-mbm-dev/kimmdy/commit/f10e53152d7a9ddc2655bb3744972d06d64eb2b2))

## [3.1.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.0.1...v3.1.0) (2023-06-30)


### Features

* changemanager [#129](https://github.com/hits-mbm-dev/kimmdy/issues/129) ([087d652](https://github.com/hits-mbm-dev/kimmdy/commit/087d652dba53857c561c8e028133d5f035e099d3))

## [3.0.1](https://github.com/hits-mbm-dev/kimmdy/compare/v3.0.0...v3.0.1) (2023-06-28)


### Bug Fixes

* dihedraltypes sections in ff/topology ([affa546](https://github.com/hits-mbm-dev/kimmdy/commit/affa5468ab537bcb37d7ae99ceded4fcf48cbc60))

## [3.0.0](https://github.com/hits-mbm-dev/kimmdy/compare/v2.4.0...v3.0.0) (2023-06-26)


### ⚠ BREAKING CHANGES

* put restrictions on the way top files are accepted
* deprecate generate_top_from_atom_list

### Features

* handle else macros ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* handle multiple moleculetype sections ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* parse top files with includes ([#124](https://github.com/hits-mbm-dev/kimmdy/issues/124)) ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* resolve includes in gromacs data dir (located relative to gmx executable) if not in cwd ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* stop if no ff found, warn if multiple found ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* topology dictionary helpers to get and set sections ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* udpate top dict with parsed ff parameters ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))


### Bug Fixes

* deprecate generate_top_from_atom_list ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* put restrictions on the way top files are accepted ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* use hat_naive instead of hat_reaction ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* use logging.warning instead of depcrecated warn ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))


### Documentation

* add tutorial section (with markdown support) ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* move qmd to preview topology to docs ([bcb6121](https://github.com/hits-mbm-dev/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))

## [2.4.0](https://github.com/hits-mbm-dev/kimmdy/compare/v2.3.0...v2.4.0) (2023-06-22)


### Features

* homolysis module and homolysis tests ([3267642](https://github.com/hits-mbm-dev/kimmdy/commit/32676421c3d087ed4e38777117b6ffe334883cf3))
* KMC module testing ([ebb771f](https://github.com/hits-mbm-dev/kimmdy/commit/ebb771fbbcd730d70a0c5fc7a879a8e0705d48b8))
* Move coordinates include time ([10bb3a7](https://github.com/hits-mbm-dev/kimmdy/commit/10bb3a78ac6a53e7f48ef3cf37113518987cedf3))
* Use proper HAT plugin ([6a4b0eb](https://github.com/hits-mbm-dev/kimmdy/commit/6a4b0eb4e13a800545b2ff234a30fe6d773439c1))


### Bug Fixes

* allow restarting from checkpoint from python ([2163173](https://github.com/hits-mbm-dev/kimmdy/commit/2163173495ea5e325b2cd3c278df1bb6c100c026))
* correct assignment of charge and mass in ff ([c5f92a6](https://github.com/hits-mbm-dev/kimmdy/commit/c5f92a672ba6932d64a1ee206fd27c13762d8fe4))
* homolysis test ([178590d](https://github.com/hits-mbm-dev/kimmdy/commit/178590d503ff0d58ad6098cc62e6ef150e197dec))
* initialize radical list ([3196747](https://github.com/hits-mbm-dev/kimmdy/commit/3196747bfd5d724e805e5c892c6c74d1964ec86a))
* initialize Topology.radicals ([5e0f335](https://github.com/hits-mbm-dev/kimmdy/commit/5e0f335c0f7472c30bf1e154913d6afdfc514b6e))
* remove deprecated imports in tests ([d078109](https://github.com/hits-mbm-dev/kimmdy/commit/d078109062a4960984591d283be9de647212a284))
* rename variables in test_kmc ([4df479d](https://github.com/hits-mbm-dev/kimmdy/commit/4df479d4a8c1f0c656368322d02606b5f74a5044))
* test_integration_move_top ([e0f71b6](https://github.com/hits-mbm-dev/kimmdy/commit/e0f71b6506ecab1c14742568f8c47c7a30a15418))


### Documentation

* **general:** add make preview command ([9806e63](https://github.com/hits-mbm-dev/kimmdy/commit/9806e63171cec8d4fc3532bbbab13499d7084089))
* render docs ([f306a6c](https://github.com/hits-mbm-dev/kimmdy/commit/f306a6c5a09176703c6fe12d3df8a51646fe5a81))

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
