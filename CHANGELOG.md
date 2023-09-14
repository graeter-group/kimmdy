# Changelog

## [4.4.1](https://github.com/hits-mbm-dev/kimmdy/compare/v4.4.0...v4.4.1) (2023-09-13)


### Bug Fixes

* **test:** increase hypothesis deadline ([f7d8867](https://github.com/hits-mbm-dev/kimmdy/commit/f7d8867560600c04099da99a28be12311cb3b0af))

## [4.4.0](https://github.com/hits-mbm-dev/kimmdy/compare/v4.3.0...v4.4.0) (2023-09-13)


### Features

* **dummy:** this is not a feature, just a release trigger ([6aa1190](https://github.com/hits-mbm-dev/kimmdy/commit/6aa1190fb0eca93388a8e9d7e95c8b2f5945a55d))


### Bug Fixes

* **analysis:** remove y axis label for steps panel in plot_energy ([a2bc88a](https://github.com/hits-mbm-dev/kimmdy/commit/a2bc88ac4b8405ae940f0f241fe4778c2039552e))
* **test:** increase hypothesis deadline ([4b9823d](https://github.com/hits-mbm-dev/kimmdy/commit/4b9823d45265aecfada2d7bd1f4b90640164d31a))
* **test:** increase hypothesis deadline ([9810e84](https://github.com/hits-mbm-dev/kimmdy/commit/9810e8478375e43311e0458d71cc9748ae80d684))

## [4.3.0](https://github.com/hits-mbm-dev/kimmdy/compare/v4.2.1...v4.3.0) (2023-09-08)


### Features

* unify output location ([#255](https://github.com/hits-mbm-dev/kimmdy/issues/255))! ([cbd4d5d](https://github.com/hits-mbm-dev/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))


### Bug Fixes

* closes  [#265](https://github.com/hits-mbm-dev/kimmdy/issues/265) , cleanup with flake8 ([d79adbc](https://github.com/hits-mbm-dev/kimmdy/commit/d79adbcd798cbe56ded35a96c09679eaf5e4b2c2))
* jobscript and logging ([cbd4d5d](https://github.com/hits-mbm-dev/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))
* place_atom ([0e0dc1d](https://github.com/hits-mbm-dev/kimmdy/commit/0e0dc1dbf4537f12b5af3cca32aba11c7a1ba4fc))
* properly merge mds general '.*' section and general options with defaults ([cbd4d5d](https://github.com/hits-mbm-dev/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))
* set defaults on non-existent sections. fixes [#264](https://github.com/hits-mbm-dev/kimmdy/issues/264) ([cbd4d5d](https://github.com/hits-mbm-dev/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))

## [4.2.1](https://github.com/hits-mbm-dev/kimmdy/compare/v4.2.0...v4.2.1) (2023-09-07)


### Bug Fixes

* apply_recipe naming ([6573886](https://github.com/hits-mbm-dev/kimmdy/commit/65738861ccc4c9418f4a1c7da3fd5f56afd33d14))
* handle `Relax` in plot_rates ([3099210](https://github.com/hits-mbm-dev/kimmdy/commit/3099210335a38504dced4aa5ceffaeaa23617dd3))

## [4.2.0](https://github.com/hits-mbm-dev/kimmdy/compare/v4.1.0...v4.2.0) (2023-09-04)


### Features

* **analysis:** open-plot and open-vmd options ([d306e66](https://github.com/hits-mbm-dev/kimmdy/commit/d306e66b893a72280bffbe12e26ab05c8777a8a9))
* plot multiple energy terms ([d306e66](https://github.com/hits-mbm-dev/kimmdy/commit/d306e66b893a72280bffbe12e26ab05c8777a8a9))
* review analysis ([#245](https://github.com/hits-mbm-dev/kimmdy/issues/245)) ([d306e66](https://github.com/hits-mbm-dev/kimmdy/commit/d306e66b893a72280bffbe12e26ab05c8777a8a9))


### Bug Fixes

* [#236](https://github.com/hits-mbm-dev/kimmdy/issues/236) rm concat cmd option ([#240](https://github.com/hits-mbm-dev/kimmdy/issues/240)) ([675ded5](https://github.com/hits-mbm-dev/kimmdy/commit/675ded537a6b241b59b02560c286945e76916f22))
* missing test file ([1e805c1](https://github.com/hits-mbm-dev/kimmdy/commit/1e805c19a13d4f590a0c8f8e7841428f64ef1c71))
* no grappa in tox tests ([feaff46](https://github.com/hits-mbm-dev/kimmdy/commit/feaff46e0d44316cf5cbf46c8045d9e0507850ed))
* Recipe aggregation broken, closes [#252](https://github.com/hits-mbm-dev/kimmdy/issues/252) ([#253](https://github.com/hits-mbm-dev/kimmdy/issues/253)) ([a6177e3](https://github.com/hits-mbm-dev/kimmdy/commit/a6177e3d99af093700167d03311718c9ce571f63))

## [4.1.0](https://github.com/hits-mbm-dev/kimmdy/compare/v4.0.1...v4.1.0) (2023-08-28)


### Features

* --debug for post-mortem debugging ([#214](https://github.com/hits-mbm-dev/kimmdy/issues/214)) ([91360d1](https://github.com/hits-mbm-dev/kimmdy/commit/91360d156a86094bb46111b8a3858404bc3d8daf))
* add remove_hydrogen tool ([#213](https://github.com/hits-mbm-dev/kimmdy/issues/213)) ([b7f7e78](https://github.com/hits-mbm-dev/kimmdy/commit/b7f7e78273458cc89532f5db8c17e58f1d6c6e50))
* Automatic Radical Parameterization ([#144](https://github.com/hits-mbm-dev/kimmdy/issues/144)) ([4b4641b](https://github.com/hits-mbm-dev/kimmdy/commit/4b4641bce00812be0859acd8200fce658577a316))
* automatically stop and restart KIMMDY for joblength-restricted HPC clusters ([#162](https://github.com/hits-mbm-dev/kimmdy/issues/162)) ([602d864](https://github.com/hits-mbm-dev/kimmdy/commit/602d86483aff2f6bf039f4bd74d018e63cd9b7c4))
* cleanup cmd.py ([#228](https://github.com/hits-mbm-dev/kimmdy/issues/228)) ([f77b6dc](https://github.com/hits-mbm-dev/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* config: use_plumed in md section ([#234](https://github.com/hits-mbm-dev/kimmdy/issues/234)) ([22dfacb](https://github.com/hits-mbm-dev/kimmdy/commit/22dfacbe38bf167ab79b391b40ba1883cf090881))
* multimolecule topology ([#211](https://github.com/hits-mbm-dev/kimmdy/issues/211)) ([153af8c](https://github.com/hits-mbm-dev/kimmdy/commit/153af8ce6780dbff928160c89ff17a5a3c373969))
* single reaction tasks in config ([#215](https://github.com/hits-mbm-dev/kimmdy/issues/215)) ([1c07689](https://github.com/hits-mbm-dev/kimmdy/commit/1c076897f429673dbe0d7dfec9caaf86a7018a71))
* workflow for local testing ([092b89e](https://github.com/hits-mbm-dev/kimmdy/commit/092b89e5fcb1addcd2a3d09f9104c0ddc2b3c3eb))


### Bug Fixes

* 212 ([016aa56](https://github.com/hits-mbm-dev/kimmdy/commit/016aa56a630f64b953cb6cab7c9a25584c47a692))
* allow reading topology with an empty forcefield ([#238](https://github.com/hits-mbm-dev/kimmdy/issues/238)) ([a0824c3](https://github.com/hits-mbm-dev/kimmdy/commit/a0824c37379f0ebdf87bbf73dd50a0b114bbeafe))
* also fix [#200](https://github.com/hits-mbm-dev/kimmdy/issues/200) ([f2b2be2](https://github.com/hits-mbm-dev/kimmdy/commit/f2b2be2169bc7f91eeb925cf0e125a321733cb1c))
* catch plugin loading exception ([268f77f](https://github.com/hits-mbm-dev/kimmdy/commit/268f77fdda7abee65d8347f6a14b8cfb5be57274))
* **ci:** render docs ([a62f351](https://github.com/hits-mbm-dev/kimmdy/commit/a62f35162bbbd4b10f9f7f39138a52734cff65c7))
* consolidate cmd interface with previous changes by [@ehhartmann](https://github.com/ehhartmann) ([f77b6dc](https://github.com/hits-mbm-dev/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* **docs:** getting started molstar resources ([f77b6dc](https://github.com/hits-mbm-dev/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* don't increment logfile if starting from a checkpoint ([f77b6dc](https://github.com/hits-mbm-dev/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* don't run tasks twice ([c756793](https://github.com/hits-mbm-dev/kimmdy/commit/c7567932bce51ded317d689e2bec7d74e2dbae4d))
* examples use new plumed field ([9493a45](https://github.com/hits-mbm-dev/kimmdy/commit/9493a4582b57ce6cba065aaa11930a701d5d494c))
* fix docstring references after https://github.com/hits-mbm-dev/kimmdy/pull/229 ([7c4ba75](https://github.com/hits-mbm-dev/kimmdy/commit/7c4ba75e2e72430f0691a4fa9019a87a2f507ee4))
* install order and tox rquirements ([e7d5957](https://github.com/hits-mbm-dev/kimmdy/commit/e7d5957ca1bf757298c6a4eff5f9b4160f3290ef))
* logger instead of logging ([272dd41](https://github.com/hits-mbm-dev/kimmdy/commit/272dd411d31ffbd7bb49ed3a1bd31364e9b4ad0f))
* missing import ([0ac234f](https://github.com/hits-mbm-dev/kimmdy/commit/0ac234f8ce88279e0d48843dc7c71c51d28fd847))
* morse transition rate homolysis tests fails on ci but not locally [#197](https://github.com/hits-mbm-dev/kimmdy/issues/197)  ([#201](https://github.com/hits-mbm-dev/kimmdy/issues/201)) ([f2b2be2](https://github.com/hits-mbm-dev/kimmdy/commit/f2b2be2169bc7f91eeb925cf0e125a321733cb1c))
* new build example ([4c6b93e](https://github.com/hits-mbm-dev/kimmdy/commit/4c6b93e813acc7eedcd45cbc2f86acc2aefb4fb8))
* plugins use task-logger ([#207](https://github.com/hits-mbm-dev/kimmdy/issues/207)) ([d5c419d](https://github.com/hits-mbm-dev/kimmdy/commit/d5c419d1dfe741958c0052465ad83c28b58ded65))
* proper usage of __str__ and __repr__ for topology ([#210](https://github.com/hits-mbm-dev/kimmdy/issues/210)) ([6fd070f](https://github.com/hits-mbm-dev/kimmdy/commit/6fd070f6ec2dcafcb0d5329a7ff60a431eee743c))
* removed need for strict bind break order ([#231](https://github.com/hits-mbm-dev/kimmdy/issues/231)) ([a7859ec](https://github.com/hits-mbm-dev/kimmdy/commit/a7859ec1100d8bd89eddbfa4816e2be0ada358e7))
* requirements typo ([2f998d2](https://github.com/hits-mbm-dev/kimmdy/commit/2f998d232b62dbc802ba871fcb91efdb4b72f507))
* runmng logger getting overwritten ([27baee9](https://github.com/hits-mbm-dev/kimmdy/commit/27baee9fd3a716930e99f3e8cb9a677e44e95047))
* show plugins ([dcb92c8](https://github.com/hits-mbm-dev/kimmdy/commit/dcb92c8737b71efaaa8429cf2343d53000bef39c))


### Documentation

* improve layout and linking ([f77b6dc](https://github.com/hits-mbm-dev/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* more expressive docstrings ([f77b6dc](https://github.com/hits-mbm-dev/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))

## [4.0.1](https://github.com/hits-mbm-dev/kimmdy/compare/v4.0.0...v4.0.1) (2023-08-11)


### Documentation

* **ci:** add docs build workflow ([#190](https://github.com/hits-mbm-dev/kimmdy/issues/190)) ([#198](https://github.com/hits-mbm-dev/kimmdy/issues/198)) ([5fdec0b](https://github.com/hits-mbm-dev/kimmdy/commit/5fdec0b4903a578abd2fe75cfadfc9b842168be1))
* **dev:** add information for developers about conventional commits ([5fdec0b](https://github.com/hits-mbm-dev/kimmdy/commit/5fdec0b4903a578abd2fe75cfadfc9b842168be1))

## [4.0.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.6.0...v4.0.0) (2023-08-11)


### ⚠ BREAKING CHANGES

* **ci:** does this work?

### Features

* **ci:** update encode to support unicode ([d9f16ad](https://github.com/hits-mbm-dev/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))
* **ci:** what if this doesn't have a body? ([d9f16ad](https://github.com/hits-mbm-dev/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))
* this is a release-please test ([d9f16ad](https://github.com/hits-mbm-dev/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))


### Bug Fixes

* **ci:** or this? ([d9f16ad](https://github.com/hits-mbm-dev/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))
* **ci:** release please adds all the things ([d9f16ad](https://github.com/hits-mbm-dev/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))

## [3.6.0](https://github.com/hits-mbm-dev/kimmdy/compare/v3.5.1...v3.6.0) (2023-08-11)


### Features

* Analysis Module ([#176](https://github.com/hits-mbm-dev/kimmdy/issues/176)) ([f37d6c8](https://github.com/hits-mbm-dev/kimmdy/commit/f37d6c8fecc97c4f256c4992d63ceaf68633629c))
* **ci:** run test on PR with testthis label ([072b4bf](https://github.com/hits-mbm-dev/kimmdy/commit/072b4bf2d9b8a0d8f9bbddede28701d59b79f6bf))


### Bug Fixes

* [#185](https://github.com/hits-mbm-dev/kimmdy/issues/185) find gromacs data dir in docker container ([#191](https://github.com/hits-mbm-dev/kimmdy/issues/191)) ([7c0b90c](https://github.com/hits-mbm-dev/kimmdy/commit/7c0b90cec39c625f780024ec0e766481113a274d))
* **ci:** use correct event context variable ([33003a9](https://github.com/hits-mbm-dev/kimmdy/commit/33003a93be90470306410071f42a5b3cd880afe0))
* don't run tests on pr twice ([7957269](https://github.com/hits-mbm-dev/kimmdy/commit/79572694e360db4e528a122afecfd2f61b106f6d))
* fix [#187](https://github.com/hits-mbm-dev/kimmdy/issues/187) ([587c426](https://github.com/hits-mbm-dev/kimmdy/commit/587c426621cb9b403bed23ff008f5e3039e628f4))
* merge issues ([5880dbb](https://github.com/hits-mbm-dev/kimmdy/commit/5880dbb67ae029a827dca77222b84a0c9e653cf1))
* plumed_mod updating ([3f815a2](https://github.com/hits-mbm-dev/kimmdy/commit/3f815a21b00aa0b7e3556ad9159c8f84e0b95b2b))
* register pytest mark properly ([6c56ed9](https://github.com/hits-mbm-dev/kimmdy/commit/6c56ed9892a9fe72b83e21b9da5c286289ad8df1))
* remove unused tests ([8632f8b](https://github.com/hits-mbm-dev/kimmdy/commit/8632f8b5bf7fcd3f3da8b5090b01f2d115d897a9))
* test releases in separate workflow to get easy link for a badge ([0a3e11c](https://github.com/hits-mbm-dev/kimmdy/commit/0a3e11cd6d08d35dfd7b8c781ecf9f07e3614b01))
* tests compatible with local gh act ([71d8204](https://github.com/hits-mbm-dev/kimmdy/commit/71d8204c60e9134af21663c8fc6abe7a05a06c8c))
* trigger workflow ([e86dccd](https://github.com/hits-mbm-dev/kimmdy/commit/e86dccd13fd06be4f27a4265df1b0d4b0fb7f399))

## [3.5.1](https://github.com/hits-mbm-dev/kimmdy/compare/v3.5.0...v3.5.1) (2023-08-08)


### Bug Fixes

* **ci:** testing ([e9b4b19](https://github.com/hits-mbm-dev/kimmdy/commit/e9b4b19b5801f342f0e75e5ccd819e5ae55a9473))
* **ci:** trigger another release PR ([9c824a8](https://github.com/hits-mbm-dev/kimmdy/commit/9c824a890628c7995fa40cca77c974d8fa35644d))

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
