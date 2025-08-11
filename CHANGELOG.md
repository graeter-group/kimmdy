# Changelog


## [8.2.0](https://github.com/graeter-group/kimmdy/compare/v8.1.0...v8.2.0) (2025-08-11)


### Features

* include kimmdy-dimerization in plugins extra ([0c4c3c0](https://github.com/graeter-group/kimmdy/commit/0c4c3c020fd0d4fc65f9c30df0304aa837db79c5))
* use all plugins from workspace ([04ce27f](https://github.com/graeter-group/kimmdy/commit/04ce27f8fccfb8e69427397c266bb3f0fc63fbe5))


### Bug Fixes

* **ci:** fix release permissions ([cb9bfe7](https://github.com/graeter-group/kimmdy/commit/cb9bfe7631a882b89b2cbcdd985a69002ae8b8c8))
* pin python version in pyproject.toml ([5d14663](https://github.com/graeter-group/kimmdy/commit/5d14663a526596301e642c0ecbe9ffa54409795d))
* rm libs of previous python version during tox run ([#562](https://github.com/graeter-group/kimmdy/issues/562)) ([03218c7](https://github.com/graeter-group/kimmdy/commit/03218c7628852142e9421ccbc7f0afa51dd8c30d))
* **testing:** TypeError raised during test_parser_invertible fixes [#557](https://github.com/graeter-group/kimmdy/issues/557) ([#558](https://github.com/graeter-group/kimmdy/issues/558)) ([dfe6b03](https://github.com/graeter-group/kimmdy/commit/dfe6b036ea31079297166fd3726b37dae3ff7647))


### Documentation

* tell uv tool which python version to use ([#563](https://github.com/graeter-group/kimmdy/issues/563)) ([be30ae5](https://github.com/graeter-group/kimmdy/commit/be30ae5578b84951b389376e3431a93f80576da1))
* update installation instructions ([#560](https://github.com/graeter-group/kimmdy/issues/560)) ([93f03b9](https://github.com/graeter-group/kimmdy/commit/93f03b91017926cd928e16eabcb424d61c75f4b6))

## [8.1.0](https://github.com/graeter-group/kimmdy/compare/v8.0.0...v8.1.0) (2025-07-21)


### Features

* buse git submodules for official plugins and uv as recommended installation([#554](https://github.com/graeter-group/kimmdy/issues/554)) ([bb26f3c](https://github.com/graeter-group/kimmdy/commit/bb26f3cdb72692b2ebe203447bdc3e25e6c41b05))
* support gromacs 2025 ([#549](https://github.com/graeter-group/kimmdy/issues/549)) ([cd1e58d](https://github.com/graeter-group/kimmdy/commit/cd1e58d06412e1853f74dadb1396b0d2025f8065))


### Bug Fixes

* allow not finding improper dihedrals for newly bound atoms ([#556](https://github.com/graeter-group/kimmdy/issues/556)) ([6379f5b](https://github.com/graeter-group/kimmdy/commit/6379f5bde66bae7701e8ef15024ed1f29db4238c))
* avoid overwriting gro at reaction time ([#551](https://github.com/graeter-group/kimmdy/issues/551)) ([35e5d9a](https://github.com/graeter-group/kimmdy/commit/35e5d9acaff7d455feda9970076cfe9c7973f09d))

## [8.0.0](https://github.com/graeter-group/kimmdy/compare/v7.1.3...v8.0.0) (2025-06-24)


### ⚠ BREAKING CHANGES

* more sensible default to not remove bonds from plumed on breaking
* deprecate edr toplevel config input
* use chosen time as single source of truth
* don't copy input files into setup task dir

### Features

* add config.changer.topology.break_plumed option ([0b25103](https://github.com/graeter-group/kimmdy/commit/0b2510399e08be5c4aa84dbd81b41c205f633d6a))
* add parser for mdp files ([fd406b5](https://github.com/graeter-group/kimmdy/commit/fd406b5a2d2da80da2e03832aa5270b85f23745d))
* add pure xtc mode ([0d3af49](https://github.com/graeter-group/kimmdy/commit/0d3af49de96b4759a01c76f245236ddfdce7f5bf))
* add Topology.from_path classmethod ([43fee89](https://github.com/graeter-group/kimmdy/commit/43fee89230adb8e73a70f7ff732008c0c503a31d))
* always log grompp command as info ([35508a4](https://github.com/graeter-group/kimmdy/commit/35508a43d47f039118981281450aa06ffcc0c52d))
* **analysis:** chose group for centering the output of trjcat ([d7491a8](https://github.com/graeter-group/kimmdy/commit/d7491a83cb9ba15115dcad827cc80027253a4e4c))
* **config:** add config.slurm.partition ([4f288b1](https://github.com/graeter-group/kimmdy/commit/4f288b14c3fb13bcd0c6fbb5c384c68ceb5083d4))
* copy all md mdps to setup ([ad71893](https://github.com/graeter-group/kimmdy/commit/ad71893e691582d7f27d16f43c96f2ff35d0c107))
* copy small config files such as mdp to 0_setup and write mdps to md ([733604d](https://github.com/graeter-group/kimmdy/commit/733604d7c3d575763302bf2f6c78aa4031a5db00))
* don't copy input files into setup task dir ([7d1fdc6](https://github.com/graeter-group/kimmdy/commit/7d1fdc6b1fe44c7f48cacf773ca132325976331a))
* helper functions to use a harmonic potential for energy barriers instead of morse ([10da947](https://github.com/graeter-group/kimmdy/commit/10da947d92a9c5fb5deb225e5fa915134cac87de))
* **hpc:** use name of output directory as part of the jobscript name ([e8145c2](https://github.com/graeter-group/kimmdy/commit/e8145c229b2b313b16fa388b16debdb1d3973883))
* just_sample mode (see config) without recipe application ([fac4b42](https://github.com/graeter-group/kimmdy/commit/fac4b42bf92ef7d7f56109abfbccd730089009a0))
* more sensible default to not remove bonds from plumed on breaking ([dd6defc](https://github.com/graeter-group/kimmdy/commit/dd6defc50690a85a0dd18903b30fd1a2c4ff152b))
* parse mdp values into python datatypes and make keys canonical ([e9e0f38](https://github.com/graeter-group/kimmdy/commit/e9e0f385219ee18cfd5c881bd46483bdf02a8364))
* round chosen time of kmc result to nearest trr frame. closes [#522](https://github.com/graeter-group/kimmdy/issues/522) ([d732f9a](https://github.com/graeter-group/kimmdy/commit/d732f9ac8128f0a36d6c860cb818ccb42f533aff))
* save harmonic forces to bondstats ([0ce0a3c](https://github.com/graeter-group/kimmdy/commit/0ce0a3c04ed5c7f10eabc293405f069c16760815))
* specify mdp file for input trr or infer from trr name ([34fd3f7](https://github.com/graeter-group/kimmdy/commit/34fd3f767b78948b21926eb6bba35888645f6dab))
* unify hydrolysis and homolysis resources ([90ae1b8](https://github.com/graeter-group/kimmdy/commit/90ae1b8cc462cdc8982daaf7102de136391ef7f4))
* validate mdp file for slow growth relax md fixes [#523](https://github.com/graeter-group/kimmdy/issues/523) ([a8470f5](https://github.com/graeter-group/kimmdy/commit/a8470f5bd8010080b64347df71e561caa5288496))


### Bug Fixes

* allow reading recipes with more rates from csv (increase field_size_limit of csv reader) ([52119a8](https://github.com/graeter-group/kimmdy/commit/52119a83ed05d63b36b1fef9ad48af715ba7f6f4))
* apply amber fix globally to atoms affected ([#536](https://github.com/graeter-group/kimmdy/issues/536)) ([30a4cc5](https://github.com/graeter-group/kimmdy/commit/30a4cc55cc6ef3e3f9a52a36663e50f0e10dcd72))
* deprecate edr toplevel config input ([5a341f8](https://github.com/graeter-group/kimmdy/commit/5a341f89d6b855c8114d1b276c849c60fc2fdee0))
* don't add pair in bind_bond if atoms are directly bound ([f2a424c](https://github.com/graeter-group/kimmdy/commit/f2a424c28563eb9a7910673d444671b83427104c))
* don't invent dihedrals with duplicate atoms (when rings exist) ([549e7ab](https://github.com/graeter-group/kimmdy/commit/549e7ab2f808475ec726db703aff95a4b659f09c))
* don't lowercase name of the kimmdy run (used in output path) ([4f13d3a](https://github.com/graeter-group/kimmdy/commit/4f13d3aa6655efd2173578b6bcf61db32b0ba61c))
* don't modify distances on place when calculating forces ([4ce43bd](https://github.com/graeter-group/kimmdy/commit/4ce43bd1c12a49b4cd83eaca08b62e7979623834))
* fix floating point errors when rounding time to nearest trr frame ([58830be](https://github.com/graeter-group/kimmdy/commit/58830be4e59d4f03de7298b9b2f47de8fe23a0d1))
* get correct improper dihedrals from previous or next residues when ([c382a3f](https://github.com/graeter-group/kimmdy/commit/c382a3ffbfc42d41da6bafef84991087c804c9ff))
* log grompp and mdrun commands just before actually running ([6d8b1f7](https://github.com/graeter-group/kimmdy/commit/6d8b1f7083132f0bc428ca4e0d2a9445442b1860))
* make get_is_reactive_predicate_f case-insensitive. closes [#513](https://github.com/graeter-group/kimmdy/issues/513) ([4f8407f](https://github.com/graeter-group/kimmdy/commit/4f8407fdb1c95befbdadfab8a5d321aa0bc2bddc))
* make kimmdy log symlink in cwd unique per run ([fa74235](https://github.com/graeter-group/kimmdy/commit/fa742356468a3fb55bc505719cd0ac165ab57183))
* make pulling test shorter to avoid pulling across the pbc ([2bf9e53](https://github.com/graeter-group/kimmdy/commit/2bf9e5313e1cc6ac4da3d2d33e88d466b4e02c3b))
* only check gmx version of slow growth is enabled ([d74bef9](https://github.com/graeter-group/kimmdy/commit/d74bef91a19a1db07c67b6e5952745c80acd6a49))
* put slurm in slurm outuput file name ([3c3ff73](https://github.com/graeter-group/kimmdy/commit/3c3ff7334e117a4210d3d311638e202cdf088d7f))
* remove deprecated skip_truncation config ([99a0975](https://github.com/graeter-group/kimmdy/commit/99a09751751ccceaf294221bf1852485dec7b056))
* reset loggers for multiple runs in the same python session ([5b779de](https://github.com/graeter-group/kimmdy/commit/5b779deaf751f7baf793d0be34c6051306453e44))
* skip additional headers in plumed distances.dat files ([5b71e67](https://github.com/graeter-group/kimmdy/commit/5b71e6704b1260df1d150a69dfdcd533ebf52d96))
* sort improper dihedral keys in the order that gromacs would ([5429d4a](https://github.com/graeter-group/kimmdy/commit/5429d4a97e2b800d76652e5e06ae54593e5a03f3))
* use chosen time as single source of truth ([1160296](https://github.com/graeter-group/kimmdy/commit/11602966bab4800551c4832de6c0a308857141c3))
* use the -time flag of grompp instead of .reaction trr and edr ([d547b2a](https://github.com/graeter-group/kimmdy/commit/d547b2a6724509ba82608cffdffe40a02d86c987)), closes [#531](https://github.com/graeter-group/kimmdy/issues/531) [#519](https://github.com/graeter-group/kimmdy/issues/519)
* **wip:** handle impropers and DNA residutypes when binding ([b7836ab](https://github.com/graeter-group/kimmdy/commit/b7836ab6b126c466b035e4e992a8505f09242606))
* write correct task directory to logmessage when starting task ([fb14d31](https://github.com/graeter-group/kimmdy/commit/fb14d31b86c2c0a467e61c2e8b751ed722b14244))


### Performance Improvements

* modify the "before" topology for relax ([#539](https://github.com/graeter-group/kimmdy/issues/539)) ([62a6282](https://github.com/graeter-group/kimmdy/commit/62a62829fb42ddaa09831be48d622c80963c3dda))

## [7.1.3](https://github.com/graeter-group/kimmdy/compare/v7.1.2...v7.1.3) (2025-03-07)


### Bug Fixes

* render docs on trigger ([72430ea](https://github.com/graeter-group/kimmdy/commit/72430eaca7f8cdea47645a7f00edd1a67ebfc013))
* try different ci triggers ([ec3aaf4](https://github.com/graeter-group/kimmdy/commit/ec3aaf47da450b2848bb341912b7e916b100674b))


### Documentation

* add deprecation warning rendering in docs ([1d98764](https://github.com/graeter-group/kimmdy/commit/1d987648f647eff1b67538d770deca68ee4c1212))

## [7.1.2](https://github.com/graeter-group/kimmdy/compare/v7.1.1...v7.1.2) (2025-03-05)


### Bug Fixes

* **ci:** fix building docs before running tests ([#527](https://github.com/graeter-group/kimmdy/issues/527)) ([2b9f969](https://github.com/graeter-group/kimmdy/commit/2b9f969d22956cf64ae600472f9f4cfa8a59d735))
* **ci:** only trigger docs manually ([2f3c4ee](https://github.com/graeter-group/kimmdy/commit/2f3c4ee3a585f88febc5be13886075a54be910da))

## [7.1.1](https://github.com/graeter-group/kimmdy/compare/v7.1.0...v7.1.1) (2025-03-05)


### Bug Fixes

* **ci:** build docs on publish and publish to pypi after release commit ([7c9af68](https://github.com/graeter-group/kimmdy/commit/7c9af68b8be3ecfbbbcc65c808ec3071471f5593))

## [7.1.0](https://github.com/graeter-group/kimmdy/compare/v7.0.0...v7.1.0) (2025-03-05)


### Features

* **config:** Add prefixes for running KIMMDY on HPC Clusters that require srun/mpirun for gromacs binaries ([#517](https://github.com/graeter-group/kimmdy/issues/517)) ([a5d4ace](https://github.com/graeter-group/kimmdy/commit/a5d4acea552503a1e1d8f13af5c90128a7fb15a3))

## [7.0.0](https://github.com/graeter-group/kimmdy/compare/v6.13.0...v7.0.0) (2025-03-04)


### ⚠ BREAKING CHANGES

* merge `slow_growth_pairs` option into `slow_growth` enum
* **top:** make slow growth relaxation more robust for complex topology changes (e.g. hydrolysis)

### Features

* allow deprecation notices in yaml schema ([89d9828](https://github.com/graeter-group/kimmdy/commit/89d98287c02392d2bea8142aeb2d9a9d58cd66b6))
* **analysis:** `trjcat` gained the `use_last_names` flag ([1e636bc](https://github.com/graeter-group/kimmdy/commit/1e636bca0f968585db4c21850e71ec991b2da283))
* flag to skip truncation for fast rate calculation without follow-up ([#500](https://github.com/graeter-group/kimmdy/issues/500)) ([fe031fa](https://github.com/graeter-group/kimmdy/commit/fe031fad7d4c94be21ab7e5ec1880e82c2731939))
* initial edr file config for starting directly with reactions ([985b52d](https://github.com/graeter-group/kimmdy/commit/985b52d8eb0738f38dbadb2960f0389a7e1f4624))
* The topology now keeps track of whether it needs parameterization and the focus_ids for said ([ac4a195](https://github.com/graeter-group/kimmdy/commit/ac4a1956fbe96abbe36d0e000aeee9412cb98d73))
* **top:** make slow growth relaxation more robust for complex topology changes (e.g. hydrolysis) ([772c0d1](https://github.com/graeter-group/kimmdy/commit/772c0d165a6af6b4a7dea63d903d9ce735a7bddb))
* **top:** The `atomic` types like `Atom`, `Bond` etc. gained a `comment` field ([b327268](https://github.com/graeter-group/kimmdy/commit/b32726862828f6c650885a679aa948d74a729092))
* write 1-frame edr file at reaction time ([5b2557c](https://github.com/graeter-group/kimmdy/commit/5b2557ca4a708b2e2e242881a4fbb0fa77772c3e))


### Bug Fixes

* always write reaction coordinates from trr file if it exists ([280ed83](https://github.com/graeter-group/kimmdy/commit/280ed8343baf949352e697d90645431e1598e0d4))
* correct yaml hint for restart flag ([f1edaa0](https://github.com/graeter-group/kimmdy/commit/f1edaa0d8da5bd29bbf52cc844baff8453af8c01))
* fix --use-last-names option for trjcat analysis ([901e0b0](https://github.com/graeter-group/kimmdy/commit/901e0b0860b4ed980198e3029cd7e4bb045bf7e3))
* fix [#518](https://github.com/graeter-group/kimmdy/issues/518) ([2832c04](https://github.com/graeter-group/kimmdy/commit/2832c04511f878a369bc06d7abdbd610f8b14389))
* fix jobscript for re-submission ([caa3383](https://github.com/graeter-group/kimmdy/commit/caa3383b6c07ed81495c4d79ff815235f780e513))
* fix removing failed task dir after prompt ([771277c](https://github.com/graeter-group/kimmdy/commit/771277c830fc6a6c5be2bac542e141ec9b198442))
* fix restarting and use edr file of previous MD task for starting next MD ([a8e8f5e](https://github.com/graeter-group/kimmdy/commit/a8e8f5e360fe2d469970c014db39b92729206595))
* only create helper pairs for morphing dihedrals where there are no ([5d2b61b](https://github.com/graeter-group/kimmdy/commit/5d2b61bd1f342387cb4b6ae7496cb0cb08b66713))
* re-link atomics after merging moleculetypes ([3dc445d](https://github.com/graeter-group/kimmdy/commit/3dc445d0fb220f0044cfb2510d9086a0d3a8d08c))
* restarting ([#504](https://github.com/graeter-group/kimmdy/issues/504)) ([151bb2e](https://github.com/graeter-group/kimmdy/commit/151bb2e2397764f2c88afd36ad7987a193903d50))


### Documentation

* **guide:** The topology explanation now includes a before-after of a HAT ([ee5d04b](https://github.com/graeter-group/kimmdy/commit/ee5d04bfd9a9d4d2d7a5c69193c684ea0bcb2b9d))
* **log:** `Recipe`s printed in the log now transparently use 1-based topology IDs. ([ee84094](https://github.com/graeter-group/kimmdy/commit/ee8409413b91f069b028a2b0b27b9ac46541b8d9))
* render all sub-options on the docs site ([db933e8](https://github.com/graeter-group/kimmdy/commit/db933e8a8a47b250ab67fed431831f3cffe02653))


### Code Refactoring

* merge `slow_growth_pairs` option into `slow_growth` enum ([ec16d29](https://github.com/graeter-group/kimmdy/commit/ec16d2903de2ccf1dd239c9608856affd8684b52))

## [6.13.0](https://github.com/graeter-group/kimmdy/compare/v6.12.0...v6.13.0) (2024-10-23)


### Features

* multi rfkmc ([#494](https://github.com/graeter-group/kimmdy/issues/494)) ([3e64242](https://github.com/graeter-group/kimmdy/commit/3e64242b90bd4a098d59b8bde7a66e5e982f1f37))


### Bug Fixes

* 495 ([#496](https://github.com/graeter-group/kimmdy/issues/496)) ([e1d3700](https://github.com/graeter-group/kimmdy/commit/e1d370009108e8f8ccc7c75d1f99f9a126db3210))
* delete last atom of a chain ([9ca51e3](https://github.com/graeter-group/kimmdy/commit/9ca51e300a2e9241e157b97b6cc32db55c81f3da))
* need parametrization after delete ([c3b5793](https://github.com/graeter-group/kimmdy/commit/c3b579314a5d242d2314b118fb28b13f8e3447f6))
* only check gmx version if applicable ([0d7e96c](https://github.com/graeter-group/kimmdy/commit/0d7e96c8d3b9f8bd543818e0c987d776fd1a4d3e))
* search for radicals only in REACTIVE molecule ([e0cd397](https://github.com/graeter-group/kimmdy/commit/e0cd397f9f9f3b2977e3603d5867e53061ad3f58))

## [6.12.0](https://github.com/graeter-group/kimmdy/compare/v6.11.0...v6.12.0) (2024-09-16)


### Features

* expose create_plumed_input ([#483](https://github.com/graeter-group/kimmdy/issues/483)) ([505c823](https://github.com/graeter-group/kimmdy/commit/505c823afdd523f3c7170be270f06c8aa1fc18bf))


### Bug Fixes

* copy gro to setup ([#487](https://github.com/graeter-group/kimmdy/issues/487)) ([25ac265](https://github.com/graeter-group/kimmdy/commit/25ac2650800c6a07b051b815b4b1b60308cc498e))
* do not log empty edissoc lines ([774d338](https://github.com/graeter-group/kimmdy/commit/774d33820dd70110ce6147e1a976ff0c8c726e2d))
* link residuetypedat ([#488](https://github.com/graeter-group/kimmdy/issues/488)) ([2576440](https://github.com/graeter-group/kimmdy/commit/2576440e3360f6a87257019eda56385c075ccafd))
* truncate time in plumed output to 1 fs ([#484](https://github.com/graeter-group/kimmdy/issues/484)) ([b73104a](https://github.com/graeter-group/kimmdy/commit/b73104a5e81974c7864703fa5e2162f6b7d4575d))
* use correct distances.dat when starting from existing traj ([#485](https://github.com/graeter-group/kimmdy/issues/485)) ([d0c9859](https://github.com/graeter-group/kimmdy/commit/d0c985967fbd9d879ae2f0c813788a91c9492779))

## [6.11.0](https://github.com/graeter-group/kimmdy/compare/v6.10.3...v6.11.0) (2024-09-04)


### Features

* recipes with deferred RecipeSteps ([#476](https://github.com/graeter-group/kimmdy/issues/476)) ([b953130](https://github.com/graeter-group/kimmdy/commit/b95313021c0e9f869e1772fdd50de771d9d1fadc))


### Bug Fixes

* **ci:** always mark draft PRs as not green ([#473](https://github.com/graeter-group/kimmdy/issues/473)) ([bb8ce10](https://github.com/graeter-group/kimmdy/commit/bb8ce105f7aefbcdf40a6791c6f671a228b7a025))
* **ci:** we don't actually need the fail-if-draft action ([b953130](https://github.com/graeter-group/kimmdy/commit/b95313021c0e9f869e1772fdd50de771d9d1fadc))

## [6.10.3](https://github.com/graeter-group/kimmdy/compare/v6.10.2...v6.10.3) (2024-08-20)


### Bug Fixes

* run tests on non-draft PRs, render docs manually with label ([#470](https://github.com/graeter-group/kimmdy/issues/470)) ([166f4d0](https://github.com/graeter-group/kimmdy/commit/166f4d08c8e889baf11afd9b75a4a299ba193df1))

## [6.10.2](https://github.com/graeter-group/kimmdy/compare/v6.10.1...v6.10.2) (2024-08-19)


### Bug Fixes

* pin griffe version for quartodoc ([d4b84ef](https://github.com/graeter-group/kimmdy/commit/d4b84ef1c5bb339aea419d832731d51ff7a84f90))

## [6.10.1](https://github.com/graeter-group/kimmdy/compare/v6.10.0...v6.10.1) (2024-08-19)


### Bug Fixes

* **ci:** build docs on release ([ca8fdc2](https://github.com/graeter-group/kimmdy/commit/ca8fdc297add568adff6d07db59a2e430d17eb1c))

## [6.10.0](https://github.com/graeter-group/kimmdy/compare/v6.9.1...v6.10.0) (2024-08-16)


### Features

* handle `Settles` and `Exclusions` when breaking and binding bonds in topology ([#458](https://github.com/graeter-group/kimmdy/issues/458)) ([4b2d96f](https://github.com/graeter-group/kimmdy/commit/4b2d96f10e79323a12ce4caeda0f629c003b029e))


### Bug Fixes

* **ci:** fix release tests ([1cf898b](https://github.com/graeter-group/kimmdy/commit/1cf898bd608a408832c5087a0e769ba961375c1d))
* copy user supplied input data to setup dir ([4b2d96f](https://github.com/graeter-group/kimmdy/commit/4b2d96f10e79323a12ce4caeda0f629c003b029e))
* simplify ci ([765ed7b](https://github.com/graeter-group/kimmdy/commit/765ed7b98dd1226e1796a31860f6786683a1bda0))

## [6.9.1](https://github.com/graeter-group/kimmdy/compare/v6.9.0...v6.9.1) (2024-08-09)


### Bug Fixes

* bug in output file discovery ([9336a91](https://github.com/graeter-group/kimmdy/commit/9336a9111cb6167ba2fa6aecd48ce5272c16cb20))
* drop dill requirement ([ff3a2c2](https://github.com/graeter-group/kimmdy/commit/ff3a2c241faa5bb9e80065e9b7dcb83f04bc7dc7))
* more stable slow growth parameters ([aedffe0](https://github.com/graeter-group/kimmdy/commit/aedffe050b1a0b34e993da91cb918ea5b0decc7a))
* prevent log spam ([2b51a0e](https://github.com/graeter-group/kimmdy/commit/2b51a0ed204caf6769c1ba899aa540a99f1595e0))
* prevent log spam ([38cbf07](https://github.com/graeter-group/kimmdy/commit/38cbf075071b3c269475d2cac21afa93b8d12195))
* proper pairs in slow-growth ([#465](https://github.com/graeter-group/kimmdy/issues/465)) ([56db2ce](https://github.com/graeter-group/kimmdy/commit/56db2ce14713b5f955f23dbf3b6ec4260a906473))
* remove unreachable code ([af0f11d](https://github.com/graeter-group/kimmdy/commit/af0f11da271ae04ba445db06504fad1965cecc70))
* update task names for restarting changed in [#422](https://github.com/graeter-group/kimmdy/issues/422) ([d31bc2e](https://github.com/graeter-group/kimmdy/commit/d31bc2e28cae7e1ad76614dc66bcb10cbb2af1fe))

## [6.9.0](https://github.com/graeter-group/kimmdy/compare/v6.8.3...v6.9.0) (2024-07-11)


### Features

* add xtc option to config ([95160df](https://github.com/graeter-group/kimmdy/commit/95160df3b38e4328e67f6f14f6672163c0745eb0))
* always setup tasks on init of RunManager ([95160df](https://github.com/graeter-group/kimmdy/commit/95160df3b38e4328e67f6f14f6672163c0745eb0))
* configure logger from init of Config ([95160df](https://github.com/graeter-group/kimmdy/commit/95160df3b38e4328e67f6f14f6672163c0745eb0))
* link nrexcl of reactive moleculetype to toplevel topology ([#451](https://github.com/graeter-group/kimmdy/issues/451)) ([b100117](https://github.com/graeter-group/kimmdy/commit/b10011724eff1d0ae7572b5d5d9b32b2947ae356))


### Bug Fixes

* close task logger handler after taks has run. fixes [#452](https://github.com/graeter-group/kimmdy/issues/452) ([451f3fd](https://github.com/graeter-group/kimmdy/commit/451f3fd02c2bf942939c2576fd13465916be4463))
* polymer_test ([7485011](https://github.com/graeter-group/kimmdy/commit/7485011fb532b6e8e51813bd1da4f3ebeffe9160))
* properly compare to empty dicts ([95160df](https://github.com/graeter-group/kimmdy/commit/95160df3b38e4328e67f6f14f6672163c0745eb0))
* recip parsing without eval ([#448](https://github.com/graeter-group/kimmdy/issues/448)) ([0a5a213](https://github.com/graeter-group/kimmdy/commit/0a5a213254d5c762a3b377305af33f0a62d9c6ab))
* resolve all non-relative paths in config relative to config.cwd ([95160df](https://github.com/graeter-group/kimmdy/commit/95160df3b38e4328e67f6f14f6672163c0745eb0))


### Documentation

* update quarto version ([95160df](https://github.com/graeter-group/kimmdy/commit/95160df3b38e4328e67f6f14f6672163c0745eb0))

## [6.8.1](https://github.com/graeter-group/kimmdy/compare/v6.8.0...v6.8.1) (2024-06-17)

## [6.8.3](https://github.com/graeter-group/kimmdy/compare/v6.8.2...v6.8.3) (2024-06-25)

### Bug Fixes

* more robust output checking for gmx check ([#445](https://github.com/graeter-group/kimmdy/issues/445)) ([0a25161](https://github.com/graeter-group/kimmdy/commit/0a2516158110b3ab16421e513d203006a6cfc076))
* name ([595df32](https://github.com/graeter-group/kimmdy/commit/595df32480d4dd1d2902cef0e5fa9dd0289392bd))

## [6.8.1](https://github.com/graeter-group/kimmdy/compare/v6.8.0...v6.8.1) (2024-06-17)

### Bug Fixes

* check for residuetypes.rtp existence ([c4e831c](https://github.com/graeter-group/kimmdy/commit/c4e831c85b061b1f1916d9a8855908deb7ed6580))
* fix [#440](https://github.com/graeter-group/kimmdy/issues/440) ([c4e831c](https://github.com/graeter-group/kimmdy/commit/c4e831c85b061b1f1916d9a8855908deb7ed6580))
* onlye require xtc or trr, not both ([#439](https://github.com/graeter-group/kimmdy/issues/439)) ([c4e831c](https://github.com/graeter-group/kimmdy/commit/c4e831c85b061b1f1916d9a8855908deb7ed6580))
* use nrexcl from first reactive molecultetype ([c4e831c](https://github.com/graeter-group/kimmdy/commit/c4e831c85b061b1f1916d9a8855908deb7ed6580))

### Documentation

* make docs of kimmdy.yml much more readable ([c4e831c](https://github.com/graeter-group/kimmdy/commit/c4e831c85b061b1f1916d9a8855908deb7ed6580))

## [6.8.0](https://github.com/graeter-group/kimmdy/compare/v6.7.0...v6.8.0) (2024-06-11)


### Features

* parse coarse-grained polymer ffs ([#435](https://github.com/graeter-group/kimmdy/issues/435)) ([ab696dc](https://github.com/graeter-group/kimmdy/commit/ab696dc8ab73c464ff1076c4ec7e3a97f8b778d6))

## [6.7.0](https://github.com/graeter-group/kimmdy/compare/v6.6.1...v6.7.0) (2024-05-20)


### Features

* pass slurm options via kimmdy.yml config ([#432](https://github.com/graeter-group/kimmdy/issues/432)) ([11d3056](https://github.com/graeter-group/kimmdy/commit/11d3056aef19cd71d4ab265daf29e6eaf75457dd))

## [6.6.1](https://github.com/graeter-group/kimmdy/compare/v6.6.0...v6.6.1) (2024-05-07)


### Bug Fixes

* generate matching gro for xtc in trjconv ([800de0a](https://github.com/graeter-group/kimmdy/commit/800de0ac16f1b21af8d41c94718fd34bdf3ab4df))

## [6.6.0](https://github.com/graeter-group/kimmdy/compare/v6.5.0...v6.6.0) (2024-04-25)


### Features

* add explicit includes and excludes for modify_top ([#427](https://github.com/graeter-group/kimmdy/issues/427)) ([c1318c3](https://github.com/graeter-group/kimmdy/commit/c1318c3d55bb21e6ee0e82e32301cc45ca8e7f24))
* Expose grappa model and charge type ([#430](https://github.com/graeter-group/kimmdy/issues/430)) ([57af388](https://github.com/graeter-group/kimmdy/commit/57af38808dd510640cc53420803d2299f5e172f5))
* modify top tut ([#428](https://github.com/graeter-group/kimmdy/issues/428)) ([823b7ff](https://github.com/graeter-group/kimmdy/commit/823b7ffe631c13501266804a90ecd11f8397886b))


### Bug Fixes

* handle time windows w/ no rates in extrande_mod correctly ([f5c3f7c](https://github.com/graeter-group/kimmdy/commit/f5c3f7c76f4c18a7e98cf2fe13afbe1c47041ead))
* typo ([7130aff](https://github.com/graeter-group/kimmdy/commit/7130aff8f09dbe56231b1d2c20ee2b620db8aaca))


### Documentation

* how_to page for analysis ([#424](https://github.com/graeter-group/kimmdy/issues/424)) ([b35e1c8](https://github.com/graeter-group/kimmdy/commit/b35e1c86534dc59090150d2a94c912f9903f0c83))

## [6.5.0](https://github.com/graeter-group/kimmdy/compare/v6.4.1...v6.5.0) (2024-03-25)


### Features

* Restart from run directory([#410](https://github.com/graeter-group/kimmdy/issues/410)) ([40ee0d2](https://github.com/graeter-group/kimmdy/commit/40ee0d20a2a8c704eb2dd3a0898101a46702c695))


### Bug Fixes

* import get_task_directories from utils ([8d4ba4d](https://github.com/graeter-group/kimmdy/commit/8d4ba4de90d33dc9b5b90f450a5d417fcb7d288a))

## [6.4.1](https://github.com/graeter-group/kimmdy/compare/v6.4.0...v6.4.1) (2024-03-21)


### Bug Fixes

* allow reactions with no chosen recipe ([171c7e1](https://github.com/graeter-group/kimmdy/commit/171c7e1df818b8da71198a99e9dbcb0738b7f800))
* colormap normalization ([#415](https://github.com/graeter-group/kimmdy/issues/415)) ([1238397](https://github.com/graeter-group/kimmdy/commit/12383974b1930af38a89ccdcd7b9d1f5c74ed60a))
* energy labels sorted correctly ([855272f](https://github.com/graeter-group/kimmdy/commit/855272f44f965d9240abd2cb529449597badd8da))
* expose arrhenius equation in morse function and fix wrong default value ([#417](https://github.com/graeter-group/kimmdy/issues/417)) ([58c47b3](https://github.com/graeter-group/kimmdy/commit/58c47b31dda698f6005e52e5fc8cf4e7736a5774))
* plot_energy now handles numbers in terms correctly. ([c16888b](https://github.com/graeter-group/kimmdy/commit/c16888bf6fc08715c0819638a2869fd435831d90))
* remove pd warning ([abc4afe](https://github.com/graeter-group/kimmdy/commit/abc4afe11fa0ddb5e35bfdd7889afe1114c0a576))
* runmanager now handles empty kmcresults correctly. ([c94bd2a](https://github.com/graeter-group/kimmdy/commit/c94bd2a2a92765a4e2d88496ea377725fd11d3b5))


### Documentation

* fix docstring ([9572791](https://github.com/graeter-group/kimmdy/commit/9572791ef4e3375bcb82def24c9cf1af2d7c4678))

## [6.4.0](https://github.com/graeter-group/kimmdy/compare/v6.3.0...v6.4.0) (2024-03-07)


### Features

* add focus to grappa parameterizer ([#408](https://github.com/graeter-group/kimmdy/issues/408)) ([9abe754](https://github.com/graeter-group/kimmdy/commit/9abe7546cf0ec9b3c1283c30e0412df5902b2e85))

## [6.3.0](https://github.com/graeter-group/kimmdy/compare/v6.2.5...v6.3.0) (2024-03-07)


### Features

* allow customizing trjconv output group ([#396](https://github.com/graeter-group/kimmdy/issues/396)) ([97b5588](https://github.com/graeter-group/kimmdy/commit/97b55880a791436c22d8c474d2d88541051f831b))
* analysis radical migration tool ([#388](https://github.com/graeter-group/kimmdy/issues/388)) ([ba4cb0f](https://github.com/graeter-group/kimmdy/commit/ba4cb0f12224c6801e21a2961a52e3f3c3172fc1))


### Bug Fixes

* squiggly lines begone ([#399](https://github.com/graeter-group/kimmdy/issues/399)) ([b1210e4](https://github.com/graeter-group/kimmdy/commit/b1210e44f222204a67e1385ce8467ff7e41dabd2))

## [6.2.5](https://github.com/graeter-group/kimmdy/compare/v6.2.4...v6.2.5) (2024-02-27)


### Bug Fixes

* integration tests ([#390](https://github.com/graeter-group/kimmdy/issues/390)) ([526d888](https://github.com/graeter-group/kimmdy/commit/526d888d9198b450a654a5a9de03429a5c6f2e95))

## [6.2.4](https://github.com/graeter-group/kimmdy/compare/v6.2.3...v6.2.4) (2024-02-27)


### Bug Fixes

* commit rendered docs with force ([#391](https://github.com/graeter-group/kimmdy/issues/391)) ([144d5ac](https://github.com/graeter-group/kimmdy/commit/144d5acdfcd3e92371a19be57303eae2f715d8b9))

## [6.2.3](https://github.com/graeter-group/kimmdy/compare/v6.2.2...v6.2.3) (2024-02-27)


### Bug Fixes

* update dihedrals, make pairs unique ([#383](https://github.com/graeter-group/kimmdy/issues/383)) fixes [#376](https://github.com/graeter-group/kimmdy/issues/376) ([9283e68](https://github.com/graeter-group/kimmdy/commit/9283e68065108566c7e0162d7aec543a5b9c41f8))

## [6.2.2](https://github.com/graeter-group/kimmdy/compare/v6.2.1...v6.2.2) (2024-02-23)


### Bug Fixes

* release-please config ([b9a3622](https://github.com/graeter-group/kimmdy/commit/b9a362288ab1b650975d2e072b6f34eb8ca993e3))

## [6.2.1](https://github.com/graeter-group/kimmdy/compare/v6.2.0...v6.2.1) (2024-02-22)


### Features

* analyse reaction participation ([#365](https://github.com/graeter-group/kimmdy/issues/365)) ([c9c7898](https://github.com/graeter-group/kimmdy/commit/c9c7898d429dfbb65a1dabb26cdad7cb2aad1ada))


### Bug Fixes

* merge regression ([9f9b162](https://github.com/graeter-group/kimmdy/commit/9f9b162c40aaf91f691b3963578030592f0ec0a5))


### Documentation

* add how-to for starting with reaction ([#369](https://github.com/graeter-group/kimmdy/issues/369)) ([937c786](https://github.com/graeter-group/kimmdy/commit/937c7864c530d000fb4475d1b33c746265b8c599))
* fix make_ndx command ([69b3d74](https://github.com/graeter-group/kimmdy/commit/69b3d74172c2f50171e85ab56500abef2adc200d))

## [6.2.0](https://github.com/graeter-group/kimmdy/compare/v6.1.0...v6.2.0) (2024-02-16)


### Features

* Generalize ff interaction ([#360](https://github.com/graeter-group/kimmdy/issues/360)) ([49bdb54](https://github.com/graeter-group/kimmdy/commit/49bdb54c4f46909125df954915c7a061b5eaed11))
* improve edissoc ([#366](https://github.com/graeter-group/kimmdy/issues/366)) ([37bb8da](https://github.com/graeter-group/kimmdy/commit/37bb8da0bfe2e1a69a24e3440b702d7c39430e17))

## [6.1.0](https://github.com/graeter-group/kimmdy/compare/v6.0.1...v6.1.0) (2024-02-14)


### Features

* modify_top saves mapping of old and new indices ([#364](https://github.com/graeter-group/kimmdy/issues/364)) ([0086bbe](https://github.com/graeter-group/kimmdy/commit/0086bbe2ff7cab9e900ea879da293f7b1f02757c))

## [6.0.1](https://github.com/graeter-group/kimmdy/compare/v6.0.0...v6.0.1) (2024-01-03)


### Documentation

* minor improvement ([111112e](https://github.com/graeter-group/kimmdy/commit/111112ebcbd3ef79e5a3fea7026c1df9b64402ce))

## [6.0.0](https://github.com/graeter-group/kimmdy/compare/v5.3.2...v6.0.0) (2024-01-03)


### ⚠ BREAKING CHANGES

* Partial charge treatment ([#353](https://github.com/graeter-group/kimmdy/issues/353))

### Features

* Partial charge treatment ([#353](https://github.com/graeter-group/kimmdy/issues/353)) ([1c057b3](https://github.com/graeter-group/kimmdy/commit/1c057b36f3afe8b027eb07b215923f611e72cd0a))


### Bug Fixes

* don't truncate if there is no trajectory. ([4e6e011](https://github.com/graeter-group/kimmdy/commit/4e6e01183ae084b120bf4f0970c327cc9d7f91ee))
* don't try to truncate to time 0 ([52099e5](https://github.com/graeter-group/kimmdy/commit/52099e50adec95022d69c0001b9ddb8a9913975d))
* hotfix for aggregate_reactions hangs [#355](https://github.com/graeter-group/kimmdy/issues/355) ([b5bc285](https://github.com/graeter-group/kimmdy/commit/b5bc285524041338565b54e7dd29cc311d22eca3))

## [5.3.2](https://github.com/graeter-group/kimmdy/compare/v5.3.1...v5.3.2) (2023-12-19)


### Bug Fixes

* Improve slowgrowth ([#338](https://github.com/graeter-group/kimmdy/issues/338)) ([31a0464](https://github.com/graeter-group/kimmdy/commit/31a046447ed50d8a545ec1eda138a30a24dd939f))
* more verbose debugger in rf kmc ([97f4a9a](https://github.com/graeter-group/kimmdy/commit/97f4a9a0e20f370202d2b70a1e7209af7eecc99f))

## [5.3.1](https://github.com/graeter-group/kimmdy/compare/v5.3.0...v5.3.1) (2023-12-15)


### Bug Fixes

* don't write top section headers for empty sections ([#351](https://github.com/graeter-group/kimmdy/issues/351)) ([1632f6f](https://github.com/graeter-group/kimmdy/commit/1632f6f2f78c7c63cf97f3a77821bbef54ea0863))
* update scheme url in examples and docs ([#348](https://github.com/graeter-group/kimmdy/issues/348)) ([2a72b90](https://github.com/graeter-group/kimmdy/commit/2a72b905d6cdd9a914106554de791a1abd346bf8))

## [5.3.0](https://github.com/graeter-group/kimmdy/compare/v5.2.3...v5.3.0) (2023-12-12)


### Features

* custom step ([#346](https://github.com/graeter-group/kimmdy/issues/346)) ([5b28f04](https://github.com/graeter-group/kimmdy/commit/5b28f049dda4ca6e56b959b413778eb73ba3de0d))

## [5.2.3](https://github.com/graeter-group/kimmdy/compare/v5.2.2...v5.2.3) (2023-12-08)


### Bug Fixes

* fix an issue where merging multiples of a moleculetype incremented the residuenumbers too fast ([#344](https://github.com/graeter-group/kimmdy/issues/344)) ([7e78665](https://github.com/graeter-group/kimmdy/commit/7e78665ea1593d65b4b49daa7c1a71a99c87e94d))

## [5.2.2](https://github.com/graeter-group/kimmdy/compare/v5.2.1...v5.2.2) (2023-12-08)


### Bug Fixes

* escape file paths for shell commands ([914550a](https://github.com/graeter-group/kimmdy/commit/914550ae57cf6428845a57cb95b388b33f29ee40))
* make ff dir optional in config ([02c5d9f](https://github.com/graeter-group/kimmdy/commit/02c5d9f3cb157b8ffc1db9ceec9285f667e73d1c))
* replace spaces in name and out in config ([304de69](https://github.com/graeter-group/kimmdy/commit/304de69b4047ad56df6caa64b381f60bc8c16ce0))

## [5.2.1](https://github.com/graeter-group/kimmdy/compare/v5.2.0...v5.2.1) (2023-12-07)


### Bug Fixes

* add title to scheme ([#341](https://github.com/graeter-group/kimmdy/issues/341)) ([a1b1ae5](https://github.com/graeter-group/kimmdy/commit/a1b1ae5f429d86555c3429ee9a4e01f5df618882))
* kimmdy run bugs ([#336](https://github.com/graeter-group/kimmdy/issues/336)) ([ce558e7](https://github.com/graeter-group/kimmdy/commit/ce558e7fb528800f75ac0133e8dde871179bd70c))
* make forcefield directory optional, only warn ([#332](https://github.com/graeter-group/kimmdy/issues/332)) ([eac4131](https://github.com/graeter-group/kimmdy/commit/eac41311e2230bdf014c3df15d141993ac2617ae))
* renumber atomnrs ([#331](https://github.com/graeter-group/kimmdy/issues/331)) ([e88f8cc](https://github.com/graeter-group/kimmdy/commit/e88f8cc809056462249b828960b98565f3768908))

## [5.2.0](https://github.com/graeter-group/kimmdy/compare/v5.1.12...v5.2.0) (2023-11-14)


### Features

* add config for exclude and include in reactive moleculetype ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))
* merge react-able molecules and make multiples explicit ([#327](https://github.com/graeter-group/kimmdy/issues/327)) ([ae8d45e](https://github.com/graeter-group/kimmdy/commit/ae8d45ea714713180989a59043d8c21224697700))
* parse settles and exclusions ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))
* use parsed exclusions for top merge ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))


### Bug Fixes

* declare kimmdy as being a typed library ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))
* expected length of output folders including the new setup task ([8fd2b44](https://github.com/graeter-group/kimmdy/commit/8fd2b4462e664c0ad374eeff161388e6f1010040))
* make molecules section a list instead of dict because order is important ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))
* start iteration at -1 such that setup task has index 0 ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))
* tasklist test after integrating setup task ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))


### Documentation

* improve dark theme ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))
* more detailed ML installation instructions for upcoming plugins ([c128069](https://github.com/graeter-group/kimmdy/commit/c128069beb5940b27d2833de1135560bd7eedd5f))
* use freeze to only re-render changed qmd docs ([44121e4](https://github.com/graeter-group/kimmdy/commit/44121e4c253decb8b2af566e3c72dff10ed8b251))

## [5.1.12](https://github.com/graeter-group/kimmdy/compare/v5.1.11...v5.1.12) (2023-11-13)


### Bug Fixes

* remove show schema path option ([2285246](https://github.com/graeter-group/kimmdy/commit/2285246c355b0f52cc65f1edd24645af7f4e4e9b))

## [5.1.11](https://github.com/graeter-group/kimmdy/compare/v5.1.10...v5.1.11) (2023-11-08)


### Documentation

* clean up docs and error messages ([#325](https://github.com/graeter-group/kimmdy/issues/325)) ([7458940](https://github.com/graeter-group/kimmdy/commit/74589405f998406847f54faa514d290c2c6a4564))

## [5.1.10](https://github.com/graeter-group/kimmdy/compare/v5.1.9...v5.1.10) (2023-11-08)


### Bug Fixes

* discover plugins for building docs ([22018e4](https://github.com/graeter-group/kimmdy/commit/22018e42d0ccb3768379f82d68a5521f1985101e))

## [5.1.9](https://github.com/graeter-group/kimmdy/compare/v5.1.8...v5.1.9) (2023-11-08)


### Bug Fixes

* add dev dependencies to build docs on gh pages ([9bd83df](https://github.com/graeter-group/kimmdy/commit/9bd83dfa8eb17bcfebbca3c0e27462af0e095a50))
* tox, ci ([1e3826a](https://github.com/graeter-group/kimmdy/commit/1e3826afdad54d2e6bb38c7d3b617d74198fc669))

## [5.1.8](https://github.com/graeter-group/kimmdy/compare/v5.1.7...v5.1.8) (2023-11-07)


### Bug Fixes

* building docs ([99506c5](https://github.com/graeter-group/kimmdy/commit/99506c592d0a8e1e25eb9ff2566542f23dd40628))
* tests with plugins ([218adf6](https://github.com/graeter-group/kimmdy/commit/218adf6b66455ec34c1dde5c33394cf5d81e0a05))

## [5.1.7](https://github.com/graeter-group/kimmdy/compare/v5.1.6...v5.1.7) (2023-11-06)


### Bug Fixes

* always write system and molecules top sections last ([c49aa9c](https://github.com/graeter-group/kimmdy/commit/c49aa9c552cfacd605121c39a17846226b84c0c0))

## [5.1.6](https://github.com/graeter-group/kimmdy/compare/v5.1.5...v5.1.6) (2023-10-27)


### Bug Fixes

* **ci:** set markdown content type for description ([6dc116c](https://github.com/graeter-group/kimmdy/commit/6dc116c0568f27437b55acf5c9cedcf0e7ae0eec))

## [5.1.5](https://github.com/graeter-group/kimmdy/compare/v5.1.4...v5.1.5) (2023-10-27)


### Bug Fixes

* **ci:** set markdown content type for description ([ca9963e](https://github.com/graeter-group/kimmdy/commit/ca9963e6a2c4cd5792872980687be281a34270c5))

## [5.1.4](https://github.com/graeter-group/kimmdy/compare/v5.1.3...v5.1.4) (2023-10-27)


### Bug Fixes

* **ci:** markdown in readme ([8a41e39](https://github.com/graeter-group/kimmdy/commit/8a41e39e346e5c6aaaa25dde63eebbc821f9e2b7))

## [5.1.3](https://github.com/graeter-group/kimmdy/compare/v5.1.2...v5.1.3) (2023-10-27)


### Bug Fixes

* **ci:** add testpyi ([390217b](https://github.com/graeter-group/kimmdy/commit/390217b5f29cce34e784437f053f0f8a6e1e1222))
* setup.cfg long description type ([1c31db5](https://github.com/graeter-group/kimmdy/commit/1c31db5ec26399699df0465aa21f401dbda421e5))

## [5.1.2](https://github.com/graeter-group/kimmdy/compare/v5.1.1...v5.1.2) (2023-10-27)


### Bug Fixes

* **ci:** unify release workflow ([488779f](https://github.com/graeter-group/kimmdy/commit/488779fc4271c9112136d5c1e9f4496f433fc0c1))

## [5.1.1](https://github.com/graeter-group/kimmdy/compare/v5.1.0...v5.1.1) (2023-10-27)


### Bug Fixes

* **ci:** only build package on tag push 'stable' ([3a839a3](https://github.com/graeter-group/kimmdy/commit/3a839a32121de8b9c7cf3038a7cd3b16a36c68dc))

## [5.1.0](https://github.com/graeter-group/kimmdy/compare/v5.0.0...v5.1.0) (2023-10-27)


### Features

* adjust resnum of a moving hydrogen ([#312](https://github.com/graeter-group/kimmdy/issues/312)) ([481cd84](https://github.com/graeter-group/kimmdy/commit/481cd84582620a6ee486253638a7ed1db2299934))
* Checkpoints optional ([#311](https://github.com/graeter-group/kimmdy/issues/311)) ([5e93913](https://github.com/graeter-group/kimmdy/commit/5e939139c867f51a7208049a2a94bef49947ebbf))
* **ci:** publish to PyPi ([e513fb6](https://github.com/graeter-group/kimmdy/commit/e513fb645c38ca2b3e791262b8c43c9baa1f8f49))
* **dev:** visualize call graph with Pycallgraph ([#303](https://github.com/graeter-group/kimmdy/issues/303)) ([ab6ffe5](https://github.com/graeter-group/kimmdy/commit/ab6ffe527041d76cedc3cc5edeb893c297ddb47d))
* show gmx errors in log. fixes [#260](https://github.com/graeter-group/kimmdy/issues/260) ([#313](https://github.com/graeter-group/kimmdy/issues/313)) ([489eb2d](https://github.com/graeter-group/kimmdy/commit/489eb2d9c0927db2a5a52712bef2632c15f89037))


### Bug Fixes

* id in scheme ([a597bba](https://github.com/graeter-group/kimmdy/commit/a597bbab89f87e8284c8e4aa4ba188468653ce92))


### Documentation

* add missing config descriptions and remove unused config.run option ([#305](https://github.com/graeter-group/kimmdy/issues/305)) ([1696f66](https://github.com/graeter-group/kimmdy/commit/1696f6679245ca0306d3d6fe021b594cbd11a24b))

## [5.0.0](https://github.com/graeter-group/kimmdy/compare/v4.4.1...v5.0.0) (2023-10-25)


### ⚠ BREAKING CHANGES

* choose highest rate as reaction start time

### Features

* Allow having single optional file for mdrun ([#151](https://github.com/graeter-group/kimmdy/issues/151)) ([6b61a8c](https://github.com/graeter-group/kimmdy/commit/6b61a8c44db785bedcc762bc5689486adc061703))
* choose highest rate as reaction start time ([fc822b3](https://github.com/graeter-group/kimmdy/commit/fc822b32d1b3c34a5f5930574e7daa2ad8f7547a))
* display residue name for radical population ([77a8d72](https://github.com/graeter-group/kimmdy/commit/77a8d726c223d45233c15100437d641c0c1cfc65))
* extrande KMC algorithm ([d1671db](https://github.com/graeter-group/kimmdy/commit/d1671db75d6d86f34d9f925adbd8d42747f05c34))
* only plot radical occupied atoms ([#308](https://github.com/graeter-group/kimmdy/issues/308)) ([77a8d72](https://github.com/graeter-group/kimmdy/commit/77a8d726c223d45233c15100437d641c0c1cfc65))
* Runtime analysis tool ([#283](https://github.com/graeter-group/kimmdy/issues/283)) ([e302be3](https://github.com/graeter-group/kimmdy/commit/e302be358fdf77c749d8c9a369871ad7a3cc353e))
* trr without edr can now be used. ([6b61a8c](https://github.com/graeter-group/kimmdy/commit/6b61a8c44db785bedcc762bc5689486adc061703))


### Bug Fixes

* [#261](https://github.com/graeter-group/kimmdy/issues/261) improve error message if ff is not found in cwd ([#307](https://github.com/graeter-group/kimmdy/issues/307)) ([ab03cbe](https://github.com/graeter-group/kimmdy/commit/ab03cbe8b4c93b1db4e52b78e119de7fb91ece7a))
* after truncating place doesn't need a time ([233d83c](https://github.com/graeter-group/kimmdy/commit/233d83cbba88396194d989a200f9f420e862aae3))
* catch edge cases in extrande algorithm ([01d4df9](https://github.com/graeter-group/kimmdy/commit/01d4df9b0e8fba921337f24cf2b4a1927885e016))
* closes [#286](https://github.com/graeter-group/kimmdy/issues/286) don't modify task mapping inplace ([c40f360](https://github.com/graeter-group/kimmdy/commit/c40f360acd8d5282b4adc4c2eba56777be6f708b))
* consistent naming for remove hydrogen cli ([e615370](https://github.com/graeter-group/kimmdy/commit/e615370a0bec48ead5f66a351f2bc682a0234028))
* coordinates and runmanager ([#271](https://github.com/graeter-group/kimmdy/issues/271)) ([695026c](https://github.com/graeter-group/kimmdy/commit/695026cc76d2311923c6a4b78608c3abee9d45a3))
* correct error message ([facda23](https://github.com/graeter-group/kimmdy/commit/facda238da55f887635728b01afb0a05861e68ef))
* correct truncated gro ([2c49c80](https://github.com/graeter-group/kimmdy/commit/2c49c809de992bee9f6da516fbe76734460b7c2f))
* duplicate key in yaml ([6204e1e](https://github.com/graeter-group/kimmdy/commit/6204e1e52d151ca53af13b122dbce45c6e66acd3))
* error handling ([d0dc8e7](https://github.com/graeter-group/kimmdy/commit/d0dc8e7ee770711f63fa5544190d1269e3f77b2b))
* execute `place` at first time in recipe ([7234d72](https://github.com/graeter-group/kimmdy/commit/7234d721a54e25ef8d43ebb2b7444a914adb37da))
* gmx trjconv --trunc not working with xtc ([ba0b922](https://github.com/graeter-group/kimmdy/commit/ba0b92255b0aecaac67f341d1dfeb9db9a7621fe))
* improved rate plotting ([d0f6d19](https://github.com/graeter-group/kimmdy/commit/d0f6d19a42c705082c8cdd5e498885d426256bca))
* improved reliability of rate plotting ([6a355d4](https://github.com/graeter-group/kimmdy/commit/6a355d439683c44b9fcef68874ff548c6b24d8a4))
* make plugins and plugin schemas discoverable by pip install ([73da927](https://github.com/graeter-group/kimmdy/commit/73da927a01de4570bf5c0571cb1bd07f8c225b10))
* only use mpme and ntmpi in test config. fixes [#302](https://github.com/graeter-group/kimmdy/issues/302) ([#304](https://github.com/graeter-group/kimmdy/issues/304)) ([764b083](https://github.com/graeter-group/kimmdy/commit/764b0830127f0faa2181d15069c405e3385fa3f1))
* parametrize before slow growth ([ff14694](https://github.com/graeter-group/kimmdy/commit/ff146945531795b9179b220db3bf0d97b7b2d870))
* rate plot with less than 9 recipes ([eb557d4](https://github.com/graeter-group/kimmdy/commit/eb557d474b5f829b9ed67e9eb6fdc155113ceaf3))
* remove double log message ([#290](https://github.com/graeter-group/kimmdy/issues/290)) ([0182d20](https://github.com/graeter-group/kimmdy/commit/0182d204636f9ea67bd198f0acd06b687bac928a))
* rename idx to ix for moleculetype for consistency ([e615370](https://github.com/graeter-group/kimmdy/commit/e615370a0bec48ead5f66a351f2bc682a0234028))
* runtime analysis double counting ([14a875a](https://github.com/graeter-group/kimmdy/commit/14a875a3dc6b47550976f1635a74123285c45f9f))
* set needs_parametrization ([e7a2996](https://github.com/graeter-group/kimmdy/commit/e7a29963724d09b7979c5a27792d27519d182d60))
* truncate bug ([bfc175c](https://github.com/graeter-group/kimmdy/commit/bfc175c819564dd73c0c30b53ed2ce489aa9ca2c))
* type issue in log message formatting ([fdee768](https://github.com/graeter-group/kimmdy/commit/fdee7684a5756ac3f75938ecdf3775703151c834))
* write backups of truncated trajectories to the correct directory ([e615370](https://github.com/graeter-group/kimmdy/commit/e615370a0bec48ead5f66a351f2bc682a0234028))
* wrong type of top when parametrizing ([fb2eb08](https://github.com/graeter-group/kimmdy/commit/fb2eb08af86c12be3a175e3cca2e20819e705da4))


### Documentation

* document kimmdy.yml yaml language server config ([#309](https://github.com/graeter-group/kimmdy/issues/309)) ([2fa494e](https://github.com/graeter-group/kimmdy/commit/2fa494ec3f8bca78ba3e707b0091171d34eb3bf2))
* render ([8637990](https://github.com/graeter-group/kimmdy/commit/86379905c0b18086874cc76e134de3c9b059d245))

## [4.4.1](https://github.com/graeter-group/kimmdy/compare/v4.4.0...v4.4.1) (2023-09-13)


### Bug Fixes

* **test:** increase hypothesis deadline ([f7d8867](https://github.com/graeter-group/kimmdy/commit/f7d8867560600c04099da99a28be12311cb3b0af))

## [4.4.0](https://github.com/graeter-group/kimmdy/compare/v4.3.0...v4.4.0) (2023-09-13)


### Features

* **dummy:** this is not a feature, just a release trigger ([6aa1190](https://github.com/graeter-group/kimmdy/commit/6aa1190fb0eca93388a8e9d7e95c8b2f5945a55d))


### Bug Fixes

* **analysis:** remove y axis label for steps panel in plot_energy ([a2bc88a](https://github.com/graeter-group/kimmdy/commit/a2bc88ac4b8405ae940f0f241fe4778c2039552e))
* **test:** increase hypothesis deadline ([4b9823d](https://github.com/graeter-group/kimmdy/commit/4b9823d45265aecfada2d7bd1f4b90640164d31a))
* **test:** increase hypothesis deadline ([9810e84](https://github.com/graeter-group/kimmdy/commit/9810e8478375e43311e0458d71cc9748ae80d684))

## [4.3.0](https://github.com/graeter-group/kimmdy/compare/v4.2.1...v4.3.0) (2023-09-08)


### Features

* unify output location ([#255](https://github.com/graeter-group/kimmdy/issues/255))! ([cbd4d5d](https://github.com/graeter-group/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))


### Bug Fixes

* closes  [#265](https://github.com/graeter-group/kimmdy/issues/265) , cleanup with flake8 ([d79adbc](https://github.com/graeter-group/kimmdy/commit/d79adbcd798cbe56ded35a96c09679eaf5e4b2c2))
* jobscript and logging ([cbd4d5d](https://github.com/graeter-group/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))
* place_atom ([0e0dc1d](https://github.com/graeter-group/kimmdy/commit/0e0dc1dbf4537f12b5af3cca32aba11c7a1ba4fc))
* properly merge mds general '.*' section and general options with defaults ([cbd4d5d](https://github.com/graeter-group/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))
* set defaults on non-existent sections. fixes [#264](https://github.com/graeter-group/kimmdy/issues/264) ([cbd4d5d](https://github.com/graeter-group/kimmdy/commit/cbd4d5db913b38f1ea02ce37d6ccfae0c50053e5))

## [4.2.1](https://github.com/graeter-group/kimmdy/compare/v4.2.0...v4.2.1) (2023-09-07)


### Bug Fixes

* apply_recipe naming ([6573886](https://github.com/graeter-group/kimmdy/commit/65738861ccc4c9418f4a1c7da3fd5f56afd33d14))
* handle `Relax` in plot_rates ([3099210](https://github.com/graeter-group/kimmdy/commit/3099210335a38504dced4aa5ceffaeaa23617dd3))

## [4.2.0](https://github.com/graeter-group/kimmdy/compare/v4.1.0...v4.2.0) (2023-09-04)


### Features

* **analysis:** open-plot and open-vmd options ([d306e66](https://github.com/graeter-group/kimmdy/commit/d306e66b893a72280bffbe12e26ab05c8777a8a9))
* plot multiple energy terms ([d306e66](https://github.com/graeter-group/kimmdy/commit/d306e66b893a72280bffbe12e26ab05c8777a8a9))
* review analysis ([#245](https://github.com/graeter-group/kimmdy/issues/245)) ([d306e66](https://github.com/graeter-group/kimmdy/commit/d306e66b893a72280bffbe12e26ab05c8777a8a9))


### Bug Fixes

* [#236](https://github.com/graeter-group/kimmdy/issues/236) rm concat cmd option ([#240](https://github.com/graeter-group/kimmdy/issues/240)) ([675ded5](https://github.com/graeter-group/kimmdy/commit/675ded537a6b241b59b02560c286945e76916f22))
* missing test file ([1e805c1](https://github.com/graeter-group/kimmdy/commit/1e805c19a13d4f590a0c8f8e7841428f64ef1c71))
* no grappa in tox tests ([feaff46](https://github.com/graeter-group/kimmdy/commit/feaff46e0d44316cf5cbf46c8045d9e0507850ed))
* Recipe aggregation broken, closes [#252](https://github.com/graeter-group/kimmdy/issues/252) ([#253](https://github.com/graeter-group/kimmdy/issues/253)) ([a6177e3](https://github.com/graeter-group/kimmdy/commit/a6177e3d99af093700167d03311718c9ce571f63))

## [4.1.0](https://github.com/graeter-group/kimmdy/compare/v4.0.1...v4.1.0) (2023-08-28)


### Features

* --debug for post-mortem debugging ([#214](https://github.com/graeter-group/kimmdy/issues/214)) ([91360d1](https://github.com/graeter-group/kimmdy/commit/91360d156a86094bb46111b8a3858404bc3d8daf))
* add remove_hydrogen tool ([#213](https://github.com/graeter-group/kimmdy/issues/213)) ([b7f7e78](https://github.com/graeter-group/kimmdy/commit/b7f7e78273458cc89532f5db8c17e58f1d6c6e50))
* Automatic Radical Parameterization ([#144](https://github.com/graeter-group/kimmdy/issues/144)) ([4b4641b](https://github.com/graeter-group/kimmdy/commit/4b4641bce00812be0859acd8200fce658577a316))
* automatically stop and restart KIMMDY for joblength-restricted HPC clusters ([#162](https://github.com/graeter-group/kimmdy/issues/162)) ([602d864](https://github.com/graeter-group/kimmdy/commit/602d86483aff2f6bf039f4bd74d018e63cd9b7c4))
* cleanup cmd.py ([#228](https://github.com/graeter-group/kimmdy/issues/228)) ([f77b6dc](https://github.com/graeter-group/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* config: use_plumed in md section ([#234](https://github.com/graeter-group/kimmdy/issues/234)) ([22dfacb](https://github.com/graeter-group/kimmdy/commit/22dfacbe38bf167ab79b391b40ba1883cf090881))
* multimolecule topology ([#211](https://github.com/graeter-group/kimmdy/issues/211)) ([153af8c](https://github.com/graeter-group/kimmdy/commit/153af8ce6780dbff928160c89ff17a5a3c373969))
* single reaction tasks in config ([#215](https://github.com/graeter-group/kimmdy/issues/215)) ([1c07689](https://github.com/graeter-group/kimmdy/commit/1c076897f429673dbe0d7dfec9caaf86a7018a71))
* workflow for local testing ([092b89e](https://github.com/graeter-group/kimmdy/commit/092b89e5fcb1addcd2a3d09f9104c0ddc2b3c3eb))


### Bug Fixes

* 212 ([016aa56](https://github.com/graeter-group/kimmdy/commit/016aa56a630f64b953cb6cab7c9a25584c47a692))
* allow reading topology with an empty forcefield ([#238](https://github.com/graeter-group/kimmdy/issues/238)) ([a0824c3](https://github.com/graeter-group/kimmdy/commit/a0824c37379f0ebdf87bbf73dd50a0b114bbeafe))
* also fix [#200](https://github.com/graeter-group/kimmdy/issues/200) ([f2b2be2](https://github.com/graeter-group/kimmdy/commit/f2b2be2169bc7f91eeb925cf0e125a321733cb1c))
* catch plugin loading exception ([268f77f](https://github.com/graeter-group/kimmdy/commit/268f77fdda7abee65d8347f6a14b8cfb5be57274))
* **ci:** render docs ([a62f351](https://github.com/graeter-group/kimmdy/commit/a62f35162bbbd4b10f9f7f39138a52734cff65c7))
* consolidate cmd interface with previous changes by [@ehhartmann](https://github.com/ehhartmann) ([f77b6dc](https://github.com/graeter-group/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* **docs:** getting started molstar resources ([f77b6dc](https://github.com/graeter-group/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* don't increment logfile if starting from a checkpoint ([f77b6dc](https://github.com/graeter-group/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* don't run tasks twice ([c756793](https://github.com/graeter-group/kimmdy/commit/c7567932bce51ded317d689e2bec7d74e2dbae4d))
* examples use new plumed field ([9493a45](https://github.com/graeter-group/kimmdy/commit/9493a4582b57ce6cba065aaa11930a701d5d494c))
* fix docstring references after https://github.com/graeter-group/kimmdy/pull/229 ([7c4ba75](https://github.com/graeter-group/kimmdy/commit/7c4ba75e2e72430f0691a4fa9019a87a2f507ee4))
* install order and tox rquirements ([e7d5957](https://github.com/graeter-group/kimmdy/commit/e7d5957ca1bf757298c6a4eff5f9b4160f3290ef))
* logger instead of logging ([272dd41](https://github.com/graeter-group/kimmdy/commit/272dd411d31ffbd7bb49ed3a1bd31364e9b4ad0f))
* missing import ([0ac234f](https://github.com/graeter-group/kimmdy/commit/0ac234f8ce88279e0d48843dc7c71c51d28fd847))
* morse transition rate homolysis tests fails on ci but not locally [#197](https://github.com/graeter-group/kimmdy/issues/197)  ([#201](https://github.com/graeter-group/kimmdy/issues/201)) ([f2b2be2](https://github.com/graeter-group/kimmdy/commit/f2b2be2169bc7f91eeb925cf0e125a321733cb1c))
* new build example ([4c6b93e](https://github.com/graeter-group/kimmdy/commit/4c6b93e813acc7eedcd45cbc2f86acc2aefb4fb8))
* plugins use task-logger ([#207](https://github.com/graeter-group/kimmdy/issues/207)) ([d5c419d](https://github.com/graeter-group/kimmdy/commit/d5c419d1dfe741958c0052465ad83c28b58ded65))
* proper usage of __str__ and __repr__ for topology ([#210](https://github.com/graeter-group/kimmdy/issues/210)) ([6fd070f](https://github.com/graeter-group/kimmdy/commit/6fd070f6ec2dcafcb0d5329a7ff60a431eee743c))
* removed need for strict bind break order ([#231](https://github.com/graeter-group/kimmdy/issues/231)) ([a7859ec](https://github.com/graeter-group/kimmdy/commit/a7859ec1100d8bd89eddbfa4816e2be0ada358e7))
* requirements typo ([2f998d2](https://github.com/graeter-group/kimmdy/commit/2f998d232b62dbc802ba871fcb91efdb4b72f507))
* runmng logger getting overwritten ([27baee9](https://github.com/graeter-group/kimmdy/commit/27baee9fd3a716930e99f3e8cb9a677e44e95047))
* show plugins ([dcb92c8](https://github.com/graeter-group/kimmdy/commit/dcb92c8737b71efaaa8429cf2343d53000bef39c))


### Documentation

* improve layout and linking ([f77b6dc](https://github.com/graeter-group/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))
* more expressive docstrings ([f77b6dc](https://github.com/graeter-group/kimmdy/commit/f77b6dcf00e77800bd8899dec6c3f40652f200b9))

## [4.0.1](https://github.com/graeter-group/kimmdy/compare/v4.0.0...v4.0.1) (2023-08-11)


### Documentation

* **ci:** add docs build workflow ([#190](https://github.com/graeter-group/kimmdy/issues/190)) ([#198](https://github.com/graeter-group/kimmdy/issues/198)) ([5fdec0b](https://github.com/graeter-group/kimmdy/commit/5fdec0b4903a578abd2fe75cfadfc9b842168be1))
* **dev:** add information for developers about conventional commits ([5fdec0b](https://github.com/graeter-group/kimmdy/commit/5fdec0b4903a578abd2fe75cfadfc9b842168be1))

## [4.0.0](https://github.com/graeter-group/kimmdy/compare/v3.6.0...v4.0.0) (2023-08-11)


### ⚠ BREAKING CHANGES

* **ci:** does this work?

### Features

* **ci:** update encode to support unicode ([d9f16ad](https://github.com/graeter-group/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))
* **ci:** what if this doesn't have a body? ([d9f16ad](https://github.com/graeter-group/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))
* this is a release-please test ([d9f16ad](https://github.com/graeter-group/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))


### Bug Fixes

* **ci:** or this? ([d9f16ad](https://github.com/graeter-group/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))
* **ci:** release please adds all the things ([d9f16ad](https://github.com/graeter-group/kimmdy/commit/d9f16adb7470479a9e700ffd86597784caa2afe7))

## [3.6.0](https://github.com/graeter-group/kimmdy/compare/v3.5.1...v3.6.0) (2023-08-11)


### Features

* Analysis Module ([#176](https://github.com/graeter-group/kimmdy/issues/176)) ([f37d6c8](https://github.com/graeter-group/kimmdy/commit/f37d6c8fecc97c4f256c4992d63ceaf68633629c))
* **ci:** run test on PR with testthis label ([072b4bf](https://github.com/graeter-group/kimmdy/commit/072b4bf2d9b8a0d8f9bbddede28701d59b79f6bf))


### Bug Fixes

* [#185](https://github.com/graeter-group/kimmdy/issues/185) find gromacs data dir in docker container ([#191](https://github.com/graeter-group/kimmdy/issues/191)) ([7c0b90c](https://github.com/graeter-group/kimmdy/commit/7c0b90cec39c625f780024ec0e766481113a274d))
* **ci:** use correct event context variable ([33003a9](https://github.com/graeter-group/kimmdy/commit/33003a93be90470306410071f42a5b3cd880afe0))
* don't run tests on pr twice ([7957269](https://github.com/graeter-group/kimmdy/commit/79572694e360db4e528a122afecfd2f61b106f6d))
* fix [#187](https://github.com/graeter-group/kimmdy/issues/187) ([587c426](https://github.com/graeter-group/kimmdy/commit/587c426621cb9b403bed23ff008f5e3039e628f4))
* merge issues ([5880dbb](https://github.com/graeter-group/kimmdy/commit/5880dbb67ae029a827dca77222b84a0c9e653cf1))
* plumed_mod updating ([3f815a2](https://github.com/graeter-group/kimmdy/commit/3f815a21b00aa0b7e3556ad9159c8f84e0b95b2b))
* register pytest mark properly ([6c56ed9](https://github.com/graeter-group/kimmdy/commit/6c56ed9892a9fe72b83e21b9da5c286289ad8df1))
* remove unused tests ([8632f8b](https://github.com/graeter-group/kimmdy/commit/8632f8b5bf7fcd3f3da8b5090b01f2d115d897a9))
* test releases in separate workflow to get easy link for a badge ([0a3e11c](https://github.com/graeter-group/kimmdy/commit/0a3e11cd6d08d35dfd7b8c781ecf9f07e3614b01))
* tests compatible with local gh act ([71d8204](https://github.com/graeter-group/kimmdy/commit/71d8204c60e9134af21663c8fc6abe7a05a06c8c))
* trigger workflow ([e86dccd](https://github.com/graeter-group/kimmdy/commit/e86dccd13fd06be4f27a4265df1b0d4b0fb7f399))

## [3.5.1](https://github.com/graeter-group/kimmdy/compare/v3.5.0...v3.5.1) (2023-08-08)


### Bug Fixes

* **ci:** testing ([e9b4b19](https://github.com/graeter-group/kimmdy/commit/e9b4b19b5801f342f0e75e5ccd819e5ae55a9473))
* **ci:** trigger another release PR ([9c824a8](https://github.com/graeter-group/kimmdy/commit/9c824a890628c7995fa40cca77c974d8fa35644d))

## [3.5.0](https://github.com/graeter-group/kimmdy/compare/v3.4.0...v3.5.0) (2023-08-08)


### Features

* **ci:** first try! ([f2c686f](https://github.com/graeter-group/kimmdy/commit/f2c686f5336f5e3e5e0c819990c0fb89bec9f647))


### Bug Fixes

* ci ([3bcc81a](https://github.com/graeter-group/kimmdy/commit/3bcc81a3231ddc2e4a6b7b3d8f112228e2c10226))
* **ci:** almost did it first try... ([2850379](https://github.com/graeter-group/kimmdy/commit/2850379ca5c0d3d0e1e755a1e2170caa4bef38b8))

## [3.4.0](https://github.com/graeter-group/kimmdy/compare/v3.3.1...v3.4.0) (2023-08-08)


### Features

* add gmx_mdrun_flags option. fixes [#145](https://github.com/graeter-group/kimmdy/issues/145) ([#161](https://github.com/graeter-group/kimmdy/issues/161)) ([e780402](https://github.com/graeter-group/kimmdy/commit/e780402e33e796d2672fc286a81031154ecd1569))
* **config:** validation, completion and defaults for kimmdy.yml ([#153](https://github.com/graeter-group/kimmdy/issues/153)) ([de7c78f](https://github.com/graeter-group/kimmdy/commit/de7c78f30de9fb5d5075bff239e82888abc0e35d))
* example from tests ([#175](https://github.com/graeter-group/kimmdy/issues/175)) ([a096a45](https://github.com/graeter-group/kimmdy/commit/a096a45af9a3cd0f7974301a307265404b105424))


### Bug Fixes

* consistent use of `edis` in homolysis config ([#170](https://github.com/graeter-group/kimmdy/issues/170)) ([e847707](https://github.com/graeter-group/kimmdy/commit/e8477074810a523f2cde79f49bc3cc103fa869a5))
* don't allow `Move` w/o ix/id_to_move ([10c9231](https://github.com/graeter-group/kimmdy/commit/10c9231f6dd98e2d51613ad94a2c690058bb84b2))
* inconsistent inheritance ([#154](https://github.com/graeter-group/kimmdy/issues/154)) ([2bec671](https://github.com/graeter-group/kimmdy/commit/2bec67167207bf3c38e2985e5bef7ed1f3a60870))
* load plugins only once in __init__ ([#171](https://github.com/graeter-group/kimmdy/issues/171)) ([88b55ee](https://github.com/graeter-group/kimmdy/commit/88b55ee281f0e9c39b5f527aa21eed88c57e5c8a))
* merge dihedrals ([#174](https://github.com/graeter-group/kimmdy/issues/174)) ([37f6be3](https://github.com/graeter-group/kimmdy/commit/37f6be3b322cdba4175817d40fcf489be992a1f8))
* scheme keys not plugin names but entrypoints ([b878b69](https://github.com/graeter-group/kimmdy/commit/b878b6966ead7dee83a21d274d3de0a8e6aa5291))
* test input files now compatible ([eda027a](https://github.com/graeter-group/kimmdy/commit/eda027a297d0d33bef404907ac7090f5b0b29204))
* test_merge_prm_top ([9d333bc](https://github.com/graeter-group/kimmdy/commit/9d333bc37423ae007019bc742b5eba36a70192e3))
* update scheme paths after building examples ([5076ca4](https://github.com/graeter-group/kimmdy/commit/5076ca4c5196a0d1f12e8f8245f6110ed24417cb))


### Documentation

* improve config option descriptions, names and error reporting ([#163](https://github.com/graeter-group/kimmdy/issues/163)) ([3318bdc](https://github.com/graeter-group/kimmdy/commit/3318bdc0f16fff99efb7493e03d6f9a086814325))
* **types:** add custom renderer ([3318bdc](https://github.com/graeter-group/kimmdy/commit/3318bdc0f16fff99efb7493e03d6f9a086814325))

## [3.3.1](https://github.com/graeter-group/kimmdy/compare/v3.3.0...v3.3.1) (2023-07-12)


### Documentation

* add a lot more docstrings to functions, classes, modules ([#149](https://github.com/graeter-group/kimmdy/issues/149)) ([a0d2a43](https://github.com/graeter-group/kimmdy/commit/a0d2a43ca257934222227897e586ba81d7507497))

## [3.3.0](https://github.com/graeter-group/kimmdy/compare/v3.2.0...v3.3.0) (2023-07-05)


### Features

* linked ixs and ids for RecipSteps ([#137](https://github.com/graeter-group/kimmdy/issues/137)) ([8238df9](https://github.com/graeter-group/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))
* **topology:** enable multiple dihedrals (unique by their periodicity) for proper dihedrals ([8238df9](https://github.com/graeter-group/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))


### Bug Fixes

* add default funct for Bonds and BondTypes ([8238df9](https://github.com/graeter-group/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))
* add more nestable sections in topology parsing ([8238df9](https://github.com/graeter-group/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))
* match X as wildcards for atomic types ([8238df9](https://github.com/graeter-group/kimmdy/commit/8238df91ffed28568443cdd20e0ff46fdba328c6))

## [3.2.0](https://github.com/graeter-group/kimmdy/compare/v3.1.0...v3.2.0) (2023-07-03)


### Features

* quartodoc ([f10e531](https://github.com/graeter-group/kimmdy/commit/f10e53152d7a9ddc2655bb3744972d06d64eb2b2))
* quartodoc ([#142](https://github.com/graeter-group/kimmdy/issues/142)) ([f10e531](https://github.com/graeter-group/kimmdy/commit/f10e53152d7a9ddc2655bb3744972d06d64eb2b2))

## [3.1.0](https://github.com/graeter-group/kimmdy/compare/v3.0.1...v3.1.0) (2023-06-30)


### Features

* changemanager [#129](https://github.com/graeter-group/kimmdy/issues/129) ([087d652](https://github.com/graeter-group/kimmdy/commit/087d652dba53857c561c8e028133d5f035e099d3))

## [3.0.1](https://github.com/graeter-group/kimmdy/compare/v3.0.0...v3.0.1) (2023-06-28)


### Bug Fixes

* dihedraltypes sections in ff/topology ([affa546](https://github.com/graeter-group/kimmdy/commit/affa5468ab537bcb37d7ae99ceded4fcf48cbc60))

## [3.0.0](https://github.com/graeter-group/kimmdy/compare/v2.4.0...v3.0.0) (2023-06-26)


### ⚠ BREAKING CHANGES

* put restrictions on the way top files are accepted
* deprecate generate_top_from_atom_list

### Features

* handle else macros ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* handle multiple moleculetype sections ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* parse top files with includes ([#124](https://github.com/graeter-group/kimmdy/issues/124)) ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* resolve includes in gromacs data dir (located relative to gmx executable) if not in cwd ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* stop if no ff found, warn if multiple found ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* topology dictionary helpers to get and set sections ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* udpate top dict with parsed ff parameters ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))


### Bug Fixes

* deprecate generate_top_from_atom_list ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* put restrictions on the way top files are accepted ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* use hat_naive instead of hat_reaction ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* use logging.warning instead of depcrecated warn ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))


### Documentation

* add tutorial section (with markdown support) ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))
* move qmd to preview topology to docs ([bcb6121](https://github.com/graeter-group/kimmdy/commit/bcb6121e73a1e54a832cad45f2dd406ac25a26da))

## [2.4.0](https://github.com/graeter-group/kimmdy/compare/v2.3.0...v2.4.0) (2023-06-22)


### Features

* homolysis module and homolysis tests ([3267642](https://github.com/graeter-group/kimmdy/commit/32676421c3d087ed4e38777117b6ffe334883cf3))
* KMC module testing ([ebb771f](https://github.com/graeter-group/kimmdy/commit/ebb771fbbcd730d70a0c5fc7a879a8e0705d48b8))
* Move coordinates include time ([10bb3a7](https://github.com/graeter-group/kimmdy/commit/10bb3a78ac6a53e7f48ef3cf37113518987cedf3))
* Use proper HAT plugin ([6a4b0eb](https://github.com/graeter-group/kimmdy/commit/6a4b0eb4e13a800545b2ff234a30fe6d773439c1))


### Bug Fixes

* allow restarting from checkpoint from python ([2163173](https://github.com/graeter-group/kimmdy/commit/2163173495ea5e325b2cd3c278df1bb6c100c026))
* correct assignment of charge and mass in ff ([c5f92a6](https://github.com/graeter-group/kimmdy/commit/c5f92a672ba6932d64a1ee206fd27c13762d8fe4))
* homolysis test ([178590d](https://github.com/graeter-group/kimmdy/commit/178590d503ff0d58ad6098cc62e6ef150e197dec))
* initialize radical list ([3196747](https://github.com/graeter-group/kimmdy/commit/3196747bfd5d724e805e5c892c6c74d1964ec86a))
* initialize Topology.radicals ([5e0f335](https://github.com/graeter-group/kimmdy/commit/5e0f335c0f7472c30bf1e154913d6afdfc514b6e))
* remove deprecated imports in tests ([d078109](https://github.com/graeter-group/kimmdy/commit/d078109062a4960984591d283be9de647212a284))
* rename variables in test_kmc ([4df479d](https://github.com/graeter-group/kimmdy/commit/4df479d4a8c1f0c656368322d02606b5f74a5044))
* test_integration_move_top ([e0f71b6](https://github.com/graeter-group/kimmdy/commit/e0f71b6506ecab1c14742568f8c47c7a30a15418))


### Documentation

* **general:** add make preview command ([9806e63](https://github.com/graeter-group/kimmdy/commit/9806e63171cec8d4fc3532bbbab13499d7084089))
* render docs ([f306a6c](https://github.com/graeter-group/kimmdy/commit/f306a6c5a09176703c6fe12d3df8a51646fe5a81))

## [2.3.0](https://github.com/graeter-group/kimmdy/compare/v2.2.0...v2.3.0) (2023-05-11)


### Features

* add MOVE conversion for HAT ([0dbe4b1](https://github.com/graeter-group/kimmdy/commit/0dbe4b1b53272bd9c813550cae875134bc679f55))
* make periodicity part of the id for dihedrals ([0710086](https://github.com/graeter-group/kimmdy/commit/0710086473a81878d86ed8ca4daa2b31b921e868))
* parse dihedral and position restraints ([1209714](https://github.com/graeter-group/kimmdy/commit/12097148ff8fe64b4f40eb0d19ad30def6070c7f))
* test for radicals based on bondorder in top initialization ([a446a70](https://github.com/graeter-group/kimmdy/commit/a446a70348bc575b5e5212a23301c51cd9c41233))


### Bug Fixes

* more robust fix for [#93](https://github.com/graeter-group/kimmdy/issues/93) ([3cea313](https://github.com/graeter-group/kimmdy/commit/3cea31394cb7c4abb02439796f724642e4b42e2b))
* rtp parsing for aminoacids.rtp ([025e6c1](https://github.com/graeter-group/kimmdy/commit/025e6c14fc50982db7809079b8418644e3d474df))


### Documentation

* Updates reaction api documentation ([87007b4](https://github.com/graeter-group/kimmdy/commit/87007b485c23262451fe86cb1fb8a9ac0a90fce6)), closes [#114](https://github.com/graeter-group/kimmdy/issues/114)

## [2.2.0](https://github.com/graeter-group/kimmdy/compare/v2.1.0...v2.2.0) (2023-03-07)


### Features

* add --version command ([bf3e183](https://github.com/graeter-group/kimmdy/commit/bf3e183ee57c7b78b3913070bde861a198adc441))


### Bug Fixes

* update importlib syntax and and importlib dependency ([abc89e0](https://github.com/graeter-group/kimmdy/commit/abc89e00cbdc8e4bc5a235bd4a5f8e387deb866e)), closes [#93](https://github.com/graeter-group/kimmdy/issues/93)

## [2.1.0](https://github.com/graeter-group/kimmdy/compare/v2.0.0...v2.1.0) (2023-02-27)


### Features

* add kimmdy checkpoint files ([b468378](https://github.com/graeter-group/kimmdy/commit/b468378bf57219e9b55e7a15b164b82d5fe4a7d3))
* test CI ([d3dff98](https://github.com/graeter-group/kimmdy/commit/d3dff988f0810e3035a4b06c20f842f98135a390))
* write task file history to separate file ([bdefba1](https://github.com/graeter-group/kimmdy/commit/bdefba17dc21a6ad20d79031f70f6c120da1ce9e))


### Bug Fixes

* [#84](https://github.com/graeter-group/kimmdy/issues/84) ([1bef009](https://github.com/graeter-group/kimmdy/commit/1bef00994e841c77e5510ad5a5d807ec8b61a76a))
* 23 ([a6eb2b7](https://github.com/graeter-group/kimmdy/commit/a6eb2b7c4c57ad51fd56de0bc95d1c4f4219144c))
* abort kimmdy if gromacs fails ([080ecc0](https://github.com/graeter-group/kimmdy/commit/080ecc098b09feeac4b166e5639e6ccf26982650)), closes [#83](https://github.com/graeter-group/kimmdy/issues/83)
* allow existing directories of starting from checkpoint ([8370713](https://github.com/graeter-group/kimmdy/commit/83707135c42534b31824bac5393d1d24fcb1642d))
* correctly apply atomtype and residue to jumping H ([dbdd21e](https://github.com/graeter-group/kimmdy/commit/dbdd21e629a9cadac40e4573bd9dde5fc453f97e)), closes [#84](https://github.com/graeter-group/kimmdy/issues/84)
* don't create empty task files if no files are changed ([615d11c](https://github.com/graeter-group/kimmdy/commit/615d11cb582da0604d005433ee340eaa208c98e6))
* don't create logfile christmas trees. ([fca231f](https://github.com/graeter-group/kimmdy/commit/fca231f5805ffbfbdf34a043a3cf5b84fc9c7bbf)), closes [#82](https://github.com/graeter-group/kimmdy/issues/82)
* make plumed really optional ([beb09cd](https://github.com/graeter-group/kimmdy/commit/beb09cde189ce3f110537def1bc3450c75e4500f))
* minimization path to mdp ([971453b](https://github.com/graeter-group/kimmdy/commit/971453b8f2fa34f15cbad1d4256a6bd0cf8408ee))
* use output file dir for histfile name ([a069c3c](https://github.com/graeter-group/kimmdy/commit/a069c3cce53987eeb3966373361d3d3054753239))
