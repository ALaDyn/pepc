# Changelog

All notable changes to PEPC will be documented in this file (starting with
2.0.0).

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2024-04-22

### Added

- Added merging algorithm, enable/disable via preprocess flag (CPPFLAGS) in
  `makefile.frontend` [pepc-breakup]
- Extended `params_template` to tweak the properties of merging and to support
  options of processing I/O [pepc-breakup]
- Added two example `params` files as a guide for tokamak or electrode plate
  simulations [pepc-breakup]

### Changed

- Renamed obsolete frontend main files used for diagnostic purposes (streakline)
  [pepc-breakup]

## [2.1.0-DVH] - 2024-03-26

### Added

- New DVH frontend [pepc-dvh]
- Included hook to check src style and suggest commit message
- CI extended to include fprettify
- CI extended to include CB

### Fixed

- Speedup interaction backend [pepc-v]
- Improved remeshing [pepc-v], (#9)
- Correct types for checkpointing [pepc-v], (#8)
- Correct VTK output [pepc-v], (#4) (#5)
- Updated type use, solved memory corruption [pepc-v], (#1)
- Consistency check improved [pepc-benchmark]

## [2.1.0-Townsend] - 2023-05-09

### Added

- Townsend frontend (`pepc-breakup`) to simulate plasma initiation in tokamaks
via a Townsend avalanche.

## [2.0.0] - 2023-01-24

Version for the initial Zenodo entry of PEPC. Prior versions have not been
released officially though they were available on request.


<!-- vim: set ts=4 sw=4 tw=80 et :-->
