# Changelog
All notable changes to this project will be documented in this file.

[Unreleased changes]

## Version [0.4.2] - 2018-11-05

#### Changed
- minor fixes in README and documentation

## Version [0.4.1] - 2018-11-02
#### Added
- In code documentation Vegetation Prior

#### Changed
- big update on general documentation
- config file is read from `package_ressources`
- prior .vrt files are now always global



## Version [0.4] - 2018-09-01
#### Added
- command line interface to allow user to add prior data
- first implementation of coarse resolution soil moisture prior based on SMAP L4 data
- averaging and aggregation of output if multiple rasters are available for one date or variable
- logging in prior engine

#### Changed
- prior engine framework
  - sub-engine from entry points in `setup.py`
  - conventions through abstract base class implementation in prior creator
- in-code documentation
- fixed travis installation

#### Removed
- -

## Version [0.3] - 2018-03-07
#### Added
- *get\_mean\_state\_vector* returns path to prior files and routes to specific submodule for soil and vegetation related priors respectively to produce/provide information.
- Vegetation prior:
  - global vegetation trait maps as static prior implemented
- Soil moisture prior:
  - basic implementation of ESA CCI soil moisture climatology based prior

#### Changed
- -

#### Removed
- -

---
*The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).*


[Unreleased changes]: https://github.com/multiply-org/prior-engine/compare/v0.4.2...HEAD
[0.4.2]: https://github.com/multiply-org/prior-engine/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/multiply-org/prior-engine/compare/v0.4...v0.4.1
[0.4]: https://github.com/multiply-org/prior-engine/compare/v0.3...v0.4
[0.3]: https://github.com/multiply-org/prior-engine/compare/c76e059...v0.3
