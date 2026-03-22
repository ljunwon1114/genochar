# Changelog

## [0.6.3.2]

### Added
- Added `genochar setup` for automatic environment setup
- Added support for isolated Prokka and CheckM2 environments
- Added automatic CheckM2 database download

### Changed
- CheckM2 is now executed before annotation when `--check` is enabled
- Simplified annotation support (Prokka only, Bakta removed)
- Improved CLI to single-command interface

### Fixed
- Fixed CheckM2 input handling for assembly directories
- Fixed workflow order and dependency handling

---

## [0.6.3]

### Added
- Initial implementation of `setup` command
- Managed environment support for external tools

---

## [0.6.2]

### Changed
- Removed Bakta from built-in annotation workflow
- Standardized Python environment recommendation

---

## [0.6.1]

### Changed
- Renamed project branding to GenoChar

---

## [0.6.0]

### Changed
- Unified CLI (removed separate summarize/pipeline commands)

---

## [0.5.0]

### Changed
- Refactored CLI options and output schema

---

## [0.4.0]

### Added
- First fully documented release
- FASTA, GFF, and CheckM2 integration