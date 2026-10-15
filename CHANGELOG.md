# Changelog

## 1.1.0 (2026-10)

### Fixed
- The synthetic accessibility score now uses the Ertl and Schuffenhauer reference implementation distributed with RDKit. Version 1.0 computed a different, author-defined complexity score while describing it as the published method; the two correlated poorly (Spearman 0.14 on the 1,128 compounds of the ESOL data set). The former complexity terms are still reported, labelled as descriptive and not part of the score.
- `validation/benchmark_fda_drugs.py` no longer fails with a syntax error.
- The SMILES labelled propranolol in `validation/validate_lipinski.py` encoded pindolol; corrected.
- `data/test_molecules.smi` was empty; it now holds nine molecules whose SMILES were checked against their molecular formulas (the mescaline, limonene and anthracene entries in the earlier draft were incorrect).
- Example and validation scripts now write their outputs next to themselves and run from any working directory.

### Changed
- Code moved into an installable `drugability_toolkit` package (`pip install .`).
- `admet_predictor` renamed `admet_flags` (`ADMETFlags`; `ADMETPredictor` remains as an alias). Outputs are qualitative categories labelled as heuristics; invented numeric values (for example, a fixed "HIA score" of 90, 70 or 30) and references that did not describe the implemented rules were removed. Mutagenicity, previously always reported as "Low Risk", is now reported as "Not assessed".
- Unsourced bioavailability and approval-year labels removed from the bundled drug data, together with the conclusion drawn from them.
- Unused dependencies (scikit-learn, SciPy) removed.

### Added
- PAINS and Brenk structural alerts from the RDKit FilterCatalog.
- Test suite and continuous integration on Linux, macOS and Windows; `CONTRIBUTING.md`.


## Version 1.0.0

### Features

- Molecular descriptor calculation
- Drug-likeness filters
- ADMET property assessment
- Synthetic accessibility scoring
- Lead assessment guidance
- Batch screening workflows

### Validation

- Demonstration-scale FDA-approved drug analysis
- Lipinski Rule of Five validation examples
- Beyond Rule of Five demonstration examples

### Documentation

- Added publication manuscript
- Added citation metadata
- Added repository documentation improvements
