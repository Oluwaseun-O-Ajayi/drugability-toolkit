# Drugability Toolkit

**Interpretable, rule-based developability profiles for small molecules, built on RDKit.**

[![tests](https://github.com/Oluwaseun-O-Ajayi/drugability-toolkit/actions/workflows/tests.yml/badge.svg)](https://github.com/Oluwaseun-O-Ajayi/drugability-toolkit/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21652624.svg)](https://doi.org/10.5281/zenodo.21652624)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Drugability Toolkit turns a SMILES string into a single, readable profile: physicochemical descriptors, the published drug-likeness rules (Lipinski, Veber, Ghose, Egan, lead-likeness), the quantitative estimate of drug-likeness (QED), the Ertl and Schuffenhauer synthetic accessibility score, PAINS and Brenk structural alerts, qualitative ADMET flags and plain-language suggestions for which properties to adjust. Every number is traceable to an RDKit function or a stated threshold.

## Scope and limitations

- Intended for teaching, for a first look at a compound set, and for transparent triage. It does not replace experimental ADMET or pharmacokinetic studies.
- **Descriptors, rule filters, QED, the SA score and structural alerts** use RDKit implementations of published methods.
- **ADMET flags** are qualitative heuristics with thresholds chosen by the author from common medicinal-chemistry rules of thumb. They have not been validated against experimental ADMET data and are labelled as heuristics in every output. Mutagenicity is not assessed.
- Rules derived from oral small-molecule drugs do not transfer to peptides, macrocycles, natural products or other beyond-rule-of-five chemistry.

## Installation

Requires Python 3.9 or later. RDKit is installed from PyPI as a dependency.

```bash
git clone https://github.com/Oluwaseun-O-Ajayi/drugability-toolkit.git
cd drugability-toolkit
pip install .            # or: pip install ".[test]" to include the test dependencies
```

## Quick start

```python
from drugability_toolkit import (ADMETFlags, DrugLikenessFilter, LeadOptimizer,
                                 MolecularDescriptorCalculator, SyntheticAccessibilityScorer)

smiles = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin

descriptors = MolecularDescriptorCalculator().calculate_all_descriptors(smiles)
print(f"MW {descriptors['MW']:.1f}, logP {descriptors['LogP']:.2f}, QED {descriptors['QED']:.2f}")

passes, details = DrugLikenessFilter().lipinski_rule_of_five(smiles)
print("Lipinski:", "pass" if passes else "fail", f"({details['violations']} violations)")

sa = SyntheticAccessibilityScorer().calculate_sa_score(smiles)
print(f"SA score {sa['SA_score']:.2f}: {sa['interpretation']}")

flags = ADMETFlags().predict_all(smiles)
print("BBB flag:", flags["Distribution"]["BBB_penetration"],
      "| Brenk alerts:", flags["Toxicity"]["Brenk_alerts"] or "none")

LeadOptimizer().print_report(smiles, "Aspirin")
```

More complete scripts are in `examples/` (single molecules, batch screening of a small library, ADMET flags and lead assessment). `validation/` contains two scripts that profile approved oral drugs and beyond-rule-of-five compounds.

## Modules

| Module | Contents | Basis |
|---|---|---|
| `molecular_descriptors` | MW, logP, HBD, HBA, TPSA, rotatable bonds, aromatic rings, Fsp3, molar refractivity, QED, Morgan fingerprints, Tanimoto similarity | RDKit implementations |
| `druglikeness_filters` | Lipinski, Veber, Ghose, Egan, lead-likeness, a simple BBB category | Published thresholds (cited in the code) |
| `synthetic_accessibility` | SA score plus descriptive ring, stereo, size, spiro and bridgehead terms | Ertl and Schuffenhauer reference implementation in RDKit Contrib |
| `admet_flags` | Qualitative absorption, distribution, metabolism, excretion and hERG flags; PAINS and Brenk alerts | Heuristics (flags); RDKit FilterCatalog (alerts) |
| `lead_optimizer` | Property issues and suggested directions of change | Rules of thumb (Gleeson 2008; Johnson et al. 2009) |

## Testing

```bash
pip install ".[test]"
python -m pytest
```

The tests check descriptor values against independent RDKit calls, rule outcomes for molecules that sit on either side of each threshold, agreement of the SA score with the RDKit reference implementation, detection of a known PAINS motif, handling of invalid SMILES, and that every example and validation script runs. They run on Linux, macOS and Windows for every push.

## Contributing and support

Bug reports, questions and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). Please use [GitHub Issues](https://github.com/Oluwaseun-O-Ajayi/drugability-toolkit/issues); the author aims to respond within two weeks.

## Citation

If you use the toolkit, please cite the archived software (see `CITATION.cff`):

Ajayi, O. O. (2026). *Drugability Toolkit*. Zenodo. https://doi.org/10.5281/zenodo.21652624

## License

MIT. See [LICENSE](LICENSE).

## Support

Questions, bug reports and feature requests are welcome as GitHub issues. The author maintains the package and reviews issues and pull requests; see CONTRIBUTING.md for how to contribute.
