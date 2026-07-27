# Drugability Toolkit
**Open-source RDKit-based workflow for interpretable small-molecule developability assessment, including molecular descriptors, drug-likeness filters, synthetic accessibility scoring, rule-based ADMET property flags, and screening workflows.**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMjAwIiBoZWlnaHQ9IjIwMCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8cmVjdCB3aWR0aD0iMjAwIiBoZWlnaHQ9IjIwMCIgZmlsbD0iIzM4MzhmZiIvPgogIDx0ZXh0IHg9IjUwJSIgeT0iNTAlIiBmb250LXNpemU9IjYwIiBmaWxsPSJ3aGl0ZSIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZHk9Ii4zZW0iPnJkPC90ZXh0Pgo8L3N2Zz4K)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


---

## Overview

Drugability Toolkit is a Python-based cheminformatics workflow for assessing physicochemical properties and early-stage developability characteristics of small molecules

The toolkit integrates established descriptor calculations, literature-derived drug-likeness rules, synthetic accessibility scoring, and qualitative ADMET-related assessments into a single reproducible workflow

It is intended for:

- Medicinal chemistry education
- Computational chemistry research
- Preliminary compound triage
- Virtual screening workflows
- Hypothesis generation and exploratory analysis

### Key Features

- Molecular descriptor calculation (MW, LogP, TPSA, HBD/HBA, QED)

- Drug-likeness assessment using established medicinal chemistry rules

- Qualitative ADMET property flags

- Synthetic accessibility assessment

- Property-based lead assessment guidance

- Batch analysis and screening workflows

- Reproducible example datasets and validation scripts

---

### Important Note

Drugability Toolkit is intended for research support, education, and preliminary screening

It does not replace experimental ADMET testing, medicinal chemistry expertise, or pharmacokinetic evaluation

Outputs should be interpreted as computational estimates and rule-based assessments rather than experimentally validated predictions

---

## Installation

### Requirements

- Python 3.8 or higher
- RDKit 2023.3.1+
- Standard scientific Python stack

### Setup

```bash
# Clone the repository
git clone https://github.com/Oluwaseun-O-Ajayi/drugability-toolkit.git
cd drugability-toolkit

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

```
rdkit>=2023.3.1
pandas>=1.5.0
numpy>=1.24.0
matplotlib>=3.6.0
seaborn>=0.12.0
scikit-learn>=1.2.0
scipy>=1.10.0
```

---

## Quick Start

### Basic Usage

```python
from molecular_descriptors import MolecularDescriptorCalculator
from druglikeness_filters import DrugLikenessFilter
from admet_predictor import ADMETPredictor

# Initialize tools
desc_calc = MolecularDescriptorCalculator()
filter_tool = DrugLikenessFilter()
admet_pred = ADMETPredictor()

# Analyze a molecule (Aspirin example)
smiles = 'CC(=O)Oc1ccccc1C(=O)O'

# Calculate descriptors
descriptors = desc_calc.calculate_all_descriptors(smiles)
print(f"MW: {descriptors['MW']:.2f} Da")
print(f"LogP: {descriptors['LogP']:.2f}")
print(f"QED: {descriptors['QED']:.3f}")

# Apply drug-likeness filters
ro5_pass, ro5_details = filter_tool.lipinski_rule_of_five(smiles)
print(f"Lipinski's Rule of Five: {'PASS' if ro5_pass else 'FAIL'}")

# Predict ADMET properties
admet = admet_pred.predict_all(smiles)
print(f"Absorption: {admet['Absorption']['HIA_class']}")
print(f"hERG Risk: {admet['Toxicity']['hERG_liability']}")
```

### Lead Optimization

```python
from lead_optimizer import LeadOptimizer

optimizer = LeadOptimizer()

# Analyze compound and get recommendations
smiles = 'your_compound_smiles_here'
optimizer.print_report(smiles, "Compound A")
```

---

## Modules

### 1. Molecular Descriptors Calculator

Calculate essential physicochemical properties for drug molecules.

**Key Descriptors:**
- Molecular weight (MW)
- Lipophilicity (LogP)
- Hydrogen bond donors/acceptors (HBD/HBA)
- Topological polar surface area (TPSA)
- Rotatable bonds
- Quantitative Estimate of Drug-likeness (QED)
- Fraction Csp3

**Example:**
```python
from molecular_descriptors import MolecularDescriptorCalculator

calc = MolecularDescriptorCalculator()
descriptors = calc.calculate_all_descriptors('CC(=O)Oc1ccccc1C(=O)O')
calc.print_descriptors('CC(=O)Oc1ccccc1C(=O)O', 'Aspirin')
```

**Output:** Comprehensive report with all molecular descriptors and drug-likeness metrics.

---

### 2. Drug-likeness Filters

Apply validated filters from medicinal chemistry literature.

**Implemented Filters:**
- **Lipinski's Rule of Five** (oral bioavailability)
- **Veber Rules** (oral bioavailability)
- **Ghose Filter** (drug-likeness)
- **Egan Rules** (ADMET properties)
- **BBB Penetration** (CNS drugs)
- **Lead-likeness** (early-stage optimization)

**Example:**
```python
from druglikeness_filters import DrugLikenessFilter

filter_tool = DrugLikenessFilter()
filter_tool.print_report('CC(=O)Oc1ccccc1C(=O)O', 'Aspirin')
```

**Lipinski's Rule of Five Criteria:**
- MW ≤ 500 Da
- LogP ≤ 5
- HBD ≤ 5
- HBA ≤ 10
- Maximum 1 violation allowed

---

### 3. ADMET Property Assessment

Provide qualitative ADMET-related assessments using physicochemical descriptors and literature-derived heuristics

**Assessments Include:**

**Absorption:**
- Human intestinal absorption (HIA)
- Caco-2 permeability
- P-glycoprotein substrate

**Distribution:**
- Volume of distribution
- Plasma protein binding
- Blood-brain barrier penetration

**Metabolism:**
- CYP450 substrate prediction
- Metabolic stability

**Excretion:**
- Renal clearance
- Half-life estimation

**Toxicity:**
- hERG cardiac liability
- Hepatotoxicity risk
- Mutagenicity assessment

**Example:**
```python
from admet_predictor import ADMETPredictor

predictor = ADMETPredictor()
predictor.print_report('CC(=O)Oc1ccccc1C(=O)O', 'Aspirin')
```

---

### 4. Synthetic Accessibility Scorer

Estimate synthetic complexity using validated SA score algorithm.

**SA Score Interpretation:**
- 1-3: Easy to synthesize
- 3-5: Moderately easy
- 5-7: Moderately difficult
- 7-10: Very difficult

**Complexity Factors:**
- Molecular complexity (graph theory)
- Ring systems (fused, spiro, bridged)
- Stereochemistry
- Molecular size
- Structural features

**Example:**
```python
from synthetic_accessibility import SyntheticAccessibilityScorer

scorer = SyntheticAccessibilityScorer()
scorer.print_report('CC(=O)Oc1ccccc1C(=O)O', 'Aspirin')
```

---

### 5. Lead Assessment Guidance

Identify descriptor-based property liabilities and provide qualitative guidance for further investigation.

**Analysis Includes:**
- Property violations (MW, LogP, TPSA, etc.)
- Severity assessment (critical, high, medium, low)
- Specific recommendations with rationale
- Descriptor-based recommendations

**Example:**
```python
from lead_optimizer import LeadOptimizer

optimizer = LeadOptimizer()
analysis = optimizer.analyze_compound('your_smiles_here')
optimizer.print_report('your_smiles_here', 'Lead Compound')
```

**Recommendations Cover:**
- Reducing molecular weight
- Adjusting lipophilicity
- Improving permeability
- Enhancing metabolic stability
- Reducing toxicity risk

---

## Example Workflows

### Workflow 1: Screening Compound Library

```python
from molecular_descriptors import MolecularDescriptorCalculator
from druglikeness_filters import DrugLikenessFilter
import pandas as pd

# Initialize tools
calc = MolecularDescriptorCalculator()
filter_tool = DrugLikenessFilter()

# Compound library
compounds = {
    'Cmp-001': 'CC(=O)Oc1ccccc1C(=O)O',
    'Cmp-002': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    # ... more compounds
}

# Screen all compounds
results = []
for comp_id, smiles in compounds.items():
    descriptors = calc.calculate_all_descriptors(smiles)
    ro5_pass, _ = filter_tool.lipinski_rule_of_five(smiles)
    
    if descriptors:
        results.append({
            'ID': comp_id,
            'MW': descriptors['MW'],
            'LogP': descriptors['LogP'],
            'QED': descriptors['QED'],
            'RO5_Pass': ro5_pass
        })

# Create results table
df = pd.DataFrame(results)
passing = df[df['RO5_Pass'] == True].sort_values('QED', ascending=False)
print(f"Passing compounds: {len(passing)}/{len(df)}")
print(passing)
```

### Workflow 2: Lead Optimization Campaign

```python
from lead_optimizer import LeadOptimizer
from synthetic_accessibility import SyntheticAccessibilityScorer

optimizer = LeadOptimizer()
sa_scorer = SyntheticAccessibilityScorer()

# Initial lead
lead_smiles = 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O'

# Analyze
print("Initial Lead Analysis:")
optimizer.print_report(lead_smiles, "Lead-001")

# Check synthetic accessibility
sa_result = sa_scorer.calculate_sa_score(lead_smiles)
print(f"\nSynthetic Accessibility: {sa_result['SA_score']:.2f}/10")

# Get specific modifications
modifications = optimizer.suggest_modifications(lead_smiles)
for mod in modifications:
    print(f"\n{mod['modification']}")
    print(f"Rationale: {mod['rationale']}")
```

### Workflow 3: ADMET Profiling

```python
from admet_predictor import ADMETPredictor
from druglikeness_filters import DrugLikenessFilter

predictor = ADMETPredictor()
filter_tool = DrugLikenessFilter()

# Candidate molecules
candidates = {
    'Candidate A': 'smiles_string_a',
    'Candidate B': 'smiles_string_b',
    'Candidate C': 'smiles_string_c'
}

# Comprehensive profiling
for name, smiles in candidates.items():
    print(f"\n{'='*70}")
    print(f"Profiling: {name}")
    print(f"{'='*70}")
    
    # ADMET
    predictor.print_report(smiles, name)
    
    # Drug-likeness
    filter_tool.print_report(smiles, name)
```

---

## Validation

### Benchmarking Against FDA-Approved Drugs

The repository includes demonstration datasets and validation scripts using representative FDA-approved oral drugs and selected beyond-Rule-of-Five compounds. These examples are intended to verify expected behavior of rule-based filters and demonstrate reproducible workflows

**Run validation:**
```bash
python validation/validate_lipinski.py
python validation/benchmark_fda_drugs.py
```

### Validation Results

**Lipinski's Rule of Five:**
- Tested on 10 known oral drugs
- Tested on 3 "beyond RO5" molecules
- Expected behavior confirmed for the included demonstration dataset

**FDA Drug Benchmark:**
- 10 FDA-approved oral drugs analyzed
- Provides illustrative comparisons among representative FDA-approved oral drugs


**Key Findings:**
- 80-90% of FDA oral drugs pass Lipinski's RO5
- Average QED score: 0.65-0.75
- High bioavailability drugs show QED > 0.6
- SA scores range from 2-6 for marketed drugs

---

## Project Structure

```
drugability-toolkit/
├── molecular_descriptors.py       # Core descriptor calculator
├── druglikeness_filters.py        # RO5, Veber, Ghose, Egan filters
├── admet_predictor.py             # ADMET property prediction
├── synthetic_accessibility.py     # SA score calculation
├── lead_optimizer.py              # Optimization recommendations
├── toxicity_predictor.py          # Toxicity assessment
├── requirements.txt               # Python dependencies
├── CITATION.cff                   # Citation information
├── README.md                      # This file
├── LICENSE                        # MIT License
├── examples/                      # Example scripts
│   ├── example_drugs.py           # FDA drug analysis
│   ├── admet_analysis.py          # ADMET profiling
│   ├── lead_optimization.py       # Lead opt workflow
│   └── batch_screening.py         # High-throughput screening
├── validation/                    # Validation scripts
│   ├── validate_lipinski.py       # RO5 validation
│   └── benchmark_fda_drugs.py     # FDA drug benchmark
├── data/                          # Reference datasets
│   ├── fda_approved_drugs.csv     # FDA drug database
│   └── test_molecules.smi         # Test SMILES
├── results/                       # Analysis outputs
└── references/                    # Literature references
```

---

## Scientific Background

### Theoretical Basis

The toolkit implements established quantitative structure-property relationships (QSPR) and computational models from pharmaceutical sciences literature.

**Key Concepts:**

1. **Drug-likeness:** Physicochemical properties associated with successful drugs
2. **ADMET:** Pharmacokinetic and safety properties critical for drug development
3. **Synthetic Accessibility:** Computational estimation of synthetic difficulty
4. **Lead Optimization:** Systematic improvement of molecular properties

### Algorithm Validation

The repository includes demonstration-scale validation examples using:
- Representative FDA-approved oral drugs
- Example drug-likeness filters
- Demonstration screening workflows
- Literature-derived medicinal chemistry rules


### Limitations

**Important Notes:**
- Predictions are computational estimates, not experimental measurements
- Results should guide, not replace, experimental validation
- Models are based on oral small molecule drugs
- Beyond Rule of Five (bRO5) molecules may not be well-predicted
- Biologics and peptides require specialized tools

---

## For Academic Use

### Publication Quality

This toolkit is designed for publication in scientific journals and protocols:

- All algorithms based on peer-reviewed literature  
- Validated against known experimental data  
- Comprehensive documentation with citations  
- Reproducible examples and benchmarks  
- Open-source for transparency  

### Use in Research

**Appropriate Uses:**
- Virtual screening campaigns
- Lead optimization projects
- Structure-property relationship studies
- Teaching computational drug design
- Hypothesis generation

**Citation Required:**
Please cite this toolkit in publications (see Citation section below)

---

## How to Cite This Work

### Software Citation

```bibtex
@software{drugability_toolkit_2024,
  author = {Ajayi, Oluwaseun O.},
  title = {Drugability Toolkit: Computational Drug Developability Assessment},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/Oluwaseun-O-Ajayi/drugability-toolkit},
  version = {1.0.0}
}
```

### When Published (Future)

After obtaining DOI and publishing protocol:
```bibtex
@article{ajayi2024drugability,
  author = {Ajayi, Oluwaseun O.},
  title = {Computational Assessment of Small Molecule Drug Developability},
  journal = {Current Protocols},
  year = {2024},
  doi = {10.xxxx/xxxxx}
}
```

---

## Methodology

### Molecular Descriptors

Implemented using RDKit molecular descriptor calculations:
- **Molecular Weight:** Standard atomic weight calculation
- **LogP:** Wildman-Crippen method
- **TPSA:** Topological polar surface area (Ertl algorithm)
- **QED:** Quantitative estimate of drug-likeness

### Drug-likeness Filters

Based on statistical analysis of known drugs:
- **Lipinski's RO5:** Analysis of World Drug Index
- **Veber Rules:** Oral bioavailability study
- **Ghose Filter:** Combinatorial library design
- **Egan Rules:** ADMET property analysis

### ADMET Predictions

Qualitative predictions based on physicochemical properties:
- Empirical rules from literature
- Property-based decision trees
- Statistical models from drug databases

### Synthetic Accessibility

Modified implementation based on molecular complexity:
- Graph theory metrics
- Ring complexity analysis
- Stereochemistry penalty
- Fragment contribution scores

---

## Contributing

Contributions are welcome! Areas for enhancement:

**High Priority:**
- Machine learning models for ADMET
- Integration with external APIs (ChEMBL, PubChem)
- 3D structure-based predictions
- Metabolite prediction

**Medium Priority:**
- GUI interface
- Web application deployment
- Database backend
- Parallel processing for large libraries

**Documentation:**
- Tutorial notebooks
- Video tutorials
- Additional examples
- Translation to other languages

### Development Guidelines

1. Follow PEP 8 style guidelines
2. Add docstrings to all functions
3. Include unit tests for new features
4. Cite literature sources for algorithms
5. Validate against experimental data

---

## License

This project is licensed under the MIT License:

```
MIT License

Copyright (c) 2024 Oluwaseun O. Ajayi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
```

---

## Acknowledgments

- **RDKit Community** for the excellent cheminformatics toolkit
- **Medicinal Chemistry Literature** for validated algorithms
- **Open Source Contributors** for scientific Python ecosystem
- **University of Georgia** for research support

---

## Contact

**Oluwaseun O. Ajayi**  
PhD Researcher in Chemistry  
University of Georgia

- **GitHub:** [@Oluwaseun-O-Ajayi](https://github.com/Oluwaseun-O-Ajayi)
- **Academic Email:** oluwaseun.ajayi@uga.edu
- **Personal Email:** seunolanikeajayi@gmail.com

---

## Additional Resources

### Related Tools
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
- [Open Babel](http://openbabel.org/) - Chemical data conversion
- [ChEMBL](https://www.ebi.ac.uk/chembl/) - Medicinal chemistry database
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Chemical information resource

### Learning Resources
- [Medicinal Chemistry MOOC](https://www.coursera.org/learn/medicinal-chemistry)
- [Drug Discovery Resources](https://www.nature.com/subjects/drug-discovery)
- [QSAR/QSPR Methods](https://pubs.acs.org/journal/jcisd8)

### Recommended Reading
See `references/` folder for comprehensive literature list

---

## Future Development

**Planned Features:**
- Deep learning ADMET models
- Reaction prediction integration
- Retrosynthetic analysis
- Protein-ligand docking preparation
- PAINS filter implementation
- Aggregator prediction
- Fluorescence interference detection

**Community Requests:**
Submit feature requests via GitHub Issues

---

## Star This Repository

If you find this toolkit useful for your research, please:
- Star the repository
- Share with colleagues
- Cite in publications
- Report issues
- Suggest improvements

---

**Advancing computational drug design through open science**
