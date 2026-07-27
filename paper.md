---
title: "Drugability Toolkit: an open-source RDKit-based workflow for interpretable small-molecule developability assessment"
tags:
  - Python
  - RDKit
  - cheminformatics
  - drug-likeness
  - ADMET
  - synthetic accessibility
  - drug discovery
authors:
  - name: Oluwaseun O. Ajayi
    affiliation: 1
affiliations:
  - name: University of Georgia
    index: 1
date: 2026
bibliography: paper.bib
---

# Summary

Drugability Toolkit is an open-source Python package for interpretable early-stage small-molecule developability assessment. The toolkit integrates RDKit-based molecular descriptor calculation, literature-derived drug-likeness filters, rule-based ADMET property flags, synthetic accessibility scoring, and example workflows for batch screening and lead assessment. It is designed for computational chemistry education, preliminary compound triage, and reproducible cheminformatics workflows.

# Statement of need

Early-stage small-molecule discovery workflows often require rapid assessment of physicochemical properties, drug-likeness, synthetic accessibility, and potential developability liabilities. Although individual cheminformatics calculations are available through established libraries such as RDKit, students and early-stage researchers may benefit from an integrated workflow that connects these calculations into a transparent and reproducible assessment pipeline. Drugability Toolkit addresses this need by organizing commonly used descriptor, filtering, and interpretive modules into a single Python-based workflow.

# Functionality

The toolkit provides modules for molecular descriptor calculation, literature-derived drug-likeness filters, qualitative ADMET property flags, synthetic accessibility scoring, property-based lead assessment, and batch screening examples. The ADMET module should be interpreted as a rule-based flagging system, not as an experimentally calibrated pharmacokinetic prediction model.

# Demonstration and validation

The repository includes demonstration datasets and validation scripts using representative FDA-approved oral drugs and beyond-Rule-of-Five molecules. These examples verify expected behavior of rule-based filters and demonstrate reproducible compound triage. The outputs should be interpreted as computational estimates and screening guidance, not substitutes for experimental ADMET or pharmacokinetic validation.

# Limitations

Drugability Toolkit does not introduce a new molecular descriptor algorithm, experimentally trained ADMET model, or medicinal chemistry optimization engine. Its ADMET outputs are rule-based flags intended for hypothesis generation. The toolkit is primarily applicable to small-molecule drug-like compounds and may be less appropriate for peptides, biologics, macrocycles, or beyond-Rule-of-Five compounds without additional validation.

# Availability

The software is released under an open-source license and is intended to be installable as a Python package. Source code, examples, validation scripts, and documentation are provided in the project repository.

# Acknowledgements

The toolkit builds on the RDKit cheminformatics ecosystem and established medicinal chemistry literature on drug-likeness, oral bioavailability, QED, and synthetic accessibility.

# References
