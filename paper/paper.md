---
title: 'Drugability Toolkit: interpretable, rule-based developability profiles for small molecules'
tags:
  - Python
  - cheminformatics
  - drug-likeness
  - medicinal chemistry
  - RDKit
authors:
  - name: Oluwaseun O. Ajayi
    orcid: 0000-0003-0040-7217
    affiliation: 1
affiliations:
  - name: Department of Chemistry, University of Georgia, Athens, Georgia, United States
    index: 1
date: 15 October 2026
bibliography: paper.bib
---

# Summary

Before a molecule is made or bought for testing, chemists usually check a handful of properties that predict whether it could become an oral drug: its size, how fat- or water-loving it is, how many hydrogen bonds it can make, how hard it would be to synthesize, and whether it contains groups that often cause trouble in assays. Drugability Toolkit is a Python package that computes these properties from a molecule's structure, written as a SMILES string, and returns them as one readable profile. It combines descriptors, the published drug-likeness rules, the quantitative estimate of drug-likeness, a synthetic accessibility score, structural alerts, qualitative absorption, distribution, metabolism, excretion and toxicity (ADMET) flags, and plain-language suggestions about which properties to change.

# Statement of need

The individual calculations are available in RDKit [@RDKit], but using them correctly requires knowing which function implements which published rule, which thresholds apply, and how to read the result. Students and early-career researchers often reassemble these steps in notebooks, which makes results hard to compare and easy to get wrong. Web services such as SwissADME [@Daina2017] provide polished reports but are not scriptable inside a local pipeline and do not expose their logic as code.

Drugability Toolkit is designed for three uses: teaching medicinal-chemistry rules with code whose every threshold can be read and cited; a first, transparent triage of small compound sets in academic projects; and a reproducible profiling step that can be called from larger workflows. It is not intended to predict ADMET outcomes, and every heuristic output is labelled as such.

# State of the field

RDKit [@RDKit] supplies the underlying descriptors, the quantitative estimate of drug-likeness (QED) [@Bickerton2012], the reference implementation of the synthetic accessibility (SA) score [@Ertl2009] and catalogs of structural alerts. The medchem library [@medchem] offers a broad, well-tested collection of rule-based and structural filters aimed at filtering large libraries. SwissADME [@Daina2017] combines rules with predictive models in a web interface.

Drugability Toolkit was written as a separate package rather than as a contribution to medchem because its unit of output is different: one interpretable report per molecule that places descriptors, rule outcomes, the SA score, alerts and suggested changes side by side, with the reasons for each verdict. It is small and dependency-light (RDKit, pandas, NumPy and Matplotlib), so its logic can be read in full.

# Software design

The package has five modules with a common interface: each takes a SMILES string and returns a dictionary, and each returns an explicit error for an invalid structure.

- `molecular_descriptors` wraps RDKit descriptors, QED and Morgan fingerprints with Tanimoto similarity.
- `druglikeness_filters` applies the Lipinski [@Lipinski2001], Veber [@Veber2002], Ghose [@Ghose1999] and Egan [@Egan2000] rules and lead-likeness criteria, reporting each threshold and whether it was met.
- `synthetic_accessibility` returns the Ertl and Schuffenhauer SA score [@Ertl2009] from the RDKit reference implementation, together with descriptive ring, stereocentre, size, spiro and bridgehead counts that are reported separately from the score.
- `admet_flags` returns PAINS [@Baell2010] and Brenk [@Brenk2008] structural alerts from the RDKit catalogs and qualitative ADMET flags. The flags use author-chosen thresholds on a few descriptors; they have not been validated against experimental data, and each output carries that statement. Mutagenicity is reported as not assessed rather than estimated.
- `lead_optimizer` combines the other modules into a report that lists property issues and the direction in which each property would need to move.

The main design choice is transparency over predictive power: no fitted model is shipped, so every verdict can be traced to a published rule or a stated threshold. The cost is that the ADMET flags are coarse; a user who needs predictions should use a validated model.

Correctness is checked by a test suite that compares descriptors with direct RDKit calls, tests each rule on molecules on either side of its thresholds, confirms that the SA score matches the reference implementation, detects a known PAINS chemotype, and runs every example and validation script. Continuous integration runs the tests on Linux, macOS and Windows. Version 1.1 corrected an earlier SA module that described a home-made complexity score as the published method; the change is documented in the changelog.

# Research impact statement

The toolkit gives instructors and students a citable, tested reference for how the common drug-likeness rules are computed, and gives small academic projects a scriptable profiling step whose assumptions are visible. The repository includes worked examples on approved oral drugs and beyond-rule-of-five compounds. Benchmarking the ADMET flags against experimental data is planned and is required before they are used for any compound decision.

# AI usage disclosure

The author designed and wrote the original version of this software in December 2025, and wrote this paper. In October 2026, Claude Opus 5.5 (Anthropic) was used for code review, refactoring, and testing. The author directed and verified all changes and takes full responsibility for the work.

# Acknowledgements

This work was conducted independently and received no specific funding from any agency in the public, commercial, or not-for-profit sectors.

# References
