# Case studies

Five worked examples that show what Drugability Toolkit can be used for. Each one
runs on data bundled with the repository, writes its results to
`examples/case_study_results/`, and can be reproduced with a single command:

```bash
python examples/case_study_1_approved_drugs.py
python examples/case_study_2_rule_disagreement.py
python examples/case_study_3_beyond_ro5.py
python examples/case_study_4_structural_alerts.py
python examples/case_study_5_lead_optimization.py
```

The numbers below were produced with version 1.1.1. Every structure in the input
files was checked against its published molecular formula.

| Case | Application | Input data |
|---|---|---|
| 1 | Profiling a compound set against published rules | `data/fda_approved_drugs.csv` (10 approved oral drugs) |
| 2 | Explaining why rules disagree on the same molecule | `data/fda_approved_drugs.csv` |
| 3 | Recognizing compounds outside the rules' domain | `data/beyond_ro5_drugs.csv` (3 drugs) |
| 4 | Triage with structural alerts | `data/alert_reference_compounds.csv` (5 reference substructures) |
| 5 | Reading lead-optimization advice critically | Atorvastatin and lisinopril |

## Case 1: Profiling approved oral drugs

**Question.** How many marketed oral drugs pass each published drug-likeness rule?

**Result.** Of the ten drugs, 9 pass Lipinski's rule of five, 8 pass Veber, 6 pass
Ghose, 8 pass Egan, and 1 (omeprazole) meets the lead-likeness criteria.
Atorvastatin fails Lipinski with two violations (MW 558.7, LogP 6.31). Sertraline
passes with one violation (LogP 5.18), because the rule allows one.

**What it shows.** The rules describe where most oral drugs lie, not a boundary
that every oral drug respects. Lead-likeness is a stricter criterion meant for
starting points, so most finished drugs fail it.

## Case 2: Where the rules disagree

**Question.** For which drugs do the four rules give different verdicts, and which
property causes each failure?

**Result.**

| Drug | Fails | Failing property | Passes |
|---|---|---|---|
| Aspirin | Ghose | atom count | Lipinski, Veber, Egan |
| Ibuprofen | Ghose | atom count | Lipinski, Veber, Egan |
| Metformin | Ghose | LogP, MW, molar refractivity, atom count | Lipinski, Veber, Egan |
| Lisinopril | Veber | rotatable bonds (12) | Lipinski, Ghose |
| Lisinopril | Egan | TPSA | Lipinski, Ghose |

**What it shows.** Ghose sets lower limits as well as upper ones, so small drugs
fail it. Veber and Egan weigh flexibility and polarity, which Lipinski does not
consider. Because every verdict is traced to a named property and threshold, a
user can see why the rules disagree instead of receiving a single score.

## Case 3: Drugs beyond the rule of five

**Question.** How do large natural-product drugs compare with small oral drugs?

**Result.** Cyclosporine (MW 1202.6), vancomycin (MW 1449.3) and paclitaxel
(MW 853.9) fail the rule of five with 2, 3 and 2 violations. Their synthetic
accessibility scores are 5.92 to 8.21, against 1.58 to 3.67 for the ten small oral
drugs (1 = easy, 10 = very difficult).

**What it shows.** These are approved drugs, so failing the rule of five does not
make a compound useless. Cyclosporine is taken orally despite its size,
paclitaxel is given intravenously, and oral vancomycin is not absorbed, so the
rule's prediction holds for some of them and not for others. All three are
natural products, and their complex ring systems and stereocentres explain the
high SA scores. The input SMILES carry no stereochemistry, so their SA scores are
approximate.

## Case 4: Structural alerts as a triage step

**Question.** Do PAINS and Brenk alerts separate known assay-interference
substructures from approved drugs?

**Result.** All five reference compounds (an ene-rhodanine, catechol, a quinone,
an aminophenol and an azo compound) trigger a PAINS alert. None of the ten
approved drugs triggers a PAINS alert, but five of them trigger at least one Brenk
alert: aspirin (phenol ester), metformin (imine), lisinopril and amlodipine
(aliphatic long chain) and warfarin (coumarin).

**What it shows.** PAINS alerts flag substructures known to interfere with
assays. Brenk alerts are broader and flag groups that may be reactive or
undesirable, so they also fire on safe, marketed drugs. An alert is a reason to
look more closely, not a reason to discard a compound.

## Case 5: Reading optimization advice critically

**Question.** What does the lead-optimization module recommend for two approved
drugs that fail a rule?

**Result.** Atorvastatin receives high priority, with issues on molecular weight
(558.65 > 500), LogP (6.31 > 5) and rotatable bonds (12 > 10), and the advice to
reduce size, lipophilicity and flexibility. Lisinopril receives low priority, with
one issue (12 rotatable bonds). Both also receive a general note that drug-likeness
is low, because their QED values (0.16 and 0.38) are below the 0.5 cut-off.

**What it shows.** Applied to atorvastatin, one of the most prescribed drugs, the
module recommends redesign. Its advice follows from rule thresholds and does not
account for potency, transporters or formulation. It is a starting point for
discussion, not a design decision.
