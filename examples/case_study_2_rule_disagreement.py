"""Case study 2: where the published rules disagree.

Question: for which approved drugs do Lipinski, Veber, Ghose and Egan give
different verdicts, and which property causes each failure?
Input:  data/fda_approved_drugs.csv
Output: examples/case_study_results/case2_rule_disagreement.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS
from drugability_toolkit import DrugLikenessFilter

drugs = pd.read_csv(DATA / "fda_approved_drugs.csv")
flt = DrugLikenessFilter()
RULES = {"Lipinski_RO5": "Lipinski", "Veber": "Veber", "Ghose": "Ghose", "Egan": "Egan"}

rows = []
for name, smiles in zip(drugs["Drug_Name"], drugs["SMILES"]):
    result = flt.comprehensive_filter(smiles)
    verdicts = {short: result["Overall"][f"{short}_pass"] for short in RULES.values()}
    if len(set(verdicts.values())) == 1:
        continue  # all rules agree
    for key, short in RULES.items():
        if verdicts[short]:
            continue
        failed = [k[:-5] for k, v in result[key].items() if k.endswith("_pass") and v is False]
        rows.append({"Drug": name, "Failed_rule": short, "Failing_properties": ", ".join(failed),
                     "Rules_passed": ", ".join(r for r, ok in verdicts.items() if ok)})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case2_rule_disagreement.csv", index=False)
print(table.to_string(index=False))
