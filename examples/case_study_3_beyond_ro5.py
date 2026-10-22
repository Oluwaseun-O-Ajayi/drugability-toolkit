"""Case study 3: approved drugs beyond the rule of five.

Question: how do large natural-product drugs compare with small oral drugs on
the rule of five and on synthetic accessibility?
Input:  data/beyond_ro5_drugs.csv, data/fda_approved_drugs.csv
Output: examples/case_study_results/case3_beyond_ro5.csv
Note:   the beyond-rule-of-five SMILES carry no stereochemistry; the SA score
        includes a stereocentre term, so values for these three are approximate.
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, sa_value
from drugability_toolkit import DrugLikenessFilter, SyntheticAccessibilityScorer

flt = DrugLikenessFilter()
sa = SyntheticAccessibilityScorer()

rows = []
for label, path in [("Beyond Ro5", "beyond_ro5_drugs.csv"), ("Small oral drug", "fda_approved_drugs.csv")]:
    data = pd.read_csv(DATA / path)
    for name, smiles in zip(data["Drug_Name"], data["SMILES"]):
        passes, ro5 = flt.lipinski_rule_of_five(smiles)
        rows.append({"Set": label, "Drug": name, "MW": round(ro5["MW"], 1), "HBD": ro5["HBD"],
                     "HBA": ro5["HBA"], "Ro5_violations": ro5["violations"], "Lipinski": passes,
                     "SA_score": sa_value(sa.calculate_sa_score(smiles))})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case3_beyond_ro5.csv", index=False)
print(table.to_string(index=False))
print("\nSA score range by set:")
print(table.groupby("Set")["SA_score"].agg(["min", "max"]).to_string())
