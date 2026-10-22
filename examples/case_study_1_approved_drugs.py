"""Case study 1: profile ten approved oral drugs.

Question: how many marketed oral drugs pass each published drug-likeness rule?
Input:  data/fda_approved_drugs.csv
Output: examples/case_study_results/case1_approved_drugs.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS
from drugability_toolkit import DrugLikenessFilter

drugs = pd.read_csv(DATA / "fda_approved_drugs.csv")
flt = DrugLikenessFilter()

rows = []
for name, smiles in zip(drugs["Drug_Name"], drugs["SMILES"]):
    overall = flt.comprehensive_filter(smiles)["Overall"]
    _, ro5 = flt.lipinski_rule_of_five(smiles)
    rows.append({
        "Drug": name,
        "MW": round(ro5["MW"], 1),
        "LogP": round(ro5["LogP"], 2),
        "Ro5_violations": ro5["violations"],
        "Lipinski": overall["Lipinski_pass"],
        "Veber": overall["Veber_pass"],
        "Ghose": overall["Ghose_pass"],
        "Egan": overall["Egan_pass"],
        "Lead_like": overall["Lead_like"],
    })

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case1_approved_drugs.csv", index=False)
print(table.to_string(index=False))
print("\nPass counts out of", len(table))
for rule in ["Lipinski", "Veber", "Ghose", "Egan", "Lead_like"]:
    print(f"  {rule:10} {int(table[rule].sum())}")
