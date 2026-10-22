"""Case study 5: what the optimization module recommends, and its limits.

Question: what does the lead-optimization module flag for two approved drugs
that fail a rule, and how should a user read that advice?
Input:  atorvastatin and lisinopril from data/fda_approved_drugs.csv
Output: examples/case_study_results/case5_lead_optimization.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS
from drugability_toolkit import LeadOptimizer

drugs = pd.read_csv(DATA / "fda_approved_drugs.csv").set_index("Drug_Name")
opt = LeadOptimizer()

rows = []
for name in ["Atorvastatin", "Lisinopril"]:
    result = opt.analyze_compound(drugs.loc[name, "SMILES"])
    for issue in result["issues"]:
        rows.append({"Drug": name, "Priority": result["priority"], "Property": issue["property"],
                     "Value": round(issue["value"], 2), "Threshold": issue["threshold"],
                     "Severity": issue["severity"]})
    print(f"\n{name}: priority = {result['priority']}")
    print("  " + result["overall_assessment"])
    for rec in result["recommendations"]:
        print("  " + rec.replace("•", "-"))

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case5_lead_optimization.csv", index=False)
print("\n" + table.to_string(index=False))
