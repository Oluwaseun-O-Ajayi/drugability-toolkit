"""Case study 4: structural alerts as a triage step.

Question: do PAINS and Brenk alerts separate classic assay-interference
substructures from approved drugs?
Input:  data/alert_reference_compounds.csv, data/fda_approved_drugs.csv
Output: examples/case_study_results/case4_structural_alerts.csv
"""

import pandas as pd
from rdkit import Chem
from _case_study_common import DATA, RESULTS
from drugability_toolkit import structural_alerts

sets = [("Reference substructure", "alert_reference_compounds.csv", "Compound"),
        ("Approved drug", "fda_approved_drugs.csv", "Drug_Name")]

rows = []
for label, path, name_col in sets:
    data = pd.read_csv(DATA / path)
    for name, smiles in zip(data[name_col], data["SMILES"]):
        alerts = structural_alerts(Chem.MolFromSmiles(smiles))
        rows.append({"Set": label, "Compound": name,
                     "PAINS": "; ".join(alerts["PAINS"]), "Brenk": "; ".join(alerts["Brenk"]),
                     "n_PAINS": len(alerts["PAINS"]), "n_Brenk": len(alerts["Brenk"])})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case4_structural_alerts.csv", index=False)
print(table.drop(columns=["n_PAINS", "n_Brenk"]).to_string(index=False))
print("\nCompounds with at least one alert:")
print(table.groupby("Set")[["n_PAINS", "n_Brenk"]].agg(lambda s: int((s > 0).sum())).to_string())
