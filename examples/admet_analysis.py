import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from admet_predictor import ADMETPredictor
import pandas as pd

print("\n" + "="*70)
print("EXAMPLE: ADMET Profiling of Drug Candidates")
print("="*70 + "\n")

predictor = ADMETPredictor()

# Test molecules with different ADMET profiles
test_molecules = {
    'Good oral bioavailability': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen
    'CNS penetrant': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
    'Poor absorption': 'C[C@H]1[C@H]([C@H]([C@@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O',  # Sugar
    'BBB non-penetrant': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O'  # Atorvastatin
}

print("Detailed ADMET Analysis:\n")

for name, smiles in test_molecules.items():
    print(f"{'='*70}")
    print(f"{name.upper()}")
    print(f"{'='*70}")
    predictor.print_report(smiles, name)

# Comparative analysis
print("\n" + "="*70)
print("COMPARATIVE ADMET SUMMARY:")
print("="*70 + "\n")

comparison_data = []
for name, smiles in test_molecules.items():
    results = predictor.predict_all(smiles)
    
    if 'error' not in results:
        comparison_data.append({
            'Compound': name,
            'HIA': results['Absorption']['HIA_class'],
            'BBB': results['Distribution']['BBB_penetration'],
            'CYP3A4': results['Metabolism']['CYP3A4_substrate'],
            'hERG': results['Toxicity']['hERG_liability']
        })

df = pd.DataFrame(comparison_data)
print(df.to_string(index=False))

print("\n" + "="*70)
print("✅ ADMET profiling complete!")
print("\nConclusions:")
print("• Oral bioavailability correlates with HIA and TPSA")
print("• CNS penetration requires low TPSA and MW")
print("• Large molecules show poor absorption")
print("• Lipophilic molecules have higher hERG risk")
print("="*70 + "\n")