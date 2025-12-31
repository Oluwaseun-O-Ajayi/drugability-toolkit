import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from molecular_descriptors import MolecularDescriptorCalculator
from druglikeness_filters import DrugLikenessFilter
from synthetic_accessibility import SyntheticAccessibilityScorer
import pandas as pd

print("\n" + "="*70)
print("EXAMPLE: High-Throughput Virtual Screening")
print("="*70 + "\n")

# Simulate a compound library
compound_library = {
    'Cmp-001': 'CC(=O)Oc1ccccc1C(=O)O',
    'Cmp-002': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    'Cmp-003': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'Cmp-004': 'CCCCCCCCCCCCc1ccccc1',  # Fails (High LogP)
    'Cmp-005': 'C[C@H]1[C@H]([C@H]([C@@H]([C@@H](O1)O)O)O)O',  # Fails (High TPSA)
    'Cmp-006': 'c1ccc2c(c1)ccc1c2ccc2c1ccc1c2cccc1',  # Fails (Large aromatic)
    'Cmp-007': 'CC(C)CC1CCC(CC1)C(C)C',
    'Cmp-008': 'Cc1ccc(cc1)C(C)C(=O)O',
    'Cmp-009': 'COc1ccc(cc1)CCN',
    'Cmp-010': 'c1ccc(cc1)C(=O)O'
}

print(f"Screening {len(compound_library)} compounds...")
print("="*70 + "\n")

# Initialize tools
desc_calc = MolecularDescriptorCalculator()
filter_tool = DrugLikenessFilter()
sa_scorer = SyntheticAccessibilityScorer()

# Screen compounds
results = []
for comp_id, smiles in compound_library.items():
    # Calculate descriptors
    descriptors = desc_calc.calculate_all_descriptors(smiles)
    if not descriptors:
        continue
    
    # Apply filters
    ro5_pass, ro5_details = filter_tool.lipinski_rule_of_five(smiles)
    veber_pass, veber_details = filter_tool.veber_rules(smiles)
    sa_result = sa_scorer.calculate_sa_score(smiles)
    
    # Collect results
    results.append({
        'Compound_ID': comp_id,
        'SMILES': smiles,
        'MW': descriptors['MW'],
        'LogP': descriptors['LogP'],
        'HBD': descriptors['HBD'],
        'HBA': descriptors['HBA'],
        'TPSA': descriptors['TPSA'],
        'QED': descriptors['QED'],
        'SA_Score': sa_result['SA_score'] if sa_result else None,
        'RO5_Pass': 'Yes' if ro5_pass else 'No',
        'Veber_Pass': 'Yes' if veber_pass else 'No',
        'Drug_Like': 'Yes' if (ro5_pass and veber_pass) else 'No'
    })

# Create DataFrame
df = pd.DataFrame(results)

# Display results
print("SCREENING RESULTS:")
print("="*70)
print(df.to_string(index=False))

# Filter passing compounds
print("\n" + "="*70)
print("PASSING COMPOUNDS (Drug-like):")
print("="*70)

passing = df[df['Drug_Like'] == 'Yes'].copy()
passing = passing.sort_values('QED', ascending=False)

if len(passing) > 0:
    print(f"\n{len(passing)} compounds passed all filters:\n")
    print(passing[['Compound_ID', 'MW', 'LogP', 'QED', 'SA_Score']].to_string(index=False))
    
    # Best compound
    best = passing.iloc[0]
    print(f"\n🏆 TOP COMPOUND:")
    print(f"   ID: {best['Compound_ID']}")
    print(f"   SMILES: {best['SMILES']}")
    print(f"   QED: {best['QED']:.3f}")
    print(f"   SA Score: {best['SA_Score']:.2f}")
else:
    print("\n❌ No compounds passed all filters")

# Statistics
print("\n" + "="*70)
print("SCREENING STATISTICS:")
print("="*70)
print(f"\nTotal screened:     {len(df)}")
print(f"Passed RO5:         {len(df[df['RO5_Pass'] == 'Yes'])} ({len(df[df['RO5_Pass'] == 'Yes'])/len(df)*100:.1f}%)")
print(f"Passed Veber:       {len(df[df['Veber_Pass'] == 'Yes'])} ({len(df[df['Veber_Pass'] == 'Yes'])/len(df)*100:.1f}%)")
print(f"Passed All:         {len(passing)} ({len(passing)/len(df)*100:.1f}%)")

print(f"\nProperty Ranges:")
print(f"  MW:     {df['MW'].min():.1f} - {df['MW'].max():.1f} Da")
print(f"  LogP:   {df['LogP'].min():.2f} - {df['LogP'].max():.2f}")
print(f"  QED:    {df['QED'].min():.3f} - {df['QED'].max():.3f}")

# Save results
output_file = 'results/screening_results.csv'
df.to_csv(output_file, index=False)
print(f"\n✅ Results saved to: {output_file}")
print("="*70 + "\n")