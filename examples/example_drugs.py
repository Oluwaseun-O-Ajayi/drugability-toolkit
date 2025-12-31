import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from molecular_descriptors import MolecularDescriptorCalculator
from druglikeness_filters import DrugLikenessFilter
from admet_predictor import ADMETPredictor
from synthetic_accessibility import SyntheticAccessibilityScorer
from lead_optimizer import LeadOptimizer
import pandas as pd

print("\n" + "="*70)
print("EXAMPLE: Analyzing FDA-Approved Drugs")
print("="*70 + "\n")

# FDA-approved drugs with known good properties
fda_drugs = {
    'Aspirin': {
        'SMILES': 'CC(=O)Oc1ccccc1C(=O)O',
        'indication': 'Analgesic/Antiplatelet',
        'year_approved': 1899
    },
    'Ibuprofen': {
        'SMILES': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'indication': 'NSAID',
        'year_approved': 1974
    },
    'Atorvastatin (Lipitor)': {
        'SMILES': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O',
        'indication': 'Statin',
        'year_approved': 1996
    },
    'Metformin': {
        'SMILES': 'CN(C)C(=N)NC(=N)N',
        'indication': 'Antidiabetic',
        'year_approved': 1995
    },
    'Sildenafil (Viagra)': {
        'SMILES': 'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
        'indication': 'PDE5 Inhibitor',
        'year_approved': 1998
    }
}

# Initialize tools
desc_calc = MolecularDescriptorCalculator()
filter_tool = DrugLikenessFilter()
admet_pred = ADMETPredictor()
sa_scorer = SyntheticAccessibilityScorer()

print("Comprehensive Analysis of FDA-Approved Drugs:\n")
print("="*70)

for drug_name, drug_info in fda_drugs.items():
    smiles = drug_info['SMILES']
    
    print(f"\n{'='*70}")
    print(f"DRUG: {drug_name}")
    print(f"Indication: {drug_info['indication']}")
    print(f"FDA Approval: {drug_info['year_approved']}")
    print(f"{'='*70}\n")
    
    # Molecular descriptors
    print("1. MOLECULAR DESCRIPTORS:")
    print("-" * 70)
    descriptors = desc_calc.calculate_all_descriptors(smiles)
    if descriptors:
        print(f"   MW:        {descriptors['MW']:.2f} Da")
        print(f"   LogP:      {descriptors['LogP']:.2f}")
        print(f"   HBD:       {descriptors['HBD']}")
        print(f"   HBA:       {descriptors['HBA']}")
        print(f"   TPSA:      {descriptors['TPSA']:.2f} Ų")
        print(f"   QED:       {descriptors['QED']:.3f}")
    
    # Drug-likeness
    print("\n2. DRUG-LIKENESS FILTERS:")
    print("-" * 70)
    ro5_pass, ro5_details = filter_tool.lipinski_rule_of_five(smiles)
    veber_pass, veber_details = filter_tool.veber_rules(smiles)
    print(f"   Lipinski RO5:  {'✓ PASS' if ro5_pass else '✗ FAIL'} ({ro5_details['violations']} violations)")
    print(f"   Veber Rules:   {'✓ PASS' if veber_pass else '✗ FAIL'}")
    
    # ADMET
    print("\n3. ADMET PROPERTIES:")
    print("-" * 70)
    admet = admet_pred.predict_all(smiles)
    if 'Absorption' in admet:
        print(f"   Absorption:    {admet['Absorption']['HIA_class']}")
        print(f"   BBB:           {admet['Distribution']['BBB_penetration']}")
        print(f"   hERG Risk:     {admet['Toxicity']['hERG_liability']}")
    
    # Synthetic Accessibility
    print("\n4. SYNTHETIC ACCESSIBILITY:")
    print("-" * 70)
    sa_result = sa_scorer.calculate_sa_score(smiles)
    if sa_result:
        print(f"   SA Score:      {sa_result['SA_score']:.2f}/10")
        print(f"   Assessment:    {sa_result['interpretation']}")
    
    print("\n" + "="*70)

# Summary table
print("\n\n" + "="*70)
print("SUMMARY TABLE:")
print("="*70 + "\n")

summary_data = []
for drug_name, drug_info in fda_drugs.items():
    smiles = drug_info['SMILES']
    descriptors = desc_calc.calculate_all_descriptors(smiles)
    ro5_pass, ro5_details = filter_tool.lipinski_rule_of_five(smiles)
    sa_result = sa_scorer.calculate_sa_score(smiles)
    
    if descriptors and sa_result:
        summary_data.append({
            'Drug': drug_name,
            'MW': f"{descriptors['MW']:.1f}",
            'LogP': f"{descriptors['LogP']:.2f}",
            'QED': f"{descriptors['QED']:.3f}",
            'RO5': 'Pass' if ro5_pass else 'Fail',
            'SA Score': f"{sa_result['SA_score']:.2f}"
        })

df = pd.DataFrame(summary_data)
print(df.to_string(index=False))

print("\n" + "="*70)
print("✅ Analysis complete!")
print("\nKey Observations:")
print("• All FDA-approved drugs show good drug-likeness")
print("• Most pass Lipinski's Rule of Five")
print("• SA scores vary based on molecular complexity")
print("• QED scores correlate with drug success")
print("="*70 + "\n")