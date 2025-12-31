import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from druglikeness_filters import DrugLikenessFilter
import pandas as pd

print("\n" + "="*70)
print("VALIDATION: Lipinski's Rule of Five")
print("Testing Against Known Drug Molecules")
print("="*70 + "\n")

filter_tool = DrugLikenessFilter()

# Validation dataset: Known drugs that PASS RO5
passing_drugs = {
    'Aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    'Acetaminophen': 'CC(=O)Nc1ccc(O)cc1',
    'Metformin': 'CN(C)C(=N)NC(=N)N',
    'Caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'Sildenafil': 'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
    'Atenolol': 'CC(C)NCC(COc1ccc(CC(N)=O)cc1)O',
    'Amoxicillin': 'CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O',
    'Warfarin': 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O',
    'Propranolol': 'CC(C)NCC(O)COc1cccc2[nH]ccc12'
}

# Known molecules that FAIL RO5 (but may still be drugs - beyond RO5)
failing_molecules = {
    'Cyclosporine': 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C',
    'Vancomycin': 'CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)O)Cl)CO)O)O)(C)N)O',
    'Paclitaxel': 'CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C'
}

print("Testing Known Drug Molecules (Expected to PASS):")
print("="*70 + "\n")

pass_results = []
for drug, smiles in passing_drugs.items():
    passes, details = filter_tool.lipinski_rule_of_five(smiles)
    pass_results.append({
        'Drug': drug,
        'Pass': passes,
        'Violations': details['violations'],
        'MW': f"{details['MW']:.1f}",
        'LogP': f"{details['LogP']:.2f}",
        'HBD': details['HBD'],
        'HBA': details['HBA']
    })
    status = "✓ PASS" if passes else "✗ FAIL"
    print(f"{drug:20} {status} (violations: {details['violations']})")

print("\n" + "="*70)
print("Testing Molecules Beyond RO5 (Expected to FAIL):")
print("="*70 + "\n")

fail_results = []
for mol, smiles in failing_molecules.items():
    passes, details = filter_tool.lipinski_rule_of_five(smiles)
    fail_results.append({
        'Molecule': mol,
        'Pass': passes,
        'Violations': details['violations'],
        'MW': f"{details['MW']:.1f}",
        'LogP': f"{details['LogP']:.2f}",
        'HBD': details['HBD'],
        'HBA': details['HBA']
    })
    status = "✓ PASS" if passes else "✗ FAIL"
    print(f"{mol:20} {status} (violations: {details['violations']})")

# Validation summary
print("\n" + "="*70)
print("VALIDATION SUMMARY:")
print("="*70 + "\n")

expected_pass = len(passing_drugs)
actual_pass = sum(1 for r in pass_results if r['Pass'])
pass_accuracy = (actual_pass / expected_pass) * 100

expected_fail = len(failing_molecules)
actual_fail = sum(1 for r in fail_results if not r['Pass'])
fail_accuracy = (actual_fail / expected_fail) * 100

print(f"Known Drugs (should pass):")
print(f"  Expected: {expected_pass}")
print(f"  Actual:   {actual_pass}")
print(f"  Accuracy: {pass_accuracy:.1f}%\n")

print(f"Beyond RO5 Molecules (should fail):")
print(f"  Expected: {expected_fail}")
print(f"  Actual:   {actual_fail}")
print(f"  Accuracy: {fail_accuracy:.1f}%\n")

overall_accuracy = ((actual_pass + actual_fail) / (expected_pass + expected_fail)) * 100
print(f"Overall Accuracy: {overall_accuracy:.1f}%")

# Create detailed tables
print("\n" + "="*70)
print("DETAILED RESULTS - DRUG MOLECULES:")
print("="*70 + "\n")

df_pass = pd.DataFrame(pass_results)
print(df_pass.to_string(index=False))

print("\n" + "="*70)
print("DETAILED RESULTS - BEYOND RO5 MOLECULES:")
print("="*70 + "\n")

df_fail = pd.DataFrame(fail_results)
print(df_fail.to_string(index=False))

# Save results
df_pass.to_csv('validation/lipinski_validation_drugs.csv', index=False)
df_fail.to_csv('validation/lipinski_validation_beyond_ro5.csv', index=False)

print("\n" + "="*70)
if overall_accuracy >= 90:
    print("✅ VALIDATION PASSED: Lipinski filter working correctly")
else:
    print("⚠️  VALIDATION WARNING: Check implementation")
print("="*70 + "\n")