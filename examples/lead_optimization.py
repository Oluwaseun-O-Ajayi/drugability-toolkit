import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from lead_optimizer import LeadOptimizer

print("\n" + "="*70)
print("EXAMPLE: Lead Optimization Workflow")
print("="*70 + "\n")

optimizer = LeadOptimizer()

# Compounds with various issues
compounds = {
    'Lead 1 (Good profile)': {
        'SMILES': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'description': 'Well-optimized NSAID'
    },
    'Lead 2 (High MW)': {
        'SMILES': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O',
        'description': 'Needs MW reduction'
    },
    'Lead 3 (High LogP)': {
        'SMILES': 'CCCCCCCCCCc1ccc(cc1)CCCCCCCCCC',
        'description': 'Too lipophilic'
    },
    'Lead 4 (High TPSA)': {
        'SMILES': 'CC(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(=O)O)C(=O)O',
        'description': 'Poor permeability'
    }
}

print("Lead Optimization Reports:\n")

for name, info in compounds.items():
    print(f"\n{'='*70}")
    print(f"{name.upper()}")
    print(f"Description: {info['description']}")
    print(f"{'='*70}")
    
    optimizer.print_report(info['SMILES'], name)

print("\n" + "="*70)
print("OPTIMIZATION SUMMARY:")
print("="*70)
print("\nKey Lessons:")
print("• Lead 1: Excellent starting point, minimal optimization needed")
print("• Lead 2: Reduce MW by removing non-essential groups")
print("• Lead 3: Add polar groups to reduce LogP")
print("• Lead 4: Protect or methylate polar groups to reduce TPSA")
print("\nGeneral Strategy:")
print("1. Identify critical issues (RO5 violations)")
print("2. Prioritize changes impacting bioavailability")
print("3. Maintain or improve potency during optimization")
print("4. Iterate and re-test modified compounds")
print("="*70 + "\n")