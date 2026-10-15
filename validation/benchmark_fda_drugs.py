import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from drugability_toolkit import MolecularDescriptorCalculator
from drugability_toolkit import DrugLikenessFilter
from drugability_toolkit import SyntheticAccessibilityScorer
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

print("\n" + "="*70)
print("BENCHMARK: FDA-Approved Drugs Analysis")
print("Descriptor, rule and SA profiles of ten approved oral drugs")
print("="*70 + "\n")

# Benchmark dataset: FDA-approved oral drugs
fda_benchmark = {
    'Aspirin': {'SMILES': 'CC(=O)Oc1ccccc1C(=O)O', 'class': 'Analgesic'},
    'Metformin': {'SMILES': 'CN(C)C(=N)NC(=N)N', 'class': 'Antidiabetic'},
    'Atorvastatin': {'SMILES': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O', 'class': 'Statin'},
    'Lisinopril': {'SMILES': 'NCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O', 'class': 'ACE Inhibitor'},
    'Amlodipine': {'SMILES': 'CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl', 'class': 'Ca Channel Blocker'},
    'Omeprazole': {'SMILES': 'COc1ccc2[nH]c(nc2c1)[S](=O)Cc1ncc(C)c(OC)c1C', 'class': 'PPI'},
    'Warfarin': {'SMILES': 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O', 'class': 'Anticoagulant'},
    'Losartan': {'SMILES': 'CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)c1ccccc1c1nnn[nH]1', 'class': 'ARB'},
    'Sildenafil': {'SMILES': 'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C', 'class': 'PDE5 Inhibitor'},
    'Sertraline': {'SMILES': 'CN[C@H]1CC[C@H](c2ccc(Cl)c(Cl)c2)c2ccccc12', 'class': 'SSRI'}
}

# Initialize tools
desc_calc = MolecularDescriptorCalculator()
filter_tool = DrugLikenessFilter()
sa_scorer = SyntheticAccessibilityScorer()

print("Analyzing FDA-approved drugs...")
print("="*70 + "\n")

# Analyze all drugs
results = []
for drug_name, drug_data in fda_benchmark.items():
    smiles = drug_data['SMILES']
    
    # Calculate properties
    descriptors = desc_calc.calculate_all_descriptors(smiles)
    ro5_pass, ro5_details = filter_tool.lipinski_rule_of_five(smiles)
    veber_pass, veber_details = filter_tool.veber_rules(smiles)
    sa_result = sa_scorer.calculate_sa_score(smiles)
    
    if descriptors and sa_result:
        results.append({
            'Drug': drug_name,
            'Class': drug_data['class'],
            'MW': descriptors['MW'],
            'LogP': descriptors['LogP'],
            'HBD': descriptors['HBD'],
            'HBA': descriptors['HBA'],
            'TPSA': descriptors['TPSA'],
            'RotBonds': descriptors['RotBonds'],
            'QED': descriptors['QED'],
            'SA_Score': sa_result['SA_score'],
            'RO5_Pass': ro5_pass,
            'RO5_Violations': ro5_details['violations'],
            'Veber_Pass': veber_pass
        })

# Create DataFrame
df = pd.DataFrame(results)

# Display results
print("BENCHMARK RESULTS:")
print("="*70)
print(df[['Drug', 'Class', 'MW', 'LogP', 'QED', 'SA_Score', 'RO5_Pass']].to_string(index=False))

# Statistical analysis
print("\n" + "="*70)
print("STATISTICAL ANALYSIS:")
print("="*70 + "\n")

print("Property Distributions:")
print(df[['MW', 'LogP', 'TPSA', 'QED', 'SA_Score']].describe())

# RO5 compliance
print("\n" + "="*70)
print("RULE OF FIVE COMPLIANCE:")
print("="*70)

ro5_pass_count = df['RO5_Pass'].sum()
ro5_pass_percent = (ro5_pass_count / len(df)) * 100
print(f"\nPassing RO5: {ro5_pass_count}/{len(df)} ({ro5_pass_percent:.1f}%)")
print(f"\nViolation Distribution:")
print(df['RO5_Violations'].value_counts().sort_index())

# Visualizations
print("\n" + "="*70)
print("GENERATING VISUALIZATIONS...")
print("="*70 + "\n")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: MW vs LogP
ax1 = axes[0, 0]
scatter = ax1.scatter(df['MW'], df['LogP'], c=df['QED'], cmap='viridis', 
                     s=100, edgecolor='black', linewidth=1)
ax1.axvline(500, color='red', linestyle='--', alpha=0.5, label='RO5 MW limit')
ax1.axhline(5, color='red', linestyle='--', alpha=0.5, label='RO5 LogP limit')
ax1.set_xlabel('Molecular Weight (Da)', fontweight='bold')
ax1.set_ylabel('LogP', fontweight='bold')
ax1.set_title('MW vs LogP (colored by QED)', fontweight='bold')
ax1.legend()
ax1.grid(alpha=0.3)
plt.colorbar(scatter, ax=ax1, label='QED')

# Plot 2: QED distribution
ax2 = axes[0, 1]
ax2.hist(df['QED'], bins=10, color='steelblue', edgecolor='black', alpha=0.7)
ax2.axvline(df['QED'].mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Mean: {df["QED"].mean():.3f}')
ax2.set_xlabel('QED Score', fontweight='bold')
ax2.set_ylabel('Frequency', fontweight='bold')
ax2.set_title('QED Distribution', fontweight='bold')
ax2.legend()
ax2.grid(alpha=0.3)

# Plot 3: SA Score vs QED
ax3 = axes[1, 0]
ax3.scatter(df['SA_Score'], df['QED'], s=100, color='coral', 
           edgecolor='black', linewidth=1)
ax3.set_xlabel('SA Score', fontweight='bold')
ax3.set_ylabel('QED', fontweight='bold')
ax3.set_title('Synthetic Accessibility vs Drug-likeness', fontweight='bold')
ax3.grid(alpha=0.3)

# Plot 4: TPSA vs rotatable bonds with Veber limits
ax4 = axes[1, 1]
ax4.scatter(df['TPSA'], df['RotBonds'], s=100, color='seagreen', edgecolor='black', linewidth=1)
ax4.axvline(140, color='red', linestyle='--', alpha=0.5, label='Veber TPSA limit')
ax4.axhline(10, color='red', linestyle='--', alpha=0.5, label='Veber rotatable-bond limit')
ax4.set_xlabel('TPSA (A^2)', fontweight='bold')
ax4.set_ylabel('Rotatable bonds', fontweight='bold')
ax4.set_title('Veber criteria', fontweight='bold')
ax4.legend()
ax4.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(Path(__file__).parent / 'fda_benchmark_analysis.png', dpi=300, bbox_inches='tight')
print("✓ Saved: validation/fda_benchmark_analysis.png")

# Save results
df.to_csv(Path(__file__).parent / 'fda_benchmark_results.csv', index=False)
print("✓ Saved: validation/fda_benchmark_results.csv")

print("\n" + "="*70)
print("SUMMARY:")
print("="*70)
print(f"\n• {ro5_pass_percent:.0f}% of FDA-approved drugs pass Lipinski's RO5")
print(f"• Average QED score: {df['QED'].mean():.3f} (range: {df['QED'].min():.3f}-{df['QED'].max():.3f})")
print(f"• SA scores range from {df['SA_Score'].min():.1f} to {df['SA_Score'].max():.1f}")
print("\nBenchmark complete")
print("="*70 + "\n")
