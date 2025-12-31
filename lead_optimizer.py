"""
Lead Optimization Tool
======================

Suggest molecular modifications to improve drug-likeness and ADMET properties.

Author: Oluwaseun O. Ajayi
Institution: University of Georgia
Email: oluwaseun.ajayi@uga.edu

References:
-----------
1. Gleeson (2008) J Med Chem 51:817-834
   "Generation of a set of simple, interpretable ADMET rules of thumb"

2. Johnson et al. (2009) Bioorg Med Chem Lett 19:5560-5564
   "Using the Golden Triangle to optimize clearance and oral absorption"
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
from molecular_descriptors import MolecularDescriptorCalculator
from druglikeness_filters import DrugLikenessFilter
from admet_predictor import ADMETPredictor
from typing import Dict, List, Tuple
import numpy as np


class LeadOptimizer:
    """
    Analyze molecules and suggest optimizations for better drug-likeness.
    
    Provides actionable recommendations based on medicinal chemistry principles.
    """
    
    def __init__(self):
        """Initialize optimizer with component tools."""
        self.descriptor_calc = MolecularDescriptorCalculator()
        self.filter_tool = DrugLikenessFilter()
        self.admet_predictor = ADMETPredictor()
        
    def analyze_compound(self, smiles: str) -> Dict[str, any]:
        """
        Comprehensive analysis of a lead compound.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary with analysis and recommendations
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        # Get all properties
        descriptors = self.descriptor_calc.calculate_all_descriptors(smiles)
        lipinski_pass, lipinski_details = self.filter_tool.lipinski_rule_of_five(smiles)
        veber_pass, veber_details = self.filter_tool.veber_rules(smiles)
        admet = self.admet_predictor.predict_all(smiles)
        
        # Identify issues
        issues = self._identify_issues(descriptors, lipinski_details, veber_details, admet)
        
        # Generate recommendations
        recommendations = self._generate_recommendations(issues, descriptors)
        
        # Calculate optimization priority
        priority = self._calculate_priority(issues)
        
        return {
            'descriptors': descriptors,
            'lipinski': lipinski_details,
            'veber': veber_details,
            'admet': admet,
            'issues': issues,
            'recommendations': recommendations,
            'priority': priority,
            'overall_assessment': self._overall_assessment(issues, priority)
        }
    
    def _identify_issues(self, descriptors: Dict, lipinski: Dict, 
                        veber: Dict, admet: Dict) -> List[Dict]:
        """
        Identify drug-likeness issues.
        
        Returns list of issues with severity ratings.
        """
        issues = []
        
        # Molecular weight issues
        if descriptors['MW'] > 500:
            severity = 'high' if descriptors['MW'] > 600 else 'medium'
            issues.append({
                'property': 'Molecular Weight',
                'value': descriptors['MW'],
                'threshold': 500,
                'severity': severity,
                'impact': 'Poor oral absorption, increased clearance'
            })
        
        # LogP issues
        if descriptors['LogP'] > 5:
            severity = 'high' if descriptors['LogP'] > 6 else 'medium'
            issues.append({
                'property': 'LogP',
                'value': descriptors['LogP'],
                'threshold': 5,
                'severity': severity,
                'impact': 'Poor solubility, potential toxicity'
            })
        elif descriptors['LogP'] < 0:
            issues.append({
                'property': 'LogP',
                'value': descriptors['LogP'],
                'threshold': 0,
                'severity': 'low',
                'impact': 'Poor membrane permeability'
            })
        
        # TPSA issues
        if descriptors['TPSA'] > 140:
            issues.append({
                'property': 'TPSA',
                'value': descriptors['TPSA'],
                'threshold': 140,
                'severity': 'high',
                'impact': 'Poor oral bioavailability'
            })
        
        # Rotatable bonds
        if descriptors['RotBonds'] > 10:
            issues.append({
                'property': 'Rotatable Bonds',
                'value': descriptors['RotBonds'],
                'threshold': 10,
                'severity': 'medium',
                'impact': 'Poor oral bioavailability, binding entropy'
            })
        
        # Hydrogen bond donors/acceptors
        if descriptors['HBD'] > 5:
            issues.append({
                'property': 'H-Bond Donors',
                'value': descriptors['HBD'],
                'threshold': 5,
                'severity': 'medium',
                'impact': 'Poor permeability'
            })
        
        if descriptors['HBA'] > 10:
            issues.append({
                'property': 'H-Bond Acceptors',
                'value': descriptors['HBA'],
                'threshold': 10,
                'severity': 'medium',
                'impact': 'Poor permeability'
            })
        
        # Aromatic rings
        if descriptors['AromaticRings'] > 4:
            issues.append({
                'property': 'Aromatic Rings',
                'value': descriptors['AromaticRings'],
                'threshold': 4,
                'severity': 'low',
                'impact': 'Poor solubility, potential toxicity'
            })
        
        # Fraction Csp3 (too low is problematic)
        if descriptors['FractionCSP3'] < 0.25:
            issues.append({
                'property': 'Fraction Csp3',
                'value': descriptors['FractionCSP3'],
                'threshold': 0.25,
                'severity': 'low',
                'impact': 'Flat molecules, poor developability'
            })
        
        return issues
    
    def _generate_recommendations(self, issues: List[Dict], 
                                 descriptors: Dict) -> List[str]:
        """
        Generate specific optimization recommendations.
        
        Based on medicinal chemistry best practices.
        """
        recommendations = []
        
        for issue in issues:
            prop = issue['property']
            
            if prop == 'Molecular Weight':
                recommendations.append(
                    "• Reduce MW: Remove non-essential functional groups, "
                    "simplify aromatic systems, or fragment into smaller leads"
                )
            
            elif prop == 'LogP':
                if issue['value'] > 5:
                    recommendations.append(
                        "• Reduce LogP: Add polar groups (OH, NH2), "
                        "replace lipophilic groups with polar isosteres, "
                        "or add heteroatoms to aromatic rings"
                    )
                else:
                    recommendations.append(
                        "• Increase LogP: Add alkyl groups or aromatic rings, "
                        "replace polar groups with non-polar equivalents"
                    )
            
            elif prop == 'TPSA':
                recommendations.append(
                    "• Reduce TPSA: Protect or remove OH/NH groups, "
                    "convert to N-methylated or O-methylated derivatives, "
                    "or cyclize to reduce exposed polar atoms"
                )
            
            elif prop == 'Rotatable Bonds':
                recommendations.append(
                    "• Reduce flexibility: Introduce conformational constraints "
                    "(cyclize, add rings), replace flexible linkers with rigid scaffolds"
                )
            
            elif prop == 'H-Bond Donors':
                recommendations.append(
                    "• Reduce HBD: N-methylate amines, convert OH to OMe, "
                    "or replace NH with O or S"
                )
            
            elif prop == 'H-Bond Acceptors':
                recommendations.append(
                    "• Reduce HBA: Remove carbonyl groups where possible, "
                    "replace with CH2 or bioisosteres"
                )
            
            elif prop == 'Aromatic Rings':
                recommendations.append(
                    "• Reduce aromaticity: Saturate rings where possible, "
                    "replace aromatic systems with aliphatic alternatives"
                )
            
            elif prop == 'Fraction Csp3':
                recommendations.append(
                    "• Increase Fsp3: Add saturated rings, replace sp2 carbons "
                    "with sp3, introduce cyclopropyl or cyclobutyl groups"
                )
        
        # Add general recommendations
        if descriptors['QED'] < 0.5:
            recommendations.append(
                "• Overall drug-likeness is low. Consider major scaffold changes "
                "or scaffold hopping to more drug-like frameworks"
            )
        
        # Golden Triangle recommendation (Ref: Johnson et al. 2009)
        if descriptors['MW'] > 300 and descriptors['LogP'] > 3:
            recommendations.append(
                "• Consider Golden Triangle space: Optimize toward MW 200-450 "
                "and LogP -2 to 5 for better balance of potency and ADMET"
            )
        
        return recommendations
    
    def _calculate_priority(self, issues: List[Dict]) -> str:
        """
        Calculate optimization priority based on issue severity.
        
        Returns: 'critical', 'high', 'medium', or 'low'
        """
        if not issues:
            return 'none'
        
        high_severity_count = sum(1 for i in issues if i['severity'] == 'high')
        medium_severity_count = sum(1 for i in issues if i['severity'] == 'medium')
        
        if high_severity_count >= 2:
            return 'critical'
        elif high_severity_count >= 1:
            return 'high'
        elif medium_severity_count >= 2:
            return 'medium'
        else:
            return 'low'
    
    def _overall_assessment(self, issues: List[Dict], priority: str) -> str:
        """Generate overall assessment message."""
        if not issues:
            return "Excellent drug-like properties. No optimization needed."
        
        if priority == 'critical':
            return "Critical issues identified. Major optimization required."
        elif priority == 'high':
            return "Significant issues present. Optimization strongly recommended."
        elif priority == 'medium':
            return "Moderate issues. Consider optimization for better profile."
        else:
            return "Minor issues. Compound is reasonably drug-like."
    
    def suggest_modifications(self, smiles: str) -> List[Dict]:
        """
        Suggest specific molecular modifications.
        
        Args:
            smiles: SMILES string
            
        Returns:
            List of modification suggestions with rationale
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        modifications = []
        
        # Analyze functional groups
        # OH groups
        oh_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]'))
        if len(oh_matches) > 3:
            modifications.append({
                'modification': 'Methylate hydroxyl groups',
                'rationale': 'Reduce HBD count and TPSA for better permeability',
                'example': 'R-OH → R-O-CH3',
                'impact': 'Improved oral bioavailability'
            })
        
        # NH groups
        nh_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'))
        if len(nh_matches) > 2:
            modifications.append({
                'modification': 'N-methylate or acetylate amines',
                'rationale': 'Reduce HBD and improve metabolic stability',
                'example': 'R-NH2 → R-NH-CH3 or R-NH-COCH3',
                'impact': 'Better PK properties'
            })
        
        # Carboxylic acids
        cooh_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[OH]'))
        if cooh_matches:
            modifications.append({
                'modification': 'Convert carboxylic acid to ester or amide',
                'rationale': 'Improve permeability and reduce charge',
                'example': 'R-COOH → R-COOCH3 (ester) or R-CONH2 (amide)',
                'impact': 'Enhanced membrane permeability'
            })
        
        # Long aliphatic chains
        descriptors = self.descriptor_calc.calculate_all_descriptors(smiles)
        if descriptors and descriptors['LogP'] > 5:
            modifications.append({
                'modification': 'Introduce polar groups or heteroatoms',
                'rationale': 'Reduce excessive lipophilicity',
                'example': 'Replace CH2 with O, NH, or add OH/F',
                'impact': 'Improved solubility and reduced toxicity risk'
            })
        
        return modifications
    
    def print_report(self, smiles: str, name: str = "Compound"):
        """
        Print comprehensive optimization report.
        
        Args:
            smiles: SMILES string
            name: Compound name
        """
        analysis = self.analyze_compound(smiles)
        
        if 'error' in analysis:
            print(f"❌ Error: {analysis['error']}")
            return
        
        print(f"\n{'='*70}")
        print(f"LEAD OPTIMIZATION REPORT: {name}")
        print(f"{'='*70}")
        print(f"SMILES: {smiles}\n")
        
        # Priority
        priority = analysis['priority']
        priority_symbol = {
            'critical': '🔴',
            'high': '🟠',
            'medium': '🟡',
            'low': '🟢',
            'none': '✅'
        }.get(priority, '⚪')
        
        print(f"Optimization Priority: {priority_symbol} {priority.upper()}")
        print(f"Assessment: {analysis['overall_assessment']}\n")
        
        # Issues
        if analysis['issues']:
            print("Identified Issues:")
            print("-" * 70)
            for issue in analysis['issues']:
                severity_symbol = {
                    'high': '❗',
                    'medium': '⚠️ ',
                    'low': 'ℹ️ '
                }.get(issue['severity'], '•')
                
                print(f"\n{severity_symbol} {issue['property']}")
                print(f"  Current value: {issue['value']:.2f}")
                print(f"  Threshold: {issue['threshold']}")
                print(f"  Impact: {issue['impact']}")
        else:
            print("✅ No significant issues identified")
        
        # Recommendations
        if analysis['recommendations']:
            print(f"\n{'='*70}")
            print("OPTIMIZATION RECOMMENDATIONS:")
            print("="*70)
            for rec in analysis['recommendations']:
                print(f"\n{rec}")
        
        # Specific modifications
        modifications = self.suggest_modifications(smiles)
        if modifications:
            print(f"\n{'='*70}")
            print("SUGGESTED MODIFICATIONS:")
            print("="*70)
            for i, mod in enumerate(modifications, 1):
                print(f"\n{i}. {mod['modification']}")
                print(f"   Rationale: {mod['rationale']}")
                print(f"   Example: {mod['example']}")
                print(f"   Impact: {mod['impact']}")
        
        print(f"\n{'='*70}\n")


# Example usage
if __name__ == "__main__":
    print("\n" + "="*70)
    print("LEAD OPTIMIZATION TOOL")
    print("Evidence-Based Molecular Optimization")
    print("="*70 + "\n")
    
    optimizer = LeadOptimizer()
    
    # Test with compounds needing optimization
    test_compounds = {
        'Good drug-like (Aspirin)': 'CC(=O)Oc1ccccc1C(=O)O',
        'Needs optimization (High MW)': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O',
        'Poor drug-like (High LogP)': 'CCCCCCCCCCCCCCCCCCCC'
    }
    
    for name, smiles in test_compounds.items():
        optimizer.print_report(smiles, name)
    
    print("✅ Lead optimization analysis completed")
    print("="*70)