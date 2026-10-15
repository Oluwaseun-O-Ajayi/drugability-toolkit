"""
Synthetic accessibility
=======================

The synthetic accessibility (SA) score is computed with the reference
implementation of Ertl and Schuffenhauer (J Cheminform 2009, 1, 8) that is
distributed with RDKit (Contrib/SA_Score). Scores range from 1 (easy) to 10
(very difficult).

The additional structural-complexity terms returned alongside the score
(ring, stereocentre, size, spiro and bridgehead counts) are descriptive
quantities defined in this module; they are not part of the SA score.

Author: Oluwaseun O. Ajayi
License: MIT
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, AllChem
from rdkit.Chem import Fragments
import numpy as np
import math
from typing import Dict, Tuple, Optional
import os
import sys

from rdkit.Chem import RDConfig

_SA_DIR = os.path.join(RDConfig.RDContribDir, 'SA_Score')
if _SA_DIR not in sys.path:
    sys.path.append(_SA_DIR)
try:
    import sascorer  # Ertl and Schuffenhauer reference implementation shipped with RDKit
except ImportError as exc:  # pragma: no cover
    raise ImportError('RDKit Contrib/SA_Score was not found; install RDKit from PyPI or conda-forge') from exc


class SyntheticAccessibilityScorer:
    """
    Calculate synthetic accessibility scores for drug-like molecules.
    
    The SA score ranges from 1 (easy to synthesize) to 10 (very difficult).
    Uses the RDKit reference implementation of Ertl & Schuffenhauer (2009).
    """
    
    def __init__(self):
        """Initialize SA scorer with complexity parameters."""
        # Complexity penalties based on literature
        self.size_penalty_threshold = 45  # atoms
        self.stereo_penalty_weight = 0.3
        
    def calculate_sa_score(self, smiles: str) -> Optional[Dict[str, float]]:
        """
        Calculate synthetic accessibility score.
        
        SA Score interpretation:
        1-3:   Easy to synthesize
        3-5:   Moderately easy
        5-7:   Moderately difficult
        7-10:  Very difficult
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary with SA score and components
            
        Reference:
            Ertl P, Schuffenhauer A (2009) "Estimation of synthetic accessibility
            score of drug-like molecules based on molecular complexity and fragment
            contributions" J Cheminform 1:8
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        sa_score = float(sascorer.calculateScore(mol))

        return {
            'SA_score': sa_score,
            'interpretation': self._interpret_score(sa_score),
            # Descriptive complexity terms; not part of the SA score.
            'complexity': self._calculate_complexity(mol),
            'ring_complexity': self._calculate_ring_complexity(mol),
            'stereo_complexity': self._calculate_stereo_complexity(mol),
            'size_penalty': self._calculate_size_penalty(mol),
            'spiro_complexity': self._calculate_spiro_complexity(mol),
            'bridgehead_complexity': self._calculate_bridgehead_complexity(mol),
        }

    def _calculate_complexity(self, mol: Chem.Mol) -> float:
        """
        Calculate molecular complexity based on graph theory.
        
        Uses number of atoms, bonds, and their types.
        """
        num_atoms = mol.GetNumHeavyAtoms()
        num_bonds = mol.GetNumBonds()
        
        # Calculate bond complexity
        bond_complexity = 0
        for bond in mol.GetBonds():
            bond_type = bond.GetBondType()
            if bond_type == Chem.BondType.DOUBLE:
                bond_complexity += 0.5
            elif bond_type == Chem.BondType.TRIPLE:
                bond_complexity += 1.0
            elif bond.GetIsAromatic():
                bond_complexity += 0.3
        
        # Combine factors
        complexity = math.log(num_atoms + 1) + bond_complexity / 10
        
        return complexity
    
    def _calculate_ring_complexity(self, mol: Chem.Mol) -> float:
        """
        Calculate complexity from ring systems.
        
        More complex ring systems are harder to synthesize.
        """
        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        
        if num_rings == 0:
            return 0.0
        
        # Penalty for multiple rings
        ring_penalty = math.log(num_rings + 1)
        
        # Additional penalty for fused rings
        fused_rings = 0
        rings = ring_info.AtomRings()
        for i, ring1 in enumerate(rings):
            for ring2 in rings[i+1:]:
                if len(set(ring1) & set(ring2)) > 1:  # Fused if share 2+ atoms
                    fused_rings += 1
        
        fused_penalty = fused_rings * 0.5
        
        # Penalty for macrocycles (rings > 7 atoms)
        macro_penalty = 0
        for ring in rings:
            if len(ring) > 7:
                macro_penalty += (len(ring) - 7) * 0.2
        
        return ring_penalty + fused_penalty + macro_penalty
    
    def _calculate_stereo_complexity(self, mol: Chem.Mol) -> float:
        """
        Calculate complexity from stereochemistry.
        
        Chiral centers increase synthetic difficulty.
        """
        # Count chiral centers
        chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        
        if chiral_centers == 0:
            return 0.0
        
        # Penalty increases with number of stereocenters
        stereo_penalty = math.log(chiral_centers + 1) * self.stereo_penalty_weight
        
        return stereo_penalty
    
    def _calculate_size_penalty(self, mol: Chem.Mol) -> float:
        """
        Penalty for very large molecules.
        
        Larger molecules are generally harder to synthesize.
        """
        num_atoms = mol.GetNumHeavyAtoms()
        
        if num_atoms <= self.size_penalty_threshold:
            return 0.0
        
        # Logarithmic penalty for size
        excess_atoms = num_atoms - self.size_penalty_threshold
        size_penalty = math.log(excess_atoms + 1) * 0.5
        
        return size_penalty
    
    def _calculate_spiro_complexity(self, mol: Chem.Mol) -> float:
        """
        Penalty for spiro atoms (shared by two rings at single atom).
        
        Spiro compounds are synthetically challenging.
        """
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        
        spiro_count = 0
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()
            # Count how many rings this atom is in
            rings_containing = sum(1 for ring in rings if atom_idx in ring)
            if rings_containing >= 2:
                # Check if it's a spiro center (in exactly 2 rings with only this atom shared)
                spiro_count += 1
        
        return spiro_count * 0.5
    
    def _calculate_bridgehead_complexity(self, mol: Chem.Mol) -> float:
        """
        Penalty for bridgehead atoms.
        
        Bridged ring systems are synthetically complex.
        """
        ring_info = mol.GetRingInfo()
        
        # Count bridgehead atoms (atoms in 3+ rings)
        bridgehead_count = 0
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()
            num_rings = sum(1 for ring in ring_info.AtomRings() if atom_idx in ring)
            if num_rings >= 3:
                bridgehead_count += 1
        
        return bridgehead_count * 0.7
    
    def _interpret_score(self, sa_score: float) -> str:
        """
        Provide interpretation of SA score.
        
        Args:
            sa_score: SA score (1-10)
            
        Returns:
            Interpretation string
        """
        if sa_score <= 3:
            return "Easy to synthesize"
        elif sa_score <= 5:
            return "Moderately easy"
        elif sa_score <= 7:
            return "Moderately difficult"
        else:
            return "Very difficult to synthesize"
    
    def retrosynthetic_complexity(self, smiles: str) -> Dict[str, any]:
        """
        Estimate retrosynthetic complexity and suggest disconnections.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary with retrosynthetic analysis
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}
        
        # Calculate strategic bonds (rotatable bonds often good disconnection points)
        num_rotatable = Lipinski.NumRotatableBonds(mol)
        
        # Count functional groups
        num_alcohols = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]')))
        num_carbonyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C](=O)')))
        num_amines = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[N;H1,H2]')))
        
        sa_result = self.calculate_sa_score(smiles)
        sa_score = sa_result['SA_score'] if sa_result else None
        
        return {
            'rotatable_bonds': num_rotatable,
            'functional_groups': {
                'alcohols': num_alcohols,
                'carbonyls': num_carbonyls,
                'amines': num_amines
            },
            'SA_score': sa_score,
            'recommendation': self._retrosynthesis_recommendation(sa_score) if sa_score is not None else 'Unknown'
        }
    
    def _retrosynthesis_recommendation(self, sa_score: float) -> str:
        """Provide retrosynthetic strategy recommendation."""
        if sa_score <= 3:
            return "Simple synthesis - consider direct approaches"
        elif sa_score <= 5:
            return "Moderate complexity - a short linear route is likely"
        elif sa_score <= 7:
            return "Complex - consider convergent synthesis"
        else:
            return "Very complex - may require extensive route optimization"
    
    def compare_molecules(self, smiles_list: list) -> list:
        """
        Compare synthetic accessibility of multiple molecules.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            List of dictionaries sorted by SA score
        """
        results = []
        
        for smiles in smiles_list:
            sa_result = self.calculate_sa_score(smiles)
            if sa_result:
                results.append({
                    'SMILES': smiles,
                    'SA_score': sa_result['SA_score'],
                    'interpretation': sa_result['interpretation']
                })
        
        # Sort by SA score (easiest first)
        results.sort(key=lambda x: x['SA_score'])
        
        return results
    
    def print_report(self, smiles: str, name: str = "Molecule"):
        """
        Print detailed SA score report.
        
        Args:
            smiles: SMILES string
            name: Molecule name
        """
        sa_result = self.calculate_sa_score(smiles)
        
        if sa_result is None:
            print(f"❌ Invalid SMILES: {smiles}")
            return
        
        retro = self.retrosynthetic_complexity(smiles)
        
        print(f"\n{'='*70}")
        print(f"SYNTHETIC ACCESSIBILITY REPORT: {name}")
        print(f"{'='*70}")
        print(f"SMILES: {smiles}\n")
        
        print("Synthetic Accessibility Score:")
        print(f"  SA Score:              {sa_result['SA_score']:.2f} / 10")
        print(f"  Interpretation:        {sa_result['interpretation']}\n")
        
        print("Structural complexity terms (descriptive, not part of the SA score):")
        print(f"  Complexity:            {sa_result['complexity']:.2f}")
        print(f"  Ring Complexity:       {sa_result['ring_complexity']:.2f}")
        print(f"  Stereo Complexity:     {sa_result['stereo_complexity']:.2f}")
        print(f"  Size Penalty:          {sa_result['size_penalty']:.2f}")
        print(f"  Spiro Complexity:      {sa_result['spiro_complexity']:.2f}")
        print(f"  Bridgehead Complexity: {sa_result['bridgehead_complexity']:.2f}\n")
        
        print("Retrosynthetic Analysis:")
        print(f"  Rotatable Bonds:       {retro['rotatable_bonds']}")
        print(f"  Functional Groups:     {sum(retro['functional_groups'].values())}")
        print(f"  Strategy:              {retro['recommendation']}")
        
        print(f"{'='*70}\n")


# Example usage and validation
if __name__ == "__main__":
    print("\n" + "="*70)
    print("SYNTHETIC ACCESSIBILITY SCORER")
    print("Ertl and Schuffenhauer SA score (RDKit reference implementation)")
    print("="*70 + "\n")
    
    scorer = SyntheticAccessibilityScorer()
    
    # Test with molecules of varying complexity
    test_molecules = {
        'Aspirin (simple)': 'CC(=O)Oc1ccccc1C(=O)O',
        'Ibuprofen (moderate)': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'Paclitaxel (complex)': 'CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C',
        'Caffeine (easy)': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    }
    
    print("Comparing Synthetic Accessibility:\n")
    
    for name, smiles in test_molecules.items():
        scorer.print_report(smiles, name)
    
    # Comparative analysis
    print("\n" + "="*70)
    print("COMPARATIVE RANKING (Easiest to Hardest):")
    print("="*70 + "\n")
    
    comparison = scorer.compare_molecules(list(test_molecules.values()))
    names_list = list(test_molecules.keys())
    
    for i, result in enumerate(comparison, 1):
        # Find molecule name
        idx = list(test_molecules.values()).index(result['SMILES'])
        mol_name = names_list[idx]
        print(f"{i}. {mol_name}")
        print(f"   SA Score: {result['SA_score']:.2f} - {result['interpretation']}\n")
    
    print("SA scores computed for the example molecules")
    print("="*70)