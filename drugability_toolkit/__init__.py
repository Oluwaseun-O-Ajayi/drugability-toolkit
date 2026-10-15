"""Drugability Toolkit: interpretable, rule-based developability profiles for small molecules.

Author: Oluwaseun O. Ajayi
License: MIT
"""

from .admet_flags import ADMETFlags, ADMETPredictor, structural_alerts
from .druglikeness_filters import DrugLikenessFilter
from .lead_optimizer import LeadOptimizer
from .molecular_descriptors import MolecularDescriptorCalculator
from .synthetic_accessibility import SyntheticAccessibilityScorer

__version__ = "1.1.0"

__all__ = [
    "ADMETFlags",
    "ADMETPredictor",
    "DrugLikenessFilter",
    "LeadOptimizer",
    "MolecularDescriptorCalculator",
    "SyntheticAccessibilityScorer",
    "structural_alerts",
    "__version__",
]
