"""Shared helpers for the case-study scripts."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DATA = ROOT / "data"
RESULTS = Path(__file__).resolve().parent / "case_study_results"
RESULTS.mkdir(exist_ok=True)


def sa_value(result):
    """Return the numeric SA score from SyntheticAccessibilityScorer output."""
    return round(result["SA_score"], 2) if isinstance(result, dict) else round(result, 2)
