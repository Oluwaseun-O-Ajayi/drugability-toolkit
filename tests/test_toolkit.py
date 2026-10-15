"""Tests for Drugability Toolkit."""

import os
import subprocess
import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, RDConfig, rdMolDescriptors

from drugability_toolkit import (ADMETFlags, DrugLikenessFilter, LeadOptimizer,
                                 MolecularDescriptorCalculator, SyntheticAccessibilityScorer,
                                 structural_alerts)

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ATORVASTATIN = "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
CYCLOSPORINE_LIKE_LONG_ALKANE = "C" * 40  # very lipophilic, fails logP rules
DRUGS = {
    "aspirin": ASPIRIN,
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "atorvastatin": ATORVASTATIN,
    "paclitaxel": "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)"
                  "OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
}


# ------------------------------------------------------------------ descriptors
@pytest.mark.parametrize("smiles", DRUGS.values())
def test_descriptors_match_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    d = MolecularDescriptorCalculator().calculate_all_descriptors(smiles)
    assert d["MW"] == pytest.approx(Descriptors.MolWt(mol), abs=0.01)
    assert d["LogP"] == pytest.approx(Crippen.MolLogP(mol), abs=0.01)
    assert d["HBD"] == Lipinski.NumHDonors(mol)
    assert d["HBA"] == Lipinski.NumHAcceptors(mol)
    assert d["TPSA"] == pytest.approx(rdMolDescriptors.CalcTPSA(mol), abs=0.01)
    assert d["QED"] == pytest.approx(QED.qed(mol), abs=1e-3)


def test_invalid_smiles_is_handled():
    assert MolecularDescriptorCalculator().calculate_all_descriptors("not_a_smiles") is None
    assert ADMETFlags().predict_all("not_a_smiles") == {"error": "Invalid SMILES"}
    assert SyntheticAccessibilityScorer().calculate_sa_score("not_a_smiles") is None


def test_tanimoto_similarity_of_identical_molecules_is_one():
    calc = MolecularDescriptorCalculator()
    assert calc.tanimoto_similarity(ASPIRIN, ASPIRIN) == pytest.approx(1.0)
    assert calc.tanimoto_similarity(ASPIRIN, DRUGS["caffeine"]) < 0.5


# ------------------------------------------------------------------ rule filters
def test_lipinski_counts_violations():
    f = DrugLikenessFilter()
    passes, details = f.lipinski_rule_of_five(ASPIRIN)
    assert passes and details["violations"] == 0
    mol = Chem.MolFromSmiles(ATORVASTATIN)
    expected = (int(Descriptors.MolWt(mol) > 500) + int(Crippen.MolLogP(mol) > 5)
                + int(Lipinski.NumHDonors(mol) > 5) + int(Lipinski.NumHAcceptors(mol) > 10))
    passes, details = f.lipinski_rule_of_five(ATORVASTATIN)
    assert details["violations"] == expected
    assert passes == (expected <= 1)
    passes, details = f.lipinski_rule_of_five(DRUGS["paclitaxel"])
    assert not passes and details["violations"] >= 2


def test_veber_thresholds():
    f = DrugLikenessFilter()
    assert f.veber_rules(ASPIRIN)[0]
    # 12 rotatable bonds exceeds the limit of 10
    assert not f.veber_rules("CCCCCCCCCCCCCCO")[0]


def test_egan_fails_for_very_lipophilic_molecule():
    passes, details = DrugLikenessFilter().egan_rules(CYCLOSPORINE_LIKE_LONG_ALKANE)
    assert not passes
    assert details["LogP"] > 5.88


def test_comprehensive_filter_reports_every_rule():
    result = DrugLikenessFilter().comprehensive_filter(DRUGS["ibuprofen"])
    for key in ("Lipinski_RO5", "Veber", "Ghose", "Egan", "BBB_Penetration", "Lead_Likeness", "Overall"):
        assert key in result


# ------------------------------------------------------------------ synthetic accessibility
@pytest.mark.parametrize("smiles", DRUGS.values())
def test_sa_score_equals_reference_implementation(smiles):
    result = SyntheticAccessibilityScorer().calculate_sa_score(smiles)
    assert result["SA_score"] == pytest.approx(sascorer.calculateScore(Chem.MolFromSmiles(smiles)))
    assert 1.0 <= result["SA_score"] <= 10.0


def test_sa_score_orders_simple_and_complex_molecules():
    scorer = SyntheticAccessibilityScorer()
    assert scorer.calculate_sa_score(ASPIRIN)["SA_score"] < scorer.calculate_sa_score(DRUGS["paclitaxel"])["SA_score"]


# ------------------------------------------------------------------ ADMET flags and alerts
def test_admet_flags_are_qualitative_and_labelled():
    flags = ADMETFlags().predict_all(ASPIRIN)
    for section in ("Absorption", "Distribution", "Metabolism", "Excretion", "Toxicity"):
        assert "Heuristic" in flags[section]["note"]
    assert flags["Absorption"]["HIA_class"] in {"High", "Moderate", "Low"}
    assert flags["Toxicity"]["mutagenicity"] == "Not assessed"


def test_known_pains_motif_is_detected():
    # 5-benzylidene rhodanine, a well-known PAINS chemotype
    rhodanine = "O=C1N(C)C(=S)SC1=Cc1ccccc1"
    assert structural_alerts(Chem.MolFromSmiles(rhodanine))["PAINS"]
    assert structural_alerts(Chem.MolFromSmiles(DRUGS["ibuprofen"]))["PAINS"] == []


def test_lead_optimizer_flags_a_large_lipophilic_compound():
    report = LeadOptimizer().analyze_compound(ATORVASTATIN)
    assert report["issues"]
    assert LeadOptimizer().analyze_compound(ASPIRIN)["priority"] != report["priority"]


# ------------------------------------------------------------------ scripts
SCRIPTS = sorted((ROOT / "examples").glob("*.py")) + sorted((ROOT / "validation").glob("*.py"))


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_example_and_validation_scripts_run(script, tmp_path):
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    proc = subprocess.run([sys.executable, str(script)], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
