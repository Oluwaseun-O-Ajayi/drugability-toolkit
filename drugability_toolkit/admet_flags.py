"""
Rule-based ADMET flags
======================

Qualitative flags for absorption, distribution, metabolism, excretion and
toxicity (ADMET) computed from a few physicochemical descriptors and from
published structural-alert catalogs.

These flags are heuristics intended for teaching and for a first look at a
compound set. Except for the structural alerts, the thresholds were chosen by
the author to reflect common medicinal-chemistry rules of thumb (for example,
lower topological polar surface area favoring permeability and brain
penetration, and higher lipophilicity favoring protein binding and metabolic
liability). They are not quantitative models, they have not been validated
against experimental ADMET data, and they must not be used to make decisions
about specific compounds without experimental measurement.

Structural alerts use the PAINS (Baell and Holloway, J Med Chem 2010,
53, 2719-2740) and Brenk (Brenk et al., ChemMedChem 2008, 3, 435-444)
catalogs distributed with RDKit.

Author: Oluwaseun O. Ajayi
License: MIT
"""

from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, MolSurf
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

HEURISTIC_NOTE = ("Heuristic flag from physicochemical descriptors; not a validated "
                  "prediction. Confirm experimentally.")


def _catalog(*names):
    params = FilterCatalogParams()
    for name in names:
        params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, name))
    return FilterCatalog(params)


_PAINS = _catalog("PAINS")
_BRENK = _catalog("BRENK")


def structural_alerts(mol: Chem.Mol) -> Dict[str, List[str]]:
    """Return the names of PAINS and Brenk alerts matched by a molecule."""
    return {
        "PAINS": sorted({m.GetDescription() for m in _PAINS.GetMatches(mol)}),
        "Brenk": sorted({m.GetDescription() for m in _BRENK.GetMatches(mol)}),
    }


class ADMETFlags:
    """Qualitative, rule-based ADMET flags for a molecule given as SMILES."""

    def predict_all(self, smiles: str) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        return {
            "SMILES": smiles,
            "Absorption": self.predict_absorption(smiles),
            "Distribution": self.predict_distribution(smiles),
            "Metabolism": self.predict_metabolism(smiles),
            "Excretion": self.predict_excretion(smiles),
            "Toxicity": self.predict_toxicity(smiles),
            "note": HEURISTIC_NOTE,
        }

    def predict_absorption(self, smiles: str) -> Dict[str, Any]:
        """Oral absorption flags.

        HIA_class uses the TPSA <= 140 A^2 cut-off of Veber et al. (J Med Chem 2002,
        45, 2615-2623) together with molecular weight; permeability and P-glycoprotein
        flags are author-defined heuristics."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        tpsa, hbd = MolSurf.TPSA(mol), Lipinski.NumHDonors(mol)
        if tpsa <= 140 and mw <= 500:
            hia = "High"
        elif tpsa <= 140 and mw <= 600:
            hia = "Moderate"
        else:
            hia = "Low"
        if tpsa < 60 and logp > 0:
            perm = "High"
        elif tpsa < 90 and logp > -1:
            perm = "Moderate"
        else:
            perm = "Low"
        pgp = mw > 400 and logp > 3 and hbd < 2
        return {
            "HIA_class": hia,
            "Caco2_permeability": perm,
            "Pgp_substrate": "Possible" if pgp else "Not flagged",
            "TPSA": tpsa,
            "assessment": {"High": "Good", "Moderate": "Moderate", "Low": "Poor"}[hia],
            "note": HEURISTIC_NOTE,
        }

    def predict_distribution(self, smiles: str) -> Dict[str, Any]:
        """Distribution flags from logP, TPSA and molecular weight (author-defined heuristics)."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        logp, tpsa, mw = Crippen.MolLogP(mol), MolSurf.TPSA(mol), Descriptors.MolWt(mol)
        level = "High" if logp > 3 else "Moderate" if logp > 1 else "Low"
        if tpsa < 60 and mw < 400:
            bbb, bbb_note = "High", "Low polar surface area and size favor brain penetration"
        elif tpsa < 90:
            bbb, bbb_note = "Moderate", "Intermediate polar surface area"
        else:
            bbb, bbb_note = "Low", "High polar surface area disfavors brain penetration"
        return {
            "VD_class": level,
            "PPB_class": level,
            "BBB_penetration": bbb,
            "BBB_note": bbb_note,
            "LogP": logp,
            "TPSA": tpsa,
            "note": HEURISTIC_NOTE,
        }

    def predict_metabolism(self, smiles: str) -> Dict[str, Any]:
        """Metabolism flags from size, logP and aromaticity (author-defined heuristics)."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        n_arom = Lipinski.NumAromaticRings(mol)
        return {
            "CYP3A4_substrate": "Possible" if (300 <= mw <= 700 and 2 <= logp <= 5) else "Not flagged",
            "CYP2D6_substrate": "Possible" if (mw < 500 and n_arom >= 1) else "Not flagged",
            "metabolic_liability": "High" if logp > 4 else "Moderate" if logp >= 2 else "Low",
            "note": HEURISTIC_NOTE,
        }

    def predict_excretion(self, smiles: str) -> Dict[str, Any]:
        """Likely dominant elimination route from size and logP (author-defined heuristic)."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        if mw < 300 and logp < 1:
            route = "Likely renal"
        elif mw < 500 and logp < 3:
            route = "Mixed renal and hepatic"
        else:
            route = "Likely hepatic"
        return {"excretion_route": route, "note": HEURISTIC_NOTE}

    def predict_toxicity(self, smiles: str) -> Dict[str, Any]:
        """hERG flag (author-defined heuristic from size, logP and aromatic rings) and
        PAINS and Brenk structural alerts. Mutagenicity is not assessed."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        n_arom = Lipinski.NumAromaticRings(mol)
        herg_score = int(mw > 400) + int(logp > 4) + int(n_arom >= 3)
        alerts = structural_alerts(mol)
        return {
            "hERG_liability": "Elevated" if herg_score >= 2 else "Moderate" if herg_score == 1 else "Low",
            "hERG_risk_score": herg_score,
            "PAINS_alerts": alerts["PAINS"],
            "Brenk_alerts": alerts["Brenk"],
            "mutagenicity": "Not assessed",
            "note": HEURISTIC_NOTE,
        }

    def print_report(self, smiles: str, name: str = "Molecule"):
        r = self.predict_all(smiles)
        if "error" in r:
            print(f"Error: {r['error']}")
            return
        a, d, m, e, t = (r[k] for k in ("Absorption", "Distribution", "Metabolism", "Excretion", "Toxicity"))
        print(f"\n{'=' * 70}\nRULE-BASED ADMET FLAGS: {name}\n{'=' * 70}")
        print(f"SMILES: {smiles}\n")
        print(f"Absorption:    HIA {a['HIA_class']}, permeability {a['Caco2_permeability']}, "
              f"P-gp substrate {a['Pgp_substrate']}")
        print(f"Distribution:  VD {d['VD_class']}, protein binding {d['PPB_class']}, BBB {d['BBB_penetration']}")
        print(f"Metabolism:    CYP3A4 {m['CYP3A4_substrate']}, CYP2D6 {m['CYP2D6_substrate']}, "
              f"liability {m['metabolic_liability']}")
        print(f"Excretion:     {e['excretion_route']}")
        print(f"Toxicity:      hERG {t['hERG_liability']}; PAINS alerts: {', '.join(t['PAINS_alerts']) or 'none'}; "
              f"Brenk alerts: {', '.join(t['Brenk_alerts']) or 'none'}")
        print(f"\nNote: {HEURISTIC_NOTE}\n")


# Backward-compatible name used in version 1.0.
ADMETPredictor = ADMETFlags
