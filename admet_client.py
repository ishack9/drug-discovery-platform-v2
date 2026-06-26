"""
Drug Discovery Platform v2 — admet-ai Client
Predicts 50+ ADMET properties locally using admet-ai.
EGFR/cancer-focused scoring with traffic light system.
"""
import math
import config

# ── Property Definitions ──────────────────────────────────────────────────────

# Red=bad(0-3), Yellow=moderate(4-6), Green=good(7-10)
TRAFFIC_THRESHOLDS = {
    "hERG":                  {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "AMES":                  {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "DILI":                  {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "ClinTox":               {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "Carcinogens_Lagunin":   {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "Skin_Reaction":         {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "Bioavailability_Ma":    {"red": 0.4,  "yellow": 0.6,  "invert": False},
    "HIA_Hou":               {"red": 0.4,  "yellow": 0.7,  "invert": False},
    "BBB_Martins":           {"red": 0.3,  "yellow": 0.6,  "invert": False},
    "Lipinski":              {"red": 0.4,  "yellow": 0.7,  "invert": False},
    "QED":                   {"red": 0.3,  "yellow": 0.5,  "invert": False},
    "PAMPA_NCATS":           {"red": 0.3,  "yellow": 0.6,  "invert": False},
    "Pgp_Broccatelli":       {"red": 0.4,  "yellow": 0.7,  "invert": True},
}

# Critical properties and weights for EGFR/cancer
SCORING_WEIGHTS = {
    # Toxicity (negative — should be low)
    "hERG":                  (-2.5, True),   # Cardiotoxicity — critical
    "AMES":                  (-2.0, True),   # Mutagenicity
    "DILI":                  (-3.0, True),   # Liver toxicity — critical
    "ClinTox":               (-1.5, True),   # Clinical toxicity
    "Carcinogens_Lagunin":   (-1.0, True),   # Carcinogenicity
    # Bioavailability (positive — should be high)
    "Bioavailability_Ma":    (2.0,  False),  # Oral bioavailability
    "HIA_Hou":               (1.5,  False),  # Intestinal absorption
    "Lipinski":              (1.5,  False),  # Lipinski rule
    "QED":                   (1.0,  False),  # Drug-likeness
    "PAMPA_NCATS":           (1.0,  False),  # Membrane permeability
}


def _traffic_light(prop: str, value: float) -> str:
    """Returns traffic light color for a property value."""
    if prop not in TRAFFIC_THRESHOLDS:
        return "gray"
    t = TRAFFIC_THRESHOLDS[prop]
    invert = t.get("invert", False)

    if invert:
        # Low = good (for toxicity)
        if value < t["red"]:
            return "green"
        elif value < t["yellow"]:
            return "yellow"
        else:
            return "red"
    else:
        # High = good (for bioavailability)
        if value > t["yellow"]:
            return "green"
        elif value > t["red"]:
            return "yellow"
        else:
            return "red"


def _compute_admet_score(props: dict) -> float:
    """
    Computes EGFR/cancer-focused ADMET score (0-10).
    Penalizes toxicity, rewards bioavailability.
    """
    score = 7.0  # starting score

    for prop, (weight, invert) in SCORING_WEIGHTS.items():
        if prop not in props:
            continue
        val = float(props[prop])
        if math.isnan(val):
            continue

        if invert:
            # Toxicity: high value = penalty
            if val > 0.7:
                score += weight  # weight is negative
            elif val > 0.4:
                score += weight * 0.5
        else:
            # Bioavailability: high value = reward
            if val > 0.7:
                score += weight * 0.5
            elif val < 0.3:
                score += weight * (-0.5)

    return round(max(0.0, min(10.0, score)), 2)


def _format_props(props: dict) -> dict:
    """Returns important properties in clean format."""
    key_props = {
        # Toxicity
        "hERG_cardiotoxicity":      props.get("hERG", None),
        "AMES_mutagenicity":        props.get("AMES", None),
        "DILI_liver_injury":        props.get("DILI", None),
        "ClinTox_clinical":         props.get("ClinTox", None),
        "Carcinogenicity":          props.get("Carcinogens_Lagunin", None),
        # Bioavailability
        "Oral_bioavailability":     props.get("Bioavailability_Ma", None),
        "Intestinal_absorption":    props.get("HIA_Hou", None),
        "Membrane_permeability":    props.get("PAMPA_NCATS", None),
        "BBB_penetration":          props.get("BBB_Martins", None),
        "Pgp_substrate":            props.get("Pgp_Broccatelli", None),
        # Drug-likeness
        "Lipinski_rule":            props.get("Lipinski", None),
        "QED_drug_likeness":        props.get("QED", None),
        "Caco2_permeability":       props.get("Caco2_Wang", None),
        # Physicochemical
        "Aqueous_solubility":       props.get("Solubility_AqSolDB", None),
        "LogP":                     props.get("logP", None),
        "TPSA":                     props.get("tpsa", None),
        "MW":                       props.get("molecular_weight", None),
        # Metabolism
        "CYP3A4_inhibitor":         props.get("CYP3A4_Veith", None),
        "CYP2D6_inhibitor":         props.get("CYP2D6_Veith", None),
        "Half_life":                props.get("Half_Life_Obach", None),
        "LD50":                     props.get("LD50_Zhu", None),
    }
    return {k: round(float(v), 3) if v is not None and not math.isnan(float(v)) else None
            for k, v in key_props.items()}


def _traffic_summary(props: dict) -> dict:
    """Returns traffic light color for each critical property."""
    raw_keys = {
        "hERG_cardiotoxicity":   "hERG",
        "AMES_mutagenicity":     "AMES",
        "DILI_liver_injury":     "DILI",
        "ClinTox_clinical":      "ClinTox",
        "Oral_bioavailability":  "Bioavailability_Ma",
        "Intestinal_absorption": "HIA_Hou",
        "Membrane_permeability": "PAMPA_NCATS",
        "Lipinski_rule":         "Lipinski",
        "QED_drug_likeness":     "QED",
    }
    result = {}
    for display_key, raw_key in raw_keys.items():
        val = props.get(raw_key)
        if val is not None and not math.isnan(float(val)):
            result[display_key] = _traffic_light(raw_key, float(val))
        else:
            result[display_key] = "gray"
    return result


def evaluate_all(candidates: list) -> list:
    """
    Evaluates all candidates with admet-ai.

    Args:
        candidates: list of dicts containing smiles and id

    Returns:
        List with admet_score, admet_props, admet_traffic added per candidate
    """
    try:
        from admet_ai import ADMETModel
    except ImportError:
        print("  ERROR: admet-ai not installed! pip install admet-ai")
        return [{**c, "admet_score": 5.0, "admet_source": "none",
                 "admet_props": {}, "admet_traffic": {}} for c in candidates]

    print(f"  Loading admet-ai model...")
    model = ADMETModel()

    smiles_list = [c["smiles"] for c in candidates]
    total = len(smiles_list)
    batch_size = config.ADMETLAB_BATCH_SIZE
    all_results = []

    print(f"  Predicting {total} candidates ({(total + batch_size - 1) // batch_size} batches)...")

    for i in range(0, total, batch_size):
        batch_smiles = smiles_list[i:i + batch_size]
        batch_no = i // batch_size + 1
        total_b = (total + batch_size - 1) // batch_size
        print(f"  Batch {batch_no}/{total_b} ({len(batch_smiles)} candidates)...")

        try:
            df = model.predict(smiles=batch_smiles)
            for j, row in df.iterrows():
                props = row.to_dict()
                all_results.append(props)
        except Exception as e:
            print(f"  Batch {batch_no} error: {e}")
            for _ in batch_smiles:
                all_results.append({})

    results = []
    for c, props in zip(candidates, all_results):
        admet_score   = _compute_admet_score(props)
        admet_props   = _format_props(props)
        admet_traffic = _traffic_summary(props)

        results.append({
            **c,
            "admet_score":   admet_score,
            "admet_source":  "admet-ai",
            "admet_props":   admet_props,
            "admet_traffic": admet_traffic,
        })

    green_count  = sum(1 for r in results if r["admet_score"] >= 7.0)
    yellow_count = sum(1 for r in results if 5.0 <= r["admet_score"] < 7.0)
    red_count    = sum(1 for r in results if r["admet_score"] < 5.0)
    print(f"  admet-ai complete: Green={green_count} Yellow={yellow_count} Red={red_count}")

    return results
