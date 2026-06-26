"""
Drug Discovery Platform v2 — Scoring Module

Metrics:
  - binding_score:           AutoDock Vina (normalized)
  - admet_score:             admet-ai
  - manufacturability_score: RDKit SA Score (synthetic accessibility)
  - novelty_score:           Tanimoto diversity (Morgan fingerprint)
"""
import config

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.Chem.QED import qed as calc_qed
    try:
        from rdkit.Contrib.SA_Score import sascorer
        _SA_AVAILABLE = True
    except ImportError:
        _SA_AVAILABLE = False
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False
    _SA_AVAILABLE = False


# ── Manufacturability: SA Score ───────────────────────────────────────────────

def _sa_score_to_10(smiles: str) -> float:
    """
    Converts RDKit SA Score (1=easy, 10=hard) to platform scale (0-10).
    SA Score 1 → 10/10, SA Score 10 → 0/10
    """
    if not _RDKIT_AVAILABLE:
        return 5.0  # default

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0

    if _SA_AVAILABLE:
        try:
            sa = sascorer.calculateScore(mol)  # 1-10 (low = easy)
            return round((10.0 - sa) / 9.0 * 10.0, 2)
        except Exception:
            pass

    # Fallback to QED if SA Score unavailable (0-1 → 0-10)
    try:
        q = calc_qed(mol)
        return round(q * 10.0, 2)
    except Exception:
        return 5.0


# ── Novelty: Tanimoto Diversity ───────────────────────────────────────────────

def _compute_novelty_scores(candidates: list) -> dict:
    """
    Computes novelty score for each candidate.
    Uses inverse of average Tanimoto similarity to other candidates.
    High similarity → low novelty, high diversity → high novelty.

    Returns:
        {candidate_id: novelty_score (0-10)}
    """
    if not _RDKIT_AVAILABLE or len(candidates) < 2:
        return {c["id"]: 5.0 for c in candidates}

    # Compute Morgan fingerprints
    fps = {}
    for c in candidates:
        mol = Chem.MolFromSmiles(c["smiles"])
        if mol:
            fps[c["id"]] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)

    scores = {}
    for cid, fp in fps.items():
        similarities = []
        for other_id, other_fp in fps.items():
            if other_id == cid:
                continue
            sim = DataStructs.TanimotoSimilarity(fp, other_fp)
            similarities.append(sim)

        if similarities:
            avg_sim = sum(similarities) / len(similarities)
            novelty = (1.0 - avg_sim) * 10.0  # uniqueness → 0-10
        else:
            novelty = 5.0

        scores[cid] = round(min(10.0, max(0.0, novelty)), 2)

    # Assign average score to candidates without fingerprints
    avg = sum(scores.values()) / len(scores) if scores else 5.0
    for c in candidates:
        if c["id"] not in scores:
            scores[c["id"]] = avg

    return scores


# ── Composite Score ───────────────────────────────────────────────────────────

def compute_scores(candidates: list) -> list:
    """
    Computes manufacturability and novelty scores for all candidates,
    combines with binding and ADMET scores into a composite score.

    Args:
        candidates: list of dicts with binding_score and admet_score

    Returns:
        List sorted by composite score (descending)
    """
    novelty_map = _compute_novelty_scores(candidates)

    scored = []
    for c in candidates:
        cid = c["id"]

        mfg_score     = _sa_score_to_10(c.get("smiles", ""))
        novelty_score = novelty_map.get(cid, 5.0)

        composite = (
            c.get("binding_score", 0)  * config.WEIGHTS["binding_score"] +
            c.get("admet_score",   0)  * config.WEIGHTS["admet_score"] +
            mfg_score                  * config.WEIGHTS["manufacturability_score"] +
            novelty_score              * config.WEIGHTS["novelty_score"]
        )

        scored.append({
            **c,
            "manufacturability_score": mfg_score,
            "novelty_score":           novelty_score,
            "composite_score":         round(composite, 2),
        })

    scored.sort(key=lambda x: x["composite_score"], reverse=True)
    return scored


def select_lab_candidates(scored: list) -> list:
    """Applies threshold filter and returns top N candidates."""
    filtered = [
        c for c in scored
        if c.get("composite_score", 0) >= config.MIN_COMPOSITE_SCORE
        and c.get("binding_score", 0)  >= config.MIN_BINDING_SCORE
        and c.get("admet_score", 0)    >= config.MIN_ADMET_SCORE
    ]
    return filtered[:config.TOP_N_FOR_LAB]


def get_statistics(scored: list) -> dict:
    """Returns summary statistics."""
    if not scored:
        return {}
    composites = [c["composite_score"] for c in scored]
    return {
        "total":     len(scored),
        "avg":       round(sum(composites) / len(composites), 2),
        "max":       max(composites),
        "min":       min(composites),
        "excellent": sum(1 for s in composites if s >= 8.0),
        "good":      sum(1 for s in composites if 6.5 <= s < 8.0),
        "average":   sum(1 for s in composites if 5.0 <= s < 6.5),
        "poor":      sum(1 for s in composites if s < 5.0),
    }
