"""
Drug Discovery Platform v2 — Configuration
Tools: REINVENT 4 | AutoDock Vina | admet-ai
"""

# ── REINVENT 4 ────────────────────────────────────────────────────────────────
REINVENT_CONDA_ENV   = "reinvent4"              # conda environment name
REINVENT_PRIOR_PATH  = "models/reinvent.prior"  # prior model path
REINVENT_NUM_SMILES  = 60                       # number of SMILES to generate

# ── AutoDock Vina ─────────────────────────────────────────────────────────────
RECEPTOR_PDB_PATH    = "receptor/protein.pdb"   # input PDB file — SET THIS
RECEPTOR_PDBQT_PATH  = "receptor/receptor.pdbqt"  # auto-generated

BINDING_SITE = {
    "center_x":  26.5,    # ← get from PyMOL or RCSB
    "center_y":  11.5,
    "center_z":  -1.0,
    "size_x":   25.0,
    "size_y":   25.0,
    "size_z":   25.0,
}

VINA_EXHAUSTIVENESS  = 8   # higher = more accurate but slower (4-32)
VINA_NUM_POSES       = 5

# Vina score normalization bounds (kcal/mol)
VINA_BEST_SCORE      = -12.0  # maps to 10/10
VINA_WORST_SCORE     =  -2.0  # maps to 0/10

# ── ADMET API ─────────────────────────────────────────────────────────────────
ADMETLAB_API_URL     = "https://admetlab3.scbdd.com/api/admet"
ADMETLAB_TIMEOUT     = 60      # seconds
ADMETLAB_BATCH_SIZE  = 20      # SMILES per request

# ── Scoring Weights (must sum to 1.0) ────────────────────────────────────────
WEIGHTS = {
    "binding_score":           0.40,   # AutoDock Vina
    "admet_score":             0.35,   # admet-ai
    "manufacturability_score": 0.15,   # RDKit SA Score
    "novelty_score":           0.10,   # Tanimoto diversity
}

# ── Selection Thresholds ──────────────────────────────────────────────────────
MIN_COMPOSITE_SCORE  = 5.5
MIN_BINDING_SCORE    = 5.0     # ~-7 kcal/mol
MIN_ADMET_SCORE      = 5.0

TOP_N_FOR_LAB        = 10

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR           = "results"

# ── REINVENT 4 RL Mode ────────────────────────────────────────────────────────
REINVENT_AGENT_PATH  = "models/egfr_agent.prior"  # RL agent (auto-created)
RL_MAX_STEPS         = 50    # fast test: 50, full training: 200-500
RL_TIMEOUT           = 600   # 10 minutes (seconds)
USE_RL_MODE          = True  # True = RL mode, False = sampling mode
