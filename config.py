"""
İlaç Keşif Platformu v2 — Yapılandırma
Gerçek araçlar: REINVENT 4 | AutoDock Vina | ADMETlab 3.0
"""

# ── REINVENT 4 ───────────────────────────────────────────────────────────────
REINVENT_CONDA_ENV   = "reinvent4"           # conda ortam adı
REINVENT_PRIOR_PATH  = "models/reinvent.prior"  # prior model dosyası — DOLDUR
REINVENT_NUM_SMILES  = 60                    # üretilecek SMILES sayısı

# ── AutoDock Vina ─────────────────────────────────────────────────────────────
RECEPTOR_PDB_PATH    = r"C:\Users\PC\drug_discovery_v2\receptor\4HJO.pdb"   # ham PDB — DOLDUR
RECEPTOR_PDBQT_PATH  = "receptor/receptor.pdbqt"  # otomatik oluşturulur

BINDING_SITE = {
    "center_x":  26.5,    # ← PyMOL / RCSB'den al
    "center_y":  11.5,    # ← PyMOL / RCSB'den al
    "center_z":  -1.0,    # ← PyMOL / RCSB'den al
    "size_x":   25.0,
    "size_y":   25.0,
    "size_z":   25.0,
}

VINA_EXHAUSTIVENESS  = 8   # yüksek = doğru ama yavaş (4-32 arası)
VINA_NUM_POSES       = 5

# Vina skor normalizasyon sınırları (kcal/mol)
VINA_BEST_SCORE      = -12.0  # bu değer → 10/10
VINA_WORST_SCORE     =  -2.0  # bu değer → 0/10

# ── ADMETlab 3.0 API ──────────────────────────────────────────────────────────
ADMETLAB_API_URL     = "https://admetlab3.scbdd.com/api/admet"
ADMETLAB_TIMEOUT     = 60      # saniye
ADMETLAB_BATCH_SIZE  = 20      # istek başına SMILES sayısı

# ── RDKit SA Score (Üretilebilirlik) ─────────────────────────────────────────
# SA Score: 1 (kolay) → 10 (zor), platforma 0-10 olarak normalize edilir

# ── Skorlama Ağırlıkları (toplamı 1.0 olmalı) ───────────────────────────────
WEIGHTS = {
    "binding_score":           0.40,   # AutoDock Vina
    "admet_score":             0.35,   # ADMETlab 3.0
    "manufacturability_score": 0.15,   # RDKit SA Score
    "novelty_score":           0.10,   # Tanimoto çeşitliliği
}

# ── Seçim Eşikleri ────────────────────────────────────────────────────────────
MIN_COMPOSITE_SCORE  = 5.5
MIN_BINDING_SCORE    = 5.0     # ~-7 kcal/mol
MIN_ADMET_SCORE      = 5.0

TOP_N_FOR_LAB        = 10

# ── Çıktı ────────────────────────────────────────────────────────────────────
OUTPUT_DIR           = "results"

# ── REINVENT 4 RL Modu ───────────────────────────────────────────────────────
REINVENT_AGENT_PATH  = "models/egfr_agent.prior"  # RL agent (otomatik olusturulur)
RL_MAX_STEPS         = 50    # Hizli test: 50, tam egitim: 200-500
RL_TIMEOUT           = 600   # 10 dakika (saniye)
USE_RL_MODE          = True  # True = RL modu, False = sampling modu