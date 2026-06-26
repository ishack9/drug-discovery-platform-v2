"""
REINVENT 4 RL modu icin harici skorlama scripti.
SMILES alir, Vina docking skoru dondurur.
Kullanim: python rl_scorer.py smiles1 smiles2 ...
"""
import sys
import os
import json

# drug_discovery_v2 klasorunu ekle
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
os.chdir(BASE_DIR)

import config
import vina_runner

def score_smiles(smiles_list):
    """
    SMILES listesi icin Vina docking skoru hesaplar.
    Skor: 0.0 (kotu) - 1.0 (iyi) araliginda normalize edilmis.
    """
    receptor_pdbqt = config.RECEPTOR_PDBQT_PATH
    if not os.path.exists(receptor_pdbqt):
        # Receptor hazirla
        receptor_pdbqt = vina_runner.prepare_receptor()

    scores = []
    for smi in smiles_list:
        result = vina_runner.dock_single(smi.strip(), receptor_pdbqt)
        if result["success"]:
            # binding_score zaten 0-10 araliginda, 0-1'e normalize et
            scores.append(result["binding_score"] / 10.0)
        else:
            scores.append(0.0)

    return scores


if __name__ == "__main__":
    smiles_list = sys.argv[1:]
    if not smiles_list:
        sys.exit(0)

    scores = score_smiles(smiles_list)

    # REINVENT externalprocess formatinda cikti
    # Her satir bir skor
    for score in scores:
        print(f"{score:.4f}")
