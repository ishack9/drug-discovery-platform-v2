"""
External scoring script for REINVENT 4 RL mode.
Accepts SMILES, returns Vina docking scores.
Usage: python rl_scorer.py smiles1 smiles2 ...
"""
import sys
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
os.chdir(BASE_DIR)

import config
import vina_runner


def score_smiles(smiles_list):
    """
    Computes Vina docking scores for a list of SMILES.
    Score normalized to 0.0 (bad) - 1.0 (good).
    """
    receptor_pdbqt = config.RECEPTOR_PDBQT_PATH
    if not os.path.exists(receptor_pdbqt):
        receptor_pdbqt = vina_runner.prepare_receptor()

    scores = []
    for smi in smiles_list:
        result = vina_runner.dock_single(smi.strip(), receptor_pdbqt)
        if result["success"]:
            # binding_score is 0-10, normalize to 0-1
            scores.append(result["binding_score"] / 10.0)
        else:
            scores.append(0.0)

    return scores


if __name__ == "__main__":
    smiles_list = sys.argv[1:]
    if not smiles_list:
        sys.exit(0)

    scores = score_smiles(smiles_list)

    # Output in REINVENT externalprocess format — one score per line
    for score in scores:
        print(f"{score:.4f}")
