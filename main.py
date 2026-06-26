"""
Drug Discovery Platform v2 — Main Orchestrator

Usage:
  python main.py EGFR
  python main.py ACE2 --top 15
  python main.py BRAF --prepare-receptor   # prepare receptor only
  python main.py --check                   # check installation
"""
import sys
import os
import json
import argparse
from datetime import datetime

import config
import reinvent_runner

# ── Reference Drugs (Validation) ─────────────────────────────────────────────
REFERENCE_DRUGS = [
    {"id": "REF001", "name": "Erlotinib (Tarceva)",    "smiles": "COCCOc1cc2ncnc(Nc3cccc(Cl)c3F)c2cc1OCCOC",                              "ic50_nm": 2.0,     "neg_control": False},
    {"id": "REF002", "name": "Gefitinib (Iressa)",      "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",                        "ic50_nm": 33.0,    "neg_control": False},
    {"id": "REF003", "name": "Osimertinib (Tagrisso)",  "smiles": "COc1cc2c(Nc3cccc(NC(=O)/C=C/CN(C)C)c3)ncnc2cc1NC1CCN(C)CC1",            "ic50_nm": 1.0,     "neg_control": False},
    {"id": "REF004", "name": "Afatinib (Gilotrif)",     "smiles": "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCC1CCN(C)CC1",         "ic50_nm": 0.5,     "neg_control": False},
    {"id": "NEG001", "name": "Imatinib (NEG control)",  "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",           "ic50_nm": 10000.0, "neg_control": True},
]

import vina_runner
import admet_client
import scorer
import reporter


def parse_args():
    p = argparse.ArgumentParser(
        description="Drug Discovery Platform v2 — REINVENT 4 | Vina | admet-ai"
    )
    p.add_argument("protein",            nargs="?",  help="Target protein name (e.g. EGFR, ACE2)")
    p.add_argument("--top",              type=int,   default=None)
    p.add_argument("--output",           default=config.OUTPUT_DIR)
    p.add_argument("--no-pdf",           action="store_true")
    p.add_argument("--prepare-receptor", action="store_true", help="Prepare receptor PDB → PDBQT and exit")
    p.add_argument("--check",            action="store_true", help="Check tool installations")
    return p.parse_args()


def banner():
    print("\n" + "═" * 66)
    print(f"{'DRUG DISCOVERY PLATFORM v2':^66}")
    print(f"{'REINVENT 4  |  AutoDock Vina  |  admet-ai':^66}")
    print("═" * 66)


def check_installation():
    """Checks if all required tools are installed."""
    print("\nInstallation Check")
    print("─" * 40)
    ok = True

    try:
        from rdkit import Chem
        print("  ✓ RDKit")
    except ImportError:
        print("  ✗ RDKit — pip install rdkit")
        ok = False

    try:
        import meeko
        print("  ✓ meeko")
    except ImportError:
        print("  ✗ meeko — pip install meeko")
        ok = False

    vina_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vina.exe")
    if os.path.exists(vina_exe):
        print(f"  ✓ vina.exe")
    else:
        print(f"  ✗ vina.exe not found — github.com/ccsb-scripps/AutoDock-Vina/releases")
        ok = False

    try:
        import pdbfixer
        print("  ✓ pdbfixer")
    except ImportError:
        print("  ⚠  pdbfixer missing (optional) — pip install pdbfixer")

    try:
        import reportlab
        print("  ✓ reportlab")
    except ImportError:
        print("  ✗ reportlab — pip install reportlab")
        ok = False

    if os.path.exists(config.REINVENT_PRIOR_PATH):
        print(f"  ✓ Prior model: {config.REINVENT_PRIOR_PATH}")
    else:
        print(f"  ✗ Prior model not found: {config.REINVENT_PRIOR_PATH}")
        ok = False

    if os.path.exists(config.RECEPTOR_PDB_PATH):
        print(f"  ✓ Receptor PDB: {config.RECEPTOR_PDB_PATH}")
    else:
        print(f"  ✗ Receptor PDB not found: {config.RECEPTOR_PDB_PATH}")
        ok = False

    bs = config.BINDING_SITE
    if all(v != 0.0 for v in [bs["center_x"], bs["center_y"], bs["center_z"]]):
        print(f"  ✓ Binding site: ({bs['center_x']}, {bs['center_y']}, {bs['center_z']})")
    else:
        print(f"  ⚠  Binding site coordinates are zero — update config.py")

    print("─" * 40)
    print("  READY ✓" if ok else "  Missing installations — see README.md")
    print()
    return ok


def main():
    args = parse_args()
    banner()

    if args.check:
        check_installation()
        return

    if args.prepare_receptor:
        print("\n[Receptor Preparation]")
        pdbqt = vina_runner.prepare_receptor(force=True)
        print(f"  Done: {pdbqt}")
        return

    protein_name = args.protein
    if not protein_name:
        protein_name = input("\n  Target protein name: ").strip()
        if not protein_name:
            print("  ERROR: Protein name is required.")
            sys.exit(1)

    if args.top:
        config.TOP_N_FOR_LAB = args.top

    print(f"\n  Target Protein : {protein_name}")
    print(f"  Generate       : {config.REINVENT_NUM_SMILES} SMILES")
    print(f"  Lab Selection  : Top {config.TOP_N_FOR_LAB}")
    print(f"  Binding Site   : ({config.BINDING_SITE['center_x']}, "
          f"{config.BINDING_SITE['center_y']}, {config.BINDING_SITE['center_z']})")
    print("─" * 66)

    # ══ Step 1 — SMILES Generation ════════════════════════════════════════════
    if config.USE_RL_MODE:
        print(f"\n[1/5] Generating SMILES with REINVENT 4 RL mode...")
        try:
            import reinvent_rl_runner
            candidates = reinvent_rl_runner.generate_candidates_rl()
        except Exception as e:
            print(f"  RL mode failed, falling back to sampling: {e}")
            candidates = reinvent_runner.generate_candidates()
    else:
        print(f"\n[1/5] Generating SMILES with REINVENT 4...")
        try:
            candidates = reinvent_runner.generate_candidates()
        except Exception as e:
            print(f"  ERROR: {e}")
            sys.exit(1)

    # ══ Step 2 — Receptor Preparation ════════════════════════════════════════
    print(f"\n[2/5] Preparing receptor...")
    try:
        receptor_pdbqt = vina_runner.prepare_receptor()
    except Exception as e:
        print(f"  ERROR: {e}")
        sys.exit(1)

    # ══ Step 3 — Molecular Docking ════════════════════════════════════════════
    print(f"\n[3/5] AutoDock Vina docking ({len(candidates)} candidates)...")
    try:
        docked = vina_runner.dock_all(candidates, receptor_pdbqt)
    except Exception as e:
        print(f"  ERROR: {e}")
        sys.exit(1)

    # ══ Step 4 — ADMET Evaluation ════════════════════════════════════════════
    print(f"\n[4/5] ADMET evaluation with admet-ai...")
    try:
        admet_scored = admet_client.evaluate_all(docked)
    except Exception as e:
        print(f"  ERROR: {e}")
        sys.exit(1)

    # ══ Step 5 — Scoring + Reports ════════════════════════════════════════════
    print(f"\n[5/5] Computing scores and generating reports...")

    all_scored     = scorer.compute_scores(admet_scored)
    lab_candidates = scorer.select_lab_candidates(all_scored)
    stats          = scorer.get_statistics(all_scored)

    print(f"  ✓ Ranked: {len(all_scored)}  "
          f"| Avg: {stats.get('avg', 0):.2f}  "
          f"| Lab: {len(lab_candidates)}")

    reporter.print_terminal_report(protein_name, all_scored, lab_candidates)

    csv_path = reporter.save_csv(all_scored, lab_candidates, args.output, protein_name)
    print(f"  ✓ CSV  → {csv_path}")

    # ══ Validation ════════════════════════════════════════════════════════════
    validation_results = []
    if not args.no_pdf:
        print(f"\n[VAL] Validating with known reference drugs...")
        try:
            for drug in REFERENCE_DRUGS:
                dock_r = vina_runner.dock_single(drug["smiles"], receptor_pdbqt)
                validation_results.append({**drug, **dock_r})
                status = f"{dock_r['vina_score_kcal']:.2f} kcal/mol" if dock_r["success"] else "FAILED"
                print(f"  {drug['name'][:30]:<30} IC50={drug['ic50_nm']:>8.1f} nM  Vina={status}")

            pos = [r for r in validation_results if not r["neg_control"] and r.get("success")]
            neg = [r for r in validation_results if r["neg_control"] and r.get("success")]
            if pos and neg:
                avg_pos = sum(r["vina_score_kcal"] for r in pos) / len(pos)
                avg_neg = sum(r["vina_score_kcal"] for r in neg) / len(neg)
                if avg_pos < avg_neg:
                    print(f"  SUCCESS: Known EGFR inhibitors (avg {avg_pos:.1f}) outperformed negative control ({avg_neg:.1f}).")
                else:
                    print(f"  WARNING: Negative control scored unexpectedly high.")
        except Exception as e:
            print(f"  Validation error (continuing): {e}")

        pdf_path = reporter.save_pdf(protein_name, all_scored, lab_candidates, args.output, validation_results)
        print(f"  ✓ PDF  → {pdf_path}")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.makedirs(args.output, exist_ok=True)
    json_path = os.path.join(args.output, f"raw_{protein_name}_{ts}.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({"protein": protein_name, "stats": stats,
                   "all_scored": all_scored, "lab_candidates": lab_candidates},
                  f, ensure_ascii=False, indent=2)
    print(f"  ✓ JSON → {json_path}")

    print(f"\n{'═' * 66}")
    print(f"  Screening complete — {len(lab_candidates)} candidates ready for lab.")
    print(f"  Output: {os.path.abspath(args.output)}/")
    print(f"{'═' * 66}\n")


if __name__ == "__main__":
    main()
