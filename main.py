"""
İlaç Keşif Platformu v2 — Ana Orkestratör

Kullanım:
  python main.py EGFR
  python main.py ACE2 --top 15
  python main.py BRAF --prepare-receptor   # sadece receptor hazırla
  python main.py --check                   # kurulum kontrolü
"""
import sys
import os
import json
import argparse
from datetime import datetime

import config
import reinvent_runner

# ── Referans İlaçlar (Validasyon) ─────────────────────────────────────────────
REFERENCE_DRUGS = [
    {"id": "REF001", "name": "Erlotinib (Tarceva)",       "smiles": "COCCOc1cc2ncnc(Nc3cccc(Cl)c3F)c2cc1OCCOC",                                                          "ic50_nm": 2.0,     "neg_control": False},
    {"id": "REF002", "name": "Gefitinib (Iressa)",         "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",                                                    "ic50_nm": 33.0,    "neg_control": False},
    {"id": "REF003", "name": "Osimertinib (Tagrisso)",     "smiles": "COc1cc2c(Nc3cccc(NC(=O)/C=C/CN(C)C)c3)ncnc2cc1NC1CCN(C)CC1",                                       "ic50_nm": 1.0,     "neg_control": False},
    {"id": "REF004", "name": "Afatinib (Gilotrif)",        "smiles": "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCC1CCN(C)CC1",                                    "ic50_nm": 0.5,     "neg_control": False},
    {"id": "NEG001", "name": "Imatinib (NEG kontrol)",     "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",                                       "ic50_nm": 10000.0, "neg_control": True},
]
import vina_runner
import admet_client
import scorer
import reporter


def parse_args():
    p = argparse.ArgumentParser(
        description="İlaç Keşif Platformu v2 — REINVENT 4 | Vina | ADMETlab 3.0"
    )
    p.add_argument("protein",            nargs="?",
                   help="Hedef protein adı (örn: EGFR, ACE2)")
    p.add_argument("--top",              type=int, default=None)
    p.add_argument("--output",           default=config.OUTPUT_DIR)
    p.add_argument("--no-pdf",           action="store_true")
    p.add_argument("--prepare-receptor", action="store_true",
                   help="Receptor PDB → PDBQT dönüşümü yap ve çık")
    p.add_argument("--check",            action="store_true",
                   help="Araç kurulumlarını kontrol et")
    return p.parse_args()


def banner():
    print("\n" + "═" * 66)
    print(f"{'İLAÇ KEŞİF PLATFORMU v2':^66}")
    print(f"{'REINVENT 4  |  AutoDock Vina  |  ADMETlab 3.0':^66}")
    print("═" * 66)


def check_installation():
    """Tüm araçların kurulu olup olmadığını kontrol eder."""
    print("\nKurulum Kontrolü")
    print("─" * 40)
    ok = True

    # RDKit
    try:
        from rdkit import Chem
        print("  ✓ RDKit")
    except ImportError:
        print("  ✗ RDKit — pip install rdkit")
        ok = False

    # meeko
    try:
        import meeko
        print("  ✓ meeko")
    except ImportError:
        print("  ✗ meeko — pip install meeko")
        ok = False

    # vina executable
    vina_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vina.exe")
    if os.path.exists(vina_exe):
        print(f"  ✓ vina.exe")
    else:
        print(f"  ✗ vina.exe bulunamadı — github.com/ccsb-scripps/AutoDock-Vina/releases")
        ok = False

    # pdbfixer
    try:
        import pdbfixer
        print("  ✓ pdbfixer")
    except ImportError:
        print("  ⚠  pdbfixer yok (opsiyonel) — pip install pdbfixer")

    # reportlab
    try:
        import reportlab
        print("  ✓ reportlab")
    except ImportError:
        print("  ✗ reportlab — pip install reportlab")
        ok = False

    # requests
    try:
        import requests
        print("  ✓ requests")
    except ImportError:
        print("  ✗ requests — pip install requests")
        ok = False

    # REINVENT 4 conda ortamı
    import subprocess
    result = subprocess.run(
        ["conda", "run", "-n", config.REINVENT_CONDA_ENV,
         "python", "-c", "import reinvent; print('ok')"],
        capture_output=True, text=True
    )
    # libmamba uyarısını görmezden gel, sadece 'ok' var mı bak
    combined = result.stdout + result.stderr
    if "ok" in result.stdout or "ok" in result.stderr:
        print(f"  ✓ REINVENT 4 (conda env: {config.REINVENT_CONDA_ENV})")
    else:
        print(f"  ✗ REINVENT 4 — KURULUM.md Adım 2'ye bakın")
        ok = False

    # Prior model
    if os.path.exists(config.REINVENT_PRIOR_PATH):
        print(f"  ✓ Prior model: {config.REINVENT_PRIOR_PATH}")
    else:
        print(f"  ✗ Prior model bulunamadı: {config.REINVENT_PRIOR_PATH}")
        ok = False

    # Receptor PDB
    if os.path.exists(config.RECEPTOR_PDB_PATH):
        print(f"  ✓ Receptor PDB: {config.RECEPTOR_PDB_PATH}")
    else:
        print(f"  ✗ Receptor PDB bulunamadı: {config.RECEPTOR_PDB_PATH}")
        ok = False

    # Binding site
    bs = config.BINDING_SITE
    if all(v != 0.0 for v in [bs["center_x"], bs["center_y"], bs["center_z"]]):
        print(f"  ✓ Binding site: ({bs['center_x']}, {bs['center_y']}, {bs['center_z']})")
    else:
        print(f"  ⚠  Binding site koordinatları sıfır — config.py'yi güncelleyin")

    print("─" * 40)
    print("  HAZIR ✓" if ok else "  Eksik kurulum var — KURULUM.md'yi okuyun")
    print()
    return ok


def main():
    args = parse_args()
    banner()

    # ── Mod: Kurulum Kontrolü ─────────────────────────────────────────────────
    if args.check:
        check_installation()
        return

    # ── Mod: Sadece Receptor Hazırla ──────────────────────────────────────────
    if args.prepare_receptor:
        print("\n[Receptor Hazırlama]")
        pdbqt = vina_runner.prepare_receptor(force=True)
        print(f"  Tamamlandı: {pdbqt}")
        return

    # ── Ana Tarama Modu ───────────────────────────────────────────────────────
    protein_name = args.protein
    if not protein_name:
        protein_name = input("\n  Hedef protein adı: ").strip()
        if not protein_name:
            print("  HATA: Protein adı gereklidir.")
            sys.exit(1)

    if args.top:
        config.TOP_N_FOR_LAB = args.top

    print(f"\n  Hedef Protein  : {protein_name}")
    print(f"  Üretim Hedefi  : {config.REINVENT_NUM_SMILES} SMILES")
    print(f"  Lab Seçimi     : En iyi {config.TOP_N_FOR_LAB}")
    print(f"  Binding Site   : ({config.BINDING_SITE['center_x']}, "
          f"{config.BINDING_SITE['center_y']}, {config.BINDING_SITE['center_z']})")
    print("─" * 66)

    # ══ 1 — REINVENT 4: SMILES Üretimi ═══════════════════════════════════════
    if config.USE_RL_MODE:
        print(f"\n[1/5] REINVENT 4 RL modu ile EGFR-odakli SMILES uretiliyor...")
        try:
            import reinvent_rl_runner
            candidates = reinvent_rl_runner.generate_candidates_rl()
        except Exception as e:
            print(f"  RL modu hatasi, sampling moduna geciliyor: {e}")
            candidates = reinvent_runner.generate_candidates()
    else:
        print(f"\n[1/5] REINVENT 4 ile SMILES üretiliyor...")
        try:
            candidates = reinvent_runner.generate_candidates()
        except Exception as e:
            print(f"  HATA: {e}")
            sys.exit(1)

    # ══ 2 — Receptor Hazırlama ════════════════════════════════════════════════
    print(f"\n[2/5] Receptor hazırlanıyor...")
    try:
        receptor_pdbqt = vina_runner.prepare_receptor()
    except Exception as e:
        print(f"  HATA: {e}")
        sys.exit(1)

    # ══ 3 — AutoDock Vina: Docking ════════════════════════════════════════════
    print(f"\n[3/5] AutoDock Vina ile docking ({len(candidates)} aday)...")
    try:
        docked = vina_runner.dock_all(candidates, receptor_pdbqt)
    except Exception as e:
        print(f"  HATA: {e}")
        sys.exit(1)

    # ══ 4 — ADMETlab 3.0: ADMET Değerlendirmesi ══════════════════════════════
    print(f"\n[4/5] ADMETlab 3.0 ile ADMET değerlendirmesi...")
    try:
        admet_scored = admet_client.evaluate_all(docked)
    except Exception as e:
        print(f"  HATA: {e}")
        sys.exit(1)

    # ══ 5 — Skor + Raporlar ═══════════════════════════════════════════════════
    print(f"\n[5/5] Skor hesaplanıyor, raporlar oluşturuluyor...")

    all_scored     = scorer.compute_scores(admet_scored)
    lab_candidates = scorer.select_lab_candidates(all_scored)
    stats          = scorer.get_statistics(all_scored)

    print(f"  ✓ Sıralandı: {len(all_scored)}  "
          f"| Ort: {stats.get('avg', 0):.2f}  "
          f"| Lab: {len(lab_candidates)}")

    reporter.print_terminal_report(protein_name, all_scored, lab_candidates)

    csv_path = reporter.save_csv(all_scored, lab_candidates, args.output, protein_name)
    print(f"  ✓ CSV  → {csv_path}")

    # ══ Validasyon ═══════════════════════════════════════════════════════════
    validation_results = []
    if not args.no_pdf:
        print(f"\n[VAL] Referans ilaçlarla validasyon yapılıyor...")
        try:
            for drug in REFERENCE_DRUGS:
                dock_r = vina_runner.dock_single(drug["smiles"], receptor_pdbqt)
                validation_results.append({**drug, **dock_r})
                status = f"{dock_r['vina_score_kcal']:.2f} kcal/mol" if dock_r["success"] else "BASARISIZ"
                print(f"  {drug['name'][:30]:<30} IC50={drug['ic50_nm']:>8.1f} nM  Vina={status}")

            pos = [r for r in validation_results if not r["neg_control"] and r.get("success")]
            neg = [r for r in validation_results if r["neg_control"] and r.get("success")]
            if pos and neg:
                avg_pos = sum(r["vina_score_kcal"] for r in pos) / len(pos)
                avg_neg = sum(r["vina_score_kcal"] for r in neg) / len(neg)
                if avg_pos < avg_neg:
                    print(f"  BASARILI: Bilinen EGFR inhibitorleri (ort {avg_pos:.1f}) negatif kontrolden ({avg_neg:.1f}) daha iyi skoru aldi.")
                else:
                    print(f"  DIKKAT: Negatif kontrol beklenenden iyi skor aldi.")
        except Exception as e:
            print(f"  Validasyon hatasi (devam ediliyor): {e}")

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
    print(f"  Tarama tamamlandı — {len(lab_candidates)} aday laboratuvara hazır.")
    print(f"  Çıktılar: {os.path.abspath(args.output)}/")
    print(f"{'═' * 66}\n")


if __name__ == "__main__":
    main()