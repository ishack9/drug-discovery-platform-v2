"""
Ilac Kesif Platformu v2 - Validasyon Scripti

Bilinen EGFR inhibitorlerini pipeline'dan gecerir ve
gercek IC50 degerleriyle karsilastirir.

Referans IC50 degerleri (nM, EGFR kinaz):
  Erlotinib    : 2.0  nM  (Tarceva, FDA onaylı 2004)
  Gefitinib    : 33.0 nM  (Iressa, FDA onaylı 2003)
  Osimertinib  : 1.0  nM  (Tagrisso, FDA onaylı 2015, 3. nesil)
  Lapatinib    : 10.8 nM  (Tykerb, FDA onaylı 2007)
  Afatinib     : 0.5  nM  (Gilotrif, FDA onaylı 2013)
  Imatinib     : 500+ nM  (Gleevec, EGFR icin zayif - negatif kontrol)
"""
import os
import sys
import json
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import vina_runner
import admet_client
import scorer

# ── Referans Ilaçlar ─────────────────────────────────────────────────────────
REFERENCE_DRUGS = [
    {
        "id":       "REF001",
        "name":     "Erlotinib (Tarceva)",
        "smiles":   "COCCOc1cc2ncnc(Nc3cccc(Cl)c3F)c2cc1OCCOC",
        "ic50_nm":  2.0,
        "approved": True,
        "mechanism": "1. nesil EGFR inhibitoru, ATP rekabetci"
    },
    {
        "id":       "REF002",
        "name":     "Gefitinib (Iressa)",
        "smiles":   "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "ic50_nm":  33.0,
        "approved": True,
        "mechanism": "1. nesil EGFR inhibitoru"
    },
    {
        "id":       "REF003",
        "name":     "Osimertinib (Tagrisso)",
        "smiles":   "COc1cc2c(Nc3cccc(NC(=O)/C=C/CN(C)C)c3)ncnc2cc1NC1CCN(C)CC1",
        "ic50_nm":  1.0,
        "approved": True,
        "mechanism": "3. nesil EGFR inhibitoru, T790M mutasyonunu hedefler"
    },
    {
        "id":       "REF004",
        "name":     "Lapatinib (Tykerb)",
        "smiles":   "CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1",
        "ic50_nm":  10.8,
        "approved": True,
        "mechanism": "EGFR/HER2 dual inhibitor"
    },
    {
        "id":       "REF005",
        "name":     "Afatinib (Gilotrif)",
        "smiles":   "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCC1CCN(C)CC1",
        "ic50_nm":  0.5,
        "approved": True,
        "mechanism": "2. nesil kovalent EGFR inhibitoru"
    },
    {
        "id":       "NEG001",
        "name":     "Imatinib (Gleevec) - Negatif Kontrol",
        "smiles":   "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "ic50_nm":  10000.0,
        "approved": True,
        "mechanism": "BCR-ABL inhibitoru, EGFR icin zayif baglanan negatif kontrol"
    },
]


def run_validation():
    """Validasyon pipeline'ini calistirir."""
    print("\n" + "="*70)
    print(f"{'VALIDASYON: BILINEN EGFR INHIBITORLERI':^70}")
    print(f"{'Erlotinib | Gefitinib | Osimertinib | Lapatinib | Afatinib':^70}")
    print("="*70)

    # ── 1. Receptor Hazirla ───────────────────────────────────────────────────
    print("\n[1/4] Receptor hazirlaniyor...")
    try:
        receptor_pdbqt = vina_runner.prepare_receptor()
    except Exception as e:
        print(f"  HATA: {e}")
        return

    # ── 2. Docking ────────────────────────────────────────────────────────────
    print(f"\n[2/4] {len(REFERENCE_DRUGS)} referans ilac icin docking yapiliyor...")
    docked = []
    for drug in REFERENCE_DRUGS:
        print(f"  Docking: {drug['name']}...")
        result = vina_runner.dock_single(drug["smiles"], receptor_pdbqt)
        merged = {**drug, **result}
        docked.append(merged)
        if result["success"]:
            print(f"    Vina: {result['vina_score_kcal']:.2f} kcal/mol "
                  f"-> {result['binding_score']:.1f}/10  "
                  f"| Gercek IC50: {drug['ic50_nm']} nM")
        else:
            print(f"    BASARISIZ")

    # ── 3. ADMET ──────────────────────────────────────────────────────────────
    print(f"\n[3/4] ADMET degerlendirmesi...")
    admet_results = admet_client.evaluate_all(docked)

    # ── 4. Skorlama ───────────────────────────────────────────────────────────
    print(f"\n[4/4] Skorlama ve analiz...")
    all_scored = scorer.compute_scores(admet_results)

    # ── Sonuclar ──────────────────────────────────────────────────────────────
    print("\n" + "="*70)
    print(f"{'VALIDASYON SONUCLARI':^70}")
    print("="*70)
    print(f"\n{'Ilac':<30} {'IC50(nM)':>10} {'Vina':>8} {'Baglanma':>10} {'ADMET':>7} {'Kompozit':>10}")
    print("-"*70)

    # IC50'ye gore sirala (kucuk = daha iyi inhibitor)
    sorted_drugs = sorted(all_scored, key=lambda x: x.get("ic50_nm", 9999))

    for d in sorted_drugs:
        neg = " (NEG)" if d["id"].startswith("NEG") else ""
        print(f"{d['name'][:29]:<30} "
              f"{d.get('ic50_nm', 0):>10.1f} "
              f"{d.get('vina_score_kcal', 0):>8.2f} "
              f"{d.get('binding_score', 0):>10.1f} "
              f"{d.get('admet_score', 0):>7.1f} "
              f"{d.get('composite_score', 0):>10.2f}{neg}")

    print("-"*70)

    # ── Korelasyon Analizi ────────────────────────────────────────────────────
    positive = [d for d in sorted_drugs if not d["id"].startswith("NEG") and d.get("success")]
    negative = [d for d in sorted_drugs if d["id"].startswith("NEG") and d.get("success")]

    if positive and negative:
        avg_pos_vina = sum(d["vina_score_kcal"] for d in positive) / len(positive)
        neg_vina     = negative[0]["vina_score_kcal"] if negative else 0

        print(f"\n  Pozitif kontroller ort. Vina: {avg_pos_vina:.2f} kcal/mol")
        print(f"  Negatif kontrol Vina:         {neg_vina:.2f} kcal/mol")

        if avg_pos_vina < neg_vina:
            print(f"\n  BASARILI: Bilinen inhibitorler negatif kontrolden")
            print(f"  daha iyi baglanma skoru aldi.")
            print(f"  Platform EGFR inhibitorlerini dogru ayirt edebiliyor!")
        else:
            print(f"\n  DIKKAT: Negatif kontrol beklenenden iyi skor aldi.")
            print(f"  Binding site koordinatlari gozden gecirilmeli.")

    # En iyi FDA onaylı ilacin platformdaki sirasi
    active = [d for d in sorted_drugs if not d["id"].startswith("NEG")]
    if active:
        best = min(active, key=lambda x: x.get("ic50_nm", 9999))
        platform_rank = sorted(active,
                               key=lambda x: x.get("composite_score", 0),
                               reverse=True)
        rank = next((i+1 for i, d in enumerate(platform_rank)
                     if d["id"] == best["id"]), None)
        print(f"\n  En potent ilac ({best['name'][:25]}, IC50={best['ic50_nm']}nM)")
        print(f"  Platform sirasi: {rank}/{len(active)}")

    # ── JSON Kaydet ───────────────────────────────────────────────────────────
    os.makedirs("results", exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = f"results/validation_{ts}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({
            "timestamp":  ts,
            "target":     "EGFR",
            "results":    all_scored,
        }, f, ensure_ascii=False, indent=2)
    print(f"\n  Validasyon sonuclari kaydedildi: {out_path}")
    print("="*70 + "\n")


if __name__ == "__main__":
    run_validation()
