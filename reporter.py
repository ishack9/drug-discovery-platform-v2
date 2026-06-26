"""
İlaç Keşif Platformu v2 — Rapor Üreticisi
Terminal, CSV ve PDF çıktıları oluşturur.
"""
import csv
import os
from datetime import datetime

from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, PageBreak
)
from reportlab.lib.enums import TA_CENTER

import config
import scorer as sc

_NAVY   = colors.HexColor("#1a3a5c")
_BLUE   = colors.HexColor("#2c5282")
_LBLUE  = colors.HexColor("#eef2f7")
_LLBLUE = colors.HexColor("#f5f8fc")
_GREEN  = colors.HexColor("#e8f4e8")
_GRAY   = colors.HexColor("#888888")
_LGRAY  = colors.HexColor("#cccccc")


def _bar(score: float, w: int = 18) -> str:
    n = max(0, min(w, int(score / 10.0 * w)))
    return "█" * n + "░" * (w - n)


# ── Terminal ──────────────────────────────────────────────────────────────────

def print_terminal_report(protein_name: str, all_scored: list, lab_candidates: list):
    stats = sc.get_statistics(all_scored)
    W = 84
    print("\n" + "═" * W)
    print(f"{'İLAÇ KEŞİF PLATFORMU v2 — TARAMA RAPORU':^{W}}")
    print("═" * W)
    print(f"  Hedef Protein  : {protein_name}")
    print(f"  Tarih / Saat   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Araçlar        : REINVENT 4 | AutoDock Vina | ADMETlab 3.0 | RDKit")
    print(f"  Taranan Aday   : {stats.get('total', 0)}")
    print(f"  Ort. Kompozit  : {stats.get('avg', 0):.2f} / 10.00")
    print(f"  Lab'a Seçilen  : {len(lab_candidates)}")
    print("─" * W)
    print(f"  Dağılım → Müke(≥8): {stats.get('excellent', 0)}  "
          f"İyi(6.5-8): {stats.get('good', 0)}  "
          f"Orta(5-6.5): {stats.get('average', 0)}  "
          f"Zayıf(<5): {stats.get('poor', 0)}")
    print("═" * W)

    if not lab_candidates:
        print("  ⚠  Eşik kriterlerini geçen aday bulunamadı.")
        print("═" * W)
        return

    print(f"\n{'LABORATUVARA GÖNDERİLECEK EN İYİ ADAYLAR':^{W}}")
    print("─" * W)
    print(f"{'#':<4} {'ID':<8} {'İsim':<24} "
          f"{'Vina':>10} {'ADMET':>7} {'Üretim':>7} {'Yenilik':>8} "
          f"{'Kompozit':>10}  {'ADMET?':>8}")
    print("─" * W)
    for i, c in enumerate(lab_candidates, 1):
        src = c.get("admet_source", "?")[:8]
        print(f"{i:<4} {c['id']:<8} {c['name'][:23]:<24} "
              f"{c.get('vina_score_kcal', 0):>10.2f} "
              f"{c.get('admet_score', 0):>7.1f} "
              f"{c.get('manufacturability_score', 0):>7.1f} "
              f"{c.get('novelty_score', 0):>8.1f} "
              f"{c['composite_score']:>10.2f}  {src:>8}")
    print("─" * W)

    top = lab_candidates[0]
    print(f"\n{'EN İYİ ADAY — DETAY':^{W}}")
    print("─" * W)
    print(f"  ID        : {top['id']}")
    print(f"  SMILES    : {top['smiles']}")
    print(f"  Vina Skor : {top.get('vina_score_kcal', 'N/A')} kcal/mol")
    print()
    print(f"  Bağlanma  [{_bar(top.get('binding_score', 0))}] {top.get('binding_score', 0):.1f}/10")
    print(f"  ADMET     [{_bar(top.get('admet_score', 0))}] {top.get('admet_score', 0):.1f}/10")
    print(f"  Üretim    [{_bar(top.get('manufacturability_score', 0))}] {top.get('manufacturability_score', 0):.1f}/10")
    print(f"  Yenilik   [{_bar(top.get('novelty_score', 0))}] {top.get('novelty_score', 0):.1f}/10")
    print("\n" + "═" * W)


# ── CSV ───────────────────────────────────────────────────────────────────────

def save_csv(all_scored: list, lab_candidates: list, output_dir: str, protein_name: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe = protein_name.replace(" ", "_")
    path = os.path.join(output_dir, f"drug_v2_{safe}_{ts}.csv")

    lab_ids = {c["id"] for c in lab_candidates}
    fields = [
        "rank", "id", "name", "smiles",
        "vina_score_kcal", "binding_score",
        "admet_score", "admet_source",
        "manufacturability_score", "novelty_score",
        "composite_score", "selected_for_lab",
    ]

    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for rank, c in enumerate(all_scored, 1):
            row = {k: c.get(k, "") for k in fields}
            row["rank"] = rank
            row["selected_for_lab"] = "EVET" if c["id"] in lab_ids else "hayir"
            w.writerow(row)

    return path


# ── PDF ───────────────────────────────────────────────────────────────────────

def _styles():
    base = getSampleStyleSheet()
    return {
        "title":    ParagraphStyle("T", parent=base["Title"],
                                    fontSize=18, spaceAfter=4,
                                    alignment=TA_CENTER, textColor=_NAVY),
        "subtitle": ParagraphStyle("S", parent=base["Normal"],
                                    fontSize=10, spaceAfter=4,
                                    alignment=TA_CENTER, textColor=_BLUE),
        "h1":       ParagraphStyle("H1", parent=base["Heading1"],
                                    fontSize=12, textColor=_NAVY,
                                    spaceBefore=12, spaceAfter=5),
        "h2":       ParagraphStyle("H2", parent=base["Heading2"],
                                    fontSize=10, textColor=_BLUE,
                                    spaceBefore=8, spaceAfter=3),
        "body":     ParagraphStyle("B", parent=base["Normal"],
                                    fontSize=8, leading=12),
        "small":    ParagraphStyle("Sm", parent=base["Normal"],
                                    fontSize=7.5, textColor=colors.HexColor("#333")),
        "footer":   ParagraphStyle("F", parent=base["Normal"],
                                    fontSize=7, textColor=_GRAY,
                                    alignment=TA_CENTER),
    }


def save_pdf(protein_name: str, all_scored: list, lab_candidates: list, output_dir: str, validation_results: list = None) -> str:
    os.makedirs(output_dir, exist_ok=True)
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe = protein_name.replace(" ", "_")
    path = os.path.join(output_dir, f"drug_v2_{safe}_{ts}.pdf")

    doc = SimpleDocTemplate(path, pagesize=A4,
                             leftMargin=2*cm, rightMargin=2*cm,
                             topMargin=2*cm, bottomMargin=2*cm)
    S     = _styles()
    story = []
    now   = datetime.now().strftime("%d.%m.%Y %H:%M")
    stats = sc.get_statistics(all_scored)
    lab_ids = {c["id"] for c in lab_candidates}
    n = stats.get("total", 1) or 1

    # ── Başlık ────────────────────────────────────────────────────────────────
    story.append(Spacer(1, 0.5*cm))
    story.append(Paragraph("ILAC KESIF PLATFORMU v2", S["title"]))
    story.append(Paragraph(
        "REINVENT 4  |  AutoDock Vina  |  ADMETlab 3.0  |  RDKit SA Score",
        S["subtitle"]
    ))
    story.append(HRFlowable(width="100%", thickness=2, color=_NAVY))
    story.append(Spacer(1, 0.4*cm))

    # Meta
    meta = [
        ["Hedef Protein", protein_name],
        ["Rapor Tarihi",  now],
        ["Taranan Aday",  str(stats.get("total", 0))],
        ["Ort. Skor",     f"{stats.get('avg', 0):.2f} / 10.00"],
        ["Lab Secimi",    str(len(lab_candidates))],
        ["Agirliklar",
         "Baglanma(Vina) %40 | ADMET %35 | Uretim(SA) %15 | Yenilik %10"],
    ]
    mt = Table(meta, colWidths=[4.5*cm, 12.5*cm])
    mt.setStyle(TableStyle([
        ("FONTNAME", (0,0), (0,-1), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,-1), 8),
        ("TEXTCOLOR",(0,0),(0,-1), _NAVY),
        ("ROWBACKGROUNDS",(0,0),(-1,-1),[_LBLUE, colors.white]),
        ("GRID",(0,0),(-1,-1),0.5,_LGRAY),
        ("PADDING",(0,0),(-1,-1),5),
    ]))
    story.append(mt)
    story.append(Spacer(1, 0.5*cm))

    # ── Skor Dağılımı ─────────────────────────────────────────────────────────
    story.append(Paragraph("Skor Dagilimi", S["h1"]))
    dist = [
        ["Kategori","Aralik","Aday","Yuzde"],
        ["Mukemmel","8.0-10.0", str(stats.get("excellent",0)),
         f"%{stats.get('excellent',0)/n*100:.1f}"],
        ["Iyi",     "6.5-7.9",  str(stats.get("good",0)),
         f"%{stats.get('good',0)/n*100:.1f}"],
        ["Orta",    "5.0-6.4",  str(stats.get("average",0)),
         f"%{stats.get('average',0)/n*100:.1f}"],
        ["Zayif",   "0-4.9",    str(stats.get("poor",0)),
         f"%{stats.get('poor',0)/n*100:.1f}"],
        ["TOPLAM",  "—",         str(n), ""],
    ]
    dt = Table(dist, colWidths=[3.5*cm, 3*cm, 3.5*cm, 4*cm])
    dt.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),_NAVY),
        ("TEXTCOLOR", (0,0),(-1,0),colors.white),
        ("FONTNAME",  (0,0),(-1,0),"Helvetica-Bold"),
        ("FONTNAME",  (0,-1),(-1,-1),"Helvetica-Bold"),
        ("BACKGROUND",(0,-1),(-1,-1),_LBLUE),
        ("FONTSIZE",  (0,0),(-1,-1),8),
        ("ROWBACKGROUNDS",(0,1),(-1,-2),[colors.white,_LLBLUE]),
        ("GRID",(0,0),(-1,-1),0.5,_LGRAY),
        ("ALIGN",(2,0),(-1,-1),"CENTER"),
        ("PADDING",(0,0),(-1,-1),5),
    ]))
    story.append(dt)

    # ── Lab Adayları ──────────────────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph(f"Laboratuvara Secilen {len(lab_candidates)} Aday", S["h1"]))

    lh = ["#","ID","Isim","Vina(kcal)","ADMET","Uretim","Yenilik","Kompozit"]
    lr = [lh]
    for i, c in enumerate(lab_candidates, 1):
        lr.append([
            str(i), c["id"], c["name"][:28],
            f"{c.get('vina_score_kcal',0):.2f}",
            f"{c.get('admet_score',0):.1f}",
            f"{c.get('manufacturability_score',0):.1f}",
            f"{c.get('novelty_score',0):.1f}",
            f"{c['composite_score']:.2f}",
        ])
    lt = Table(lr, colWidths=[0.8*cm,1.5*cm,5.5*cm,2.2*cm,2*cm,2*cm,2*cm,2*cm],
               repeatRows=1)
    lt.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),_NAVY),
        ("TEXTCOLOR", (0,0),(-1,0),colors.white),
        ("FONTNAME",  (0,0),(-1,0),"Helvetica-Bold"),
        ("FONTSIZE",  (0,0),(-1,-1),7.5),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,_LBLUE]),
        ("GRID",(0,0),(-1,-1),0.5,_LGRAY),
        ("ALIGN",(3,0),(-1,-1),"CENTER"),
        ("PADDING",(0,0),(-1,-1),4),
    ]))
    story.append(lt)
    story.append(Spacer(1, 0.8*cm))

    # ── Bireysel Detaylar ─────────────────────────────────────────────────────
    story.append(Paragraph("Detayli Aday Analizi", S["h1"]))
    for i, c in enumerate(lab_candidates, 1):
        story.append(Paragraph(f"{i}. {c['id']} — {c['name']}", S["h2"]))
        rows = [
            ["SMILES",      c.get("smiles","")],
            ["Vina Skoru",  f"{c.get('vina_score_kcal',0):.2f} kcal/mol  "
                            f"(normalize: {c.get('binding_score',0):.1f}/10)"],
            ["ADMET Skoru", f"{c.get('admet_score',0):.1f}/10  "
                            f"(kaynak: {c.get('admet_source','?')})"],
            ["SA / Uretim", f"{c.get('manufacturability_score',0):.1f}/10"],
            ["Yenilik",     f"{c.get('novelty_score',0):.1f}/10  (Tanimoto cesitliligi)"],
            ["Kompozit",    f"{c.get('composite_score',0):.2f} / 10.00"],
        ]
        # ADMET detay tablosu
        traffic = c.get("admet_traffic", {})
        props   = c.get("admet_props", {})
        traffic_colors = {"green": colors.HexColor("#2ecc71"),
                          "yellow": colors.HexColor("#f39c12"),
                          "red":   colors.HexColor("#e74c3c"),
                          "gray":  colors.HexColor("#bdc3c7")}
        admet_rows = [["Ozellik", "Deger", "Durum"]]
        prop_display = {
            "hERG_kardiyotoksisite":  ("Kalp Ritim Bozuklugu (hERG)", True),
            "AMES_mutajenite":        ("DNA Hasari / Mutajenite (AMES)", True),
            "DILI_karaciger":         ("Karaciger Hasari Riski (DILI)", True),
            "ClinTox_klinik":         ("Klinik Toksisite Riski", True),
            "Oral_biyoyararlilik":    ("Oral Kullanim Uygunlugu", False),
            "Intestinal_absorpsiyon": ("Bagirsaktan Emilim", False),
            "Membran_gecirgenlik":    ("Hucre Zarindan Gecis", False),
            "Lipinski_kurali":        ("Ilac Benzeri Yapi (Lipinski)", False),
            "QED_drug_likeness":      ("Ilac Kalitesi Skoru (QED)", False),
            "Suda_cozunurluk":        ("Suda Cozunurluk", False),
            "LogP":                   ("Yag/Su Dagilimi (LogP)", False),
            "CYP3A4_inhibitor":       ("Karaciger Enzim Blokaji (CYP3A4)", True),
            "Yarim_omur":             ("Vucut Kalis Suresi (saat)", False),
            "LD50":                   ("Oldurucü Doz LD50 (mg/kg)", False),
        }
        for key, (label, is_tox) in prop_display.items():
            val  = props.get(key)
            traf = traffic.get(key, "gray")
            if val is None:
                val_str = "N/A"
            elif isinstance(val, float):
                val_str = f"{val:.3f}"
            else:
                val_str = str(val)
            color = traffic_colors.get(traf, traffic_colors["gray"])
            admet_rows.append([label, val_str,
                               Paragraph(
                                   "<font size='8' color='black'><b>"
                                   + ("iyi (guvenli)" if traf=="green" else "Orta (dikkat)" if traf=="yellow" else "Kotu (toksik)" if traf=="red" else "Veri Yok")
                                   + "</b></font>", S["small"]
                               )])

        admet_tbl = Table(admet_rows, colWidths=[5*cm, 3*cm, 2.5*cm])
        admet_style = [
            ("BACKGROUND", (0,0), (-1,0), _NAVY),
            ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
            ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",   (0,0), (-1,-1), 7.5),
            ("GRID",       (0,0), (-1,-1), 0.3, _LGRAY),
            ("PADDING",    (0,0), (-1,-1), 3),
        ]
        for row_i, (key, (label, is_tox)) in enumerate(prop_display.items(), 1):
            traf = traffic.get(key, "gray")
            bg   = traffic_colors.get(traf, traffic_colors["gray"])
            admet_style.append(("BACKGROUND", (2, row_i), (2, row_i), bg))
        admet_tbl.setStyle(TableStyle(admet_style))

        # ADMET baslik ve tabloyu story'ye ekle
        story.append(Paragraph("ADMET Profili", S["h2"]))
        story.append(admet_tbl)
        story.append(Spacer(1, 0.2*cm))

        para_rows = [
            [Paragraph(f"<b>{r[0]}</b>", S["body"]),
             Paragraph(r[1], S["small"])]
            for r in rows
        ]
        dt2 = Table(para_rows, colWidths=[3*cm, 14*cm])
        dt2.setStyle(TableStyle([
            ("ROWBACKGROUNDS",(0,0),(-1,-1),[_LLBLUE,colors.white]),
            ("GRID",(0,0),(-1,-1),0.3,_LGRAY),
            ("PADDING",(0,0),(-1,-1),4),
            ("VALIGN",(0,0),(-1,-1),"TOP"),
            ("BACKGROUND",(0,-1),(-1,-1),colors.HexColor("#ddeeff")),
        ]))
        story.append(dt2)
        story.append(Spacer(1, 0.25*cm))
        if i % 4 == 0 and i < len(lab_candidates):
            story.append(PageBreak())

    # ── Tam Sıralama (ilk 20) ─────────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph("Tum Adaylar — Ilk 20 Siralama", S["h1"]))

    fh = ["Sira","ID","Isim","Vina","ADMET","Uretim","Yenilik","Komp.","Lab"]
    fr = [fh]
    hl = []
    for rank, c in enumerate(all_scored[:20], 1):
        is_lab = c["id"] in lab_ids
        fr.append([
            str(rank), c["id"], c["name"][:22],
            f"{c.get('vina_score_kcal',0):.1f}",
            f"{c.get('admet_score',0):.1f}",
            f"{c.get('manufacturability_score',0):.1f}",
            f"{c.get('novelty_score',0):.1f}",
            f"{c['composite_score']:.2f}",
            "EVET" if is_lab else "",
        ])
        if is_lab:
            hl.append(("BACKGROUND",(0,rank),(-1,rank),_GREEN))

    ft = Table(fr, colWidths=[1*cm,1.5*cm,4.2*cm,1.8*cm,1.8*cm,1.8*cm,1.8*cm,1.8*cm,1.5*cm],
               repeatRows=1)
    ft.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),_NAVY),
        ("TEXTCOLOR", (0,0),(-1,0),colors.white),
        ("FONTNAME",  (0,0),(-1,0),"Helvetica-Bold"),
        ("FONTSIZE",  (0,0),(-1,-1),7.5),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,_LLBLUE]),
        ("GRID",(0,0),(-1,-1),0.5,_LGRAY),
        ("ALIGN",(3,0),(-1,-1),"CENTER"),
        ("PADDING",(0,0),(-1,-1),3.5),
        *hl,
    ]))
    story.append(ft)

    # ── Validasyon Bolumu ────────────────────────────────────────────────────
    if validation_results:
        story.append(PageBreak())
        story.append(Paragraph("Validasyon: Bilinen EGFR Inhibitorleri", S["h1"]))
        story.append(Paragraph(
            "Platform, FDA onayli EGFR inhibitorlerini negatif kontrolden ayirt edebiliyor mu?",
            S["body"]
        ))
        story.append(Spacer(1, 0.3*cm))

        val_hdr = ["Ilac Adi", "IC50 (nM)", "Vina (kcal/mol)", "Baglanma/10", "Tip"]
        val_rows = [val_hdr]
        pos_scores = []
        neg_scores = []

        for r in sorted(validation_results, key=lambda x: x.get("ic50_nm", 9999)):
            is_neg = r.get("neg_control", False)
            vina   = r.get("vina_score_kcal", 0)
            bind   = r.get("binding_score", 0)
            tip    = "NEG Kontrol" if is_neg else "FDA Onayli"
            val_rows.append([
                r["name"][:35],
                f"{r.get('ic50_nm', 0):.1f}",
                f"{vina:.2f}" if r.get("success") else "BASARISIZ",
                f"{bind:.1f}" if r.get("success") else "-",
                tip,
            ])
            if r.get("success"):
                if is_neg:
                    neg_scores.append(vina)
                else:
                    pos_scores.append(vina)

        val_tbl = Table(val_rows, colWidths=[6*cm, 2.5*cm, 3.5*cm, 2.5*cm, 3*cm], repeatRows=1)
        val_style = [
            ("BACKGROUND", (0,0), (-1,0), _NAVY),
            ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
            ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",   (0,0), (-1,-1), 8),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, _LBLUE]),
            ("GRID", (0,0), (-1,-1), 0.5, _LGRAY),
            ("ALIGN", (1,0), (-1,-1), "CENTER"),
            ("PADDING", (0,0), (-1,-1), 5),
        ]
        # Negatif kontrolü sarı yap
        for ri, r in enumerate(validation_results, 1):
            if r.get("neg_control"):
                val_style.append(("BACKGROUND", (0,ri), (-1,ri), colors.HexColor("#fff3cd")))
        val_tbl.setStyle(TableStyle(val_style))
        story.append(val_tbl)
        story.append(Spacer(1, 0.4*cm))

        # Sonuc
        if pos_scores and neg_scores:
            avg_pos = sum(pos_scores) / len(pos_scores)
            avg_neg = sum(neg_scores) / len(neg_scores)
            if avg_pos < avg_neg:
                result_text = (f"BASARILI: FDA onaylilar ort {avg_pos:.1f} kcal/mol, "
                               f"negatif kontrol {avg_neg:.1f} kcal/mol. "
                               f"Platform EGFR inhibitorlerini dogru ayirt ediyor.")
                result_color = colors.HexColor("#d4edda")
            else:
                result_text = (f"DIKKAT: Negatif kontrol ({avg_neg:.1f}) onaylilarla benzer "
                               f"skor aldi ({avg_pos:.1f}). Binding site gozden gecirilmeli.")
                result_color = colors.HexColor("#fff3cd")

            result_tbl = Table([[Paragraph(result_text, S["body"])]], colWidths=[17.5*cm])
            result_tbl.setStyle(TableStyle([
                ("BACKGROUND", (0,0), (-1,-1), result_color),
                ("GRID", (0,0), (-1,-1), 0.5, _LGRAY),
                ("PADDING", (0,0), (-1,-1), 8),
            ]))
            story.append(result_tbl)
            story.append(Spacer(1, 0.3*cm))

        # Metodolojik not — her zaman goster
        note_text = (
            "Metodolojik Not: Rigid docking skorlari molekul buyuklugunden etkilenir. "
            "Ayni iskelet sinifindaki ilaclar icin goreli siralama anlamlidir. "
            "Mutlak IC50 korelasyonu icin MM-GBSA rescoring veya FEP hesaplamalari onerilir."
        )
        note_tbl = Table([[Paragraph(note_text, S["small"])]], colWidths=[17.5*cm])
        note_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,-1), colors.HexColor("#f8f9fa")),
            ("GRID", (0,0), (-1,-1), 0.5, _LGRAY),
            ("PADDING", (0,0), (-1,-1), 8),
            ("LEFTPADDING", (0,0), (-1,-1), 10),
        ]))
        story.append(note_tbl)

    # Footer
    story.append(Spacer(1, 0.5*cm))
    story.append(HRFlowable(width="100%", thickness=1, color=_LGRAY))
    story.append(Spacer(1, 0.2*cm))
    story.append(Paragraph(
        "REINVENT 4 generatif model | AutoDock Vina moleküler docking | "
        "ADMETlab 3.0 ADMET tahmini | RDKit SA Score sentez skoru. "
        f"Olusturulma: {now}",
        S["footer"]
    ))
    doc.build(story)
    return path