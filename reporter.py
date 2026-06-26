"""
Drug Discovery Platform v2 — Report Generator
Generates terminal, CSV and PDF outputs.
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
    print(f"{'DRUG DISCOVERY PLATFORM v2 — SCREENING REPORT':^{W}}")
    print("═" * W)
    print(f"  Target Protein : {protein_name}")
    print(f"  Date / Time    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Tools          : REINVENT 4 | AutoDock Vina | admet-ai | RDKit")
    print(f"  Candidates     : {stats.get('total', 0)}")
    print(f"  Avg Composite  : {stats.get('avg', 0):.2f} / 10.00")
    print(f"  Lab Selected   : {len(lab_candidates)}")
    print("─" * W)
    print(f"  Distribution → Excellent(>=8): {stats.get('excellent', 0)}  "
          f"Good(6.5-8): {stats.get('good', 0)}  "
          f"Average(5-6.5): {stats.get('average', 0)}  "
          f"Poor(<5): {stats.get('poor', 0)}")
    print("═" * W)

    if not lab_candidates:
        print("  No candidates passed threshold criteria.")
        print("═" * W)
        return

    print(f"\n{'TOP CANDIDATES FOR LABORATORY':^{W}}")
    print("─" * W)
    print(f"{'#':<4} {'ID':<8} {'Name':<24} "
          f"{'Vina':>10} {'ADMET':>7} {'Mfg':>7} {'Novelty':>8} "
          f"{'Composite':>10}  {'ADMET?':>8}")
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
    print(f"\n{'BEST CANDIDATE — DETAIL':^{W}}")
    print("─" * W)
    print(f"  ID         : {top['id']}")
    print(f"  SMILES     : {top['smiles']}")
    print(f"  Vina Score : {top.get('vina_score_kcal', 'N/A')} kcal/mol")
    print()
    print(f"  Binding  [{_bar(top.get('binding_score', 0))}] {top.get('binding_score', 0):.1f}/10")
    print(f"  ADMET    [{_bar(top.get('admet_score', 0))}] {top.get('admet_score', 0):.1f}/10")
    print(f"  Mfg      [{_bar(top.get('manufacturability_score', 0))}] {top.get('manufacturability_score', 0):.1f}/10")
    print(f"  Novelty  [{_bar(top.get('novelty_score', 0))}] {top.get('novelty_score', 0):.1f}/10")
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
            row["selected_for_lab"] = "YES" if c["id"] in lab_ids else "no"
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


def save_pdf(protein_name: str, all_scored: list, lab_candidates: list,
             output_dir: str, validation_results: list = None) -> str:
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

    # ── Header ────────────────────────────────────────────────────────────────
    story.append(Spacer(1, 0.5*cm))
    story.append(Paragraph("DRUG DISCOVERY PLATFORM v2", S["title"]))
    story.append(Paragraph(
        "REINVENT 4  |  AutoDock Vina  |  admet-ai  |  RDKit SA Score",
        S["subtitle"]
    ))
    story.append(HRFlowable(width="100%", thickness=2, color=_NAVY))
    story.append(Spacer(1, 0.4*cm))

    # Meta table
    meta = [
        ["Target Protein", protein_name],
        ["Report Date",    now],
        ["Candidates",     str(stats.get("total", 0))],
        ["Avg Score",      f"{stats.get('avg', 0):.2f} / 10.00"],
        ["Lab Selection",  str(len(lab_candidates))],
        ["Weights",
         "Binding(Vina) 40% | ADMET 35% | Manufacturability(SA) 15% | Novelty 10%"],
    ]
    mt = Table(meta, colWidths=[4.5*cm, 12.5*cm])
    mt.setStyle(TableStyle([
        ("FONTNAME",  (0,0), (0,-1), "Helvetica-Bold"),
        ("FONTSIZE",  (0,0), (-1,-1), 8),
        ("ROWBACKGROUNDS", (0,0), (-1,-1), [_LLBLUE, colors.white]),
        ("GRID",      (0,0), (-1,-1), 0.3, _LGRAY),
        ("PADDING",   (0,0), (-1,-1), 5),
    ]))
    story.append(mt)
    story.append(Spacer(1, 0.4*cm))

    # ── Score Distribution ────────────────────────────────────────────────────
    story.append(Paragraph("Score Distribution", S["h1"]))
    dist_hdr = ["Category", "Range", "Candidates", "Percent"]
    dist_data = [
        dist_hdr,
        ["Excellent", "8.0-10.0", str(stats.get("excellent", 0)),
         f"%{stats.get('excellent',0)/n*100:.1f}"],
        ["Good",      "6.5-7.9",  str(stats.get("good", 0)),
         f"%{stats.get('good',0)/n*100:.1f}"],
        ["Average",   "5.0-6.4",  str(stats.get("average", 0)),
         f"%{stats.get('average',0)/n*100:.1f}"],
        ["Poor",      "0-4.9",    str(stats.get("poor", 0)),
         f"%{stats.get('poor',0)/n*100:.1f}"],
        ["TOTAL", "—", str(stats.get("total", 0)), ""],
    ]
    dt = Table(dist_data, colWidths=[4*cm, 3*cm, 3.5*cm, 3*cm])
    dt.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), _NAVY),
        ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
        ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTNAME",   (0,-1),(-1,-1),"Helvetica-Bold"),
        ("FONTSIZE",   (0,0), (-1,-1), 8),
        ("ROWBACKGROUNDS", (0,1), (-1,-2), [colors.white, _LLBLUE]),
        ("GRID",       (0,0), (-1,-1), 0.5, _LGRAY),
        ("ALIGN",      (2,0), (-1,-1), "CENTER"),
        ("PADDING",    (0,0), (-1,-1), 5),
    ]))
    story.append(dt)
    story.append(Spacer(1, 0.4*cm))

    # ── Lab Candidates Table ──────────────────────────────────────────────────
    story.append(Paragraph(f"Lab Candidates — Top {len(lab_candidates)}", S["h1"]))
    if not lab_candidates:
        story.append(Paragraph("No candidates passed threshold criteria.", S["body"]))
    else:
        lh = ["#", "ID", "Name", "Vina(kcal)", "ADMET", "Mfg", "Novelty", "Composite"]
        lr = [lh]
        for i, c in enumerate(lab_candidates, 1):
            lr.append([
                str(i), c["id"], c["name"][:25],
                f"{c.get('vina_score_kcal',0):.2f}",
                f"{c.get('admet_score',0):.1f}",
                f"{c.get('manufacturability_score',0):.1f}",
                f"{c.get('novelty_score',0):.1f}",
                f"{c['composite_score']:.2f}",
            ])
        lt = Table(lr, colWidths=[1*cm,1.8*cm,4.5*cm,2.5*cm,2*cm,2*cm,2*cm,2.2*cm],
                   repeatRows=1)
        lt.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), _NAVY),
            ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
            ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",   (0,0), (-1,-1), 8),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, _LLBLUE]),
            ("GRID",       (0,0), (-1,-1), 0.5, _LGRAY),
            ("ALIGN",      (3,0), (-1,-1), "CENTER"),
            ("PADDING",    (0,0), (-1,-1), 4),
        ]))
        story.append(lt)
    story.append(Spacer(1, 0.4*cm))

    # ── Detailed Candidate Analysis ───────────────────────────────────────────
    story.append(Paragraph("Detailed Candidate Analysis", S["h1"]))

    for i, c in enumerate(lab_candidates, 1):
        story.append(Paragraph(f"{i}. {c['id']} — {c['name']}", S["h2"]))

        rows = [
            ["SMILES",         c.get("smiles","")],
            ["Vina Score",     f"{c.get('vina_score_kcal',0):.2f} kcal/mol  "
                               f"(normalized: {c.get('binding_score',0):.1f}/10)"],
            ["ADMET Score",    f"{c.get('admet_score',0):.1f}/10  "
                               f"(source: {c.get('admet_source','?')})"],
            ["SA / Mfg",       f"{c.get('manufacturability_score',0):.1f}/10"],
            ["Novelty",        f"{c.get('novelty_score',0):.1f}/10  (Tanimoto diversity)"],
            ["Composite",      f"{c.get('composite_score',0):.2f} / 10.00"],
        ]

        # ADMET detail table
        traffic = c.get("admet_traffic", {})
        props   = c.get("admet_props", {})
        traffic_colors = {
            "green":  colors.HexColor("#2ecc71"),
            "yellow": colors.HexColor("#f39c12"),
            "red":    colors.HexColor("#e74c3c"),
            "gray":   colors.HexColor("#bdc3c7"),
        }
        admet_rows = [["Property", "Value", "Status"]]
        prop_display = {
            "hERG_cardiotoxicity":   ("Cardiac Toxicity (hERG)", True),
            "AMES_mutagenicity":     ("DNA Damage / Mutagenicity (AMES)", True),
            "DILI_liver_injury":     ("Drug-Induced Liver Injury (DILI)", True),
            "ClinTox_clinical":      ("Clinical Toxicity Risk", True),
            "Oral_bioavailability":  ("Oral Bioavailability", False),
            "Intestinal_absorption": ("Intestinal Absorption", False),
            "Membrane_permeability": ("Membrane Permeability", False),
            "Lipinski_rule":         ("Lipinski Rule of Five", False),
            "QED_drug_likeness":     ("Drug-Likeness Score (QED)", False),
            "Aqueous_solubility":    ("Aqueous Solubility", False),
            "LogP":                  ("Lipophilicity (LogP)", False),
            "CYP3A4_inhibitor":      ("Liver Enzyme Inhibition (CYP3A4)", True),
            "Half_life":             ("Half-Life (hours)", False),
            "LD50":                  ("Lethal Dose LD50 (mg/kg)", False),
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

            status = "Good" if traf == "green" else "Moderate" if traf == "yellow" else "Poor" if traf == "red" else "No data"
            admet_rows.append([label, val_str,
                               Paragraph(
                                   f"<font size='8' color='black'><b>{status}</b></font>",
                                   S["small"]
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

        # Add ADMET table to story
        story.append(Paragraph("ADMET Profile", S["h2"]))
        story.append(admet_tbl)
        story.append(Spacer(1, 0.2*cm))

        para_rows = [
            [Paragraph(f"<b>{r[0]}</b>", S["body"]),
             Paragraph(r[1], S["small"])]
            for r in rows
        ]
        dt2 = Table(para_rows, colWidths=[3*cm, 14*cm])
        dt2.setStyle(TableStyle([
            ("ROWBACKGROUNDS", (0,0), (-1,-1), [_LLBLUE, colors.white]),
            ("GRID",    (0,0), (-1,-1), 0.3, _LGRAY),
            ("PADDING", (0,0), (-1,-1), 4),
            ("VALIGN",  (0,0), (-1,-1), "TOP"),
            ("BACKGROUND", (0,-1), (-1,-1), colors.HexColor("#ddeeff")),
        ]))
        story.append(dt2)
        story.append(Spacer(1, 0.25*cm))
        if i % 4 == 0 and i < len(lab_candidates):
            story.append(PageBreak())

    # ── Full Ranking (top 20) ─────────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph("All Candidates — Top 20 Ranking", S["h1"]))

    fh = ["Rank", "ID", "Name", "Vina", "ADMET", "Mfg", "Novelty", "Comp.", "Lab"]
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
            "YES" if is_lab else "",
        ])
        if is_lab:
            hl.append(("BACKGROUND", (0,rank), (-1,rank), _GREEN))

    ft = Table(fr, colWidths=[1*cm,1.5*cm,4.2*cm,1.8*cm,1.8*cm,1.8*cm,1.8*cm,1.8*cm,1.5*cm],
               repeatRows=1)
    ft.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), _NAVY),
        ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
        ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE",   (0,0), (-1,-1), 7.5),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, _LLBLUE]),
        ("GRID",  (0,0), (-1,-1), 0.5, _LGRAY),
        ("ALIGN", (3,0), (-1,-1), "CENTER"),
        ("PADDING", (0,0), (-1,-1), 3.5),
        *hl,
    ]))
    story.append(ft)

    # ── Validation Section ────────────────────────────────────────────────────
    if validation_results:
        story.append(PageBreak())
        story.append(Paragraph("Validation: Known EGFR Inhibitors", S["h1"]))
        story.append(Paragraph(
            "Can the platform discriminate FDA-approved EGFR inhibitors from negative control?",
            S["body"]
        ))
        story.append(Spacer(1, 0.3*cm))

        val_hdr = ["Drug Name", "IC50 (nM)", "Vina (kcal/mol)", "Binding/10", "Type"]
        val_rows = [val_hdr]
        pos_scores = []
        neg_scores = []

        for r in sorted(validation_results, key=lambda x: x.get("ic50_nm", 9999)):
            is_neg = r.get("neg_control", False)
            vina   = r.get("vina_score_kcal", 0)
            bind   = r.get("binding_score", 0)
            tip    = "NEG Control" if is_neg else "FDA Approved"
            val_rows.append([
                r["name"][:35],
                f"{r.get('ic50_nm', 0):.1f}",
                f"{vina:.2f}" if r.get("success") else "FAILED",
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
            ("GRID",    (0,0), (-1,-1), 0.5, _LGRAY),
            ("ALIGN",   (1,0), (-1,-1), "CENTER"),
            ("PADDING", (0,0), (-1,-1), 5),
        ]
        # Highlight negative control in yellow
        for ri, r in enumerate(validation_results, 1):
            if r.get("neg_control"):
                val_style.append(("BACKGROUND", (0,ri), (-1,ri), colors.HexColor("#fff3cd")))
        val_tbl.setStyle(TableStyle(val_style))
        story.append(val_tbl)
        story.append(Spacer(1, 0.4*cm))

        # Result summary
        if pos_scores and neg_scores:
            avg_pos = sum(pos_scores) / len(pos_scores)
            avg_neg = sum(neg_scores) / len(neg_scores)
            if avg_pos < avg_neg:
                result_text = (f"SUCCESS: FDA-approved drugs avg {avg_pos:.1f} kcal/mol, "
                               f"negative control {avg_neg:.1f} kcal/mol. "
                               f"Platform correctly discriminates EGFR inhibitors.")
                result_color = colors.HexColor("#d4edda")
            else:
                result_text = (f"WARNING: Negative control ({avg_neg:.1f}) scored similarly "
                               f"to approved drugs ({avg_pos:.1f}). Review binding site coordinates.")
                result_color = colors.HexColor("#fff3cd")

            result_tbl = Table([[Paragraph(result_text, S["body"])]], colWidths=[17.5*cm])
            result_tbl.setStyle(TableStyle([
                ("BACKGROUND", (0,0), (-1,-1), result_color),
                ("GRID",    (0,0), (-1,-1), 0.5, _LGRAY),
                ("PADDING", (0,0), (-1,-1), 8),
            ]))
            story.append(result_tbl)
            story.append(Spacer(1, 0.3*cm))

        # Methodological note — always show
        note_text = (
            "Methodological Note: Rigid docking scores are influenced by molecular size. "
            "Relative ranking within the same scaffold class is meaningful. "
            "For absolute IC50 correlation, MM-GBSA rescoring or FEP calculations are recommended."
        )
        note_tbl = Table([[Paragraph(note_text, S["small"])]], colWidths=[17.5*cm])
        note_tbl.setStyle(TableStyle([
            ("BACKGROUND",   (0,0), (-1,-1), colors.HexColor("#f8f9fa")),
            ("GRID",         (0,0), (-1,-1), 0.5, _LGRAY),
            ("PADDING",      (0,0), (-1,-1), 8),
            ("LEFTPADDING",  (0,0), (-1,-1), 10),
        ]))
        story.append(note_tbl)

    # Footer
    story.append(Spacer(1, 0.5*cm))
    story.append(HRFlowable(width="100%", thickness=1, color=_LGRAY))
    story.append(Spacer(1, 0.2*cm))
    story.append(Paragraph(
        "REINVENT 4 generative model | AutoDock Vina molecular docking | "
        "admet-ai ADMET prediction | RDKit SA Score synthesizability. "
        f"Generated: {now}",
        S["footer"]
    ))
    doc.build(story)
    return path
