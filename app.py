"""
Ilac Kesif Platformu v2 - Gradio Masaustu Arayuzu
Kullanim: python app.py
"""
import gradio as gr
import os
import sys
import json
import glob
import shutil
import threading
import time
from datetime import datetime

# Platform klasorunu ekle
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
os.chdir(BASE_DIR)

# Hugging Face / Linux ortamında otomatik kurulum
try:
    import setup_hf
    setup_hf.setup()
except Exception as e:
    print(f"Setup atlandı: {e}")

# Dil sozlugu
LANG = {
    "tr": {
        "title": "Ilac Kesif Platformu v2",
        "subtitle": "REINVENT 4 | AutoDock Vina | admet-ai | RDKit",
        "tab_run": "Tarama",
        "tab_results": "Sonuclar",
        "tab_settings": "Ayarlar",
        "tab_about": "Hakkinda",
        "protein_name": "Hedef Protein Adi",
        "protein_ph": "ornek: EGFR, ACE2, BRAF",
        "pdb_file": "Protein PDB Dosyasi",
        "coord_x": "Binding Site X",
        "coord_y": "Binding Site Y",
        "coord_z": "Binding Site Z",
        "box_size": "Arama Kutusu Boyutu (A)",
        "rl_mode": "RL Modu (EGFR odakli)",
        "rl_steps": "RL Adim Sayisi",
        "top_n": "Lab icin En Iyi N Aday",
        "start_btn": "Taramayi Baslat",
        "stop_btn": "Durdur",
        "download_pdf": "PDF Raporu Indir",
        "download_csv": "CSV Indir",
        "log_title": "Islem Logu",
        "progress_title": "Ilerleme",
        "results_title": "Sonuclar",
        "status_idle": "Hazir",
        "status_running": "Calisıyor...",
        "status_done": "Tamamlandi",
        "status_error": "Hata",
        "step1": "1/5 REINVENT 4 ile SMILES uretimi",
        "step2": "2/5 Receptor hazırlama",
        "step3": "3/5 AutoDock Vina docking",
        "step4": "4/5 ADMET degerlendirme",
        "step5": "5/5 Skorlama ve rapor",
        "stepval": "Validasyon",
        "no_results": "Henuz sonuc yok. Taramayi baslatin.",
        "settings_weights": "Agirlik Ayarlari",
        "binding_w": "Baglanma Agirligi",
        "admet_w": "ADMET Agirligi",
        "mfg_w": "Uretilebilirlik Agirligi",
        "novelty_w": "Yenilik Agirligi",
        "about_text": """
## Ilac Kesif Platformu v2

Bu platform, hedef proteine yonelik yeni ilac adaylarini otomatik olarak uretir,
degerlendirir ve en iyi adaylari secer.

### Araclar
- **REINVENT 4** — AstraZeneca'nin generatif molekul tasarim araci
- **AutoDock Vina** — Molekuler docking
- **admet-ai** — ADMET tahmin modelleri (50+ ozellik)
- **RDKit** — Kimyasal bilgisayim

### Is Akisi
1. REINVENT 4 ile SMILES uretimi
2. Receptor PDB hazırlama (pdbfixer + ADFR)
3. AutoDock Vina ile docking
4. admet-ai ile ADMET degerlendirme
5. Agirlikli skorlama ve PDF raporu

### Referanslar
- Eberhardt et al., J. Chem. Inf. Model. (2021) - AutoDock Vina
- Loeffler et al., J. Cheminform. (2024) - REINVENT 4
        """,
    },
    "en": {
        "title": "Drug Discovery Platform v2",
        "subtitle": "REINVENT 4 | AutoDock Vina | admet-ai | RDKit",
        "tab_run": "Run",
        "tab_results": "Results",
        "tab_settings": "Settings",
        "tab_about": "About",
        "protein_name": "Target Protein Name",
        "protein_ph": "e.g. EGFR, ACE2, BRAF",
        "pdb_file": "Protein PDB File",
        "coord_x": "Binding Site X",
        "coord_y": "Binding Site Y",
        "coord_z": "Binding Site Z",
        "box_size": "Search Box Size (A)",
        "rl_mode": "RL Mode (EGFR-focused)",
        "rl_steps": "RL Steps",
        "top_n": "Top N Candidates for Lab",
        "start_btn": "Start Screening",
        "stop_btn": "Stop",
        "download_pdf": "Download PDF Report",
        "download_csv": "Download CSV",
        "log_title": "Process Log",
        "progress_title": "Progress",
        "results_title": "Results",
        "status_idle": "Ready",
        "status_running": "Running...",
        "status_done": "Completed",
        "status_error": "Error",
        "step1": "1/5 SMILES generation with REINVENT 4",
        "step2": "2/5 Receptor preparation",
        "step3": "3/5 AutoDock Vina docking",
        "step4": "4/5 ADMET evaluation",
        "step5": "5/5 Scoring and report",
        "stepval": "Validation",
        "no_results": "No results yet. Start a screening run.",
        "settings_weights": "Weight Settings",
        "binding_w": "Binding Weight",
        "admet_w": "ADMET Weight",
        "mfg_w": "Manufacturability Weight",
        "novelty_w": "Novelty Weight",
        "about_text": """
## Drug Discovery Platform v2

This platform automatically generates, evaluates, and selects the best
drug candidates targeting a specific protein.

### Tools
- **REINVENT 4** — AstraZeneca's generative molecule design tool
- **AutoDock Vina** — Molecular docking
- **admet-ai** — ADMET prediction models (50+ properties)
- **RDKit** — Cheminformatics

### Workflow
1. SMILES generation with REINVENT 4
2. Receptor PDB preparation (pdbfixer + ADFR)
3. Docking with AutoDock Vina
4. ADMET evaluation with admet-ai
5. Weighted scoring and PDF report

### References
- Eberhardt et al., J. Chem. Inf. Model. (2021) - AutoDock Vina
- Loeffler et al., J. Cheminform. (2024) - REINVENT 4
        """,
    }
}

# Global state
current_lang = "tr"
run_state = {"running": False, "log": [], "progress": 0, "step": "", "results": None, "pdf": None, "csv": None}


def t(key):
    return LANG[current_lang].get(key, key)


def run_pipeline(protein_name, pdb_file, cx, cy, cz, box_size, rl_mode, rl_steps, top_n,
                 b_weight, a_weight, m_weight, n_weight):
    """Pipeline'i ayri thread'de calistirir, log ve progress gunceller."""
    global run_state

    run_state["running"] = True
    run_state["log"] = []
    run_state["progress"] = 0
    run_state["results"] = None
    run_state["pdf"] = None
    run_state["csv"] = None

    def log(msg):
        ts = datetime.now().strftime("%H:%M:%S")
        run_state["log"].append(f"[{ts}] {msg}")

    try:
        import config

        # PDB dosyasını kopyala
        if pdb_file is not None:
            dest = os.path.join(BASE_DIR, "receptor", "protein.pdb")
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copy2(pdb_file, dest)
            config.RECEPTOR_PDB_PATH = dest
            # Eski PDBQT'yi sil
            pdbqt = os.path.join(BASE_DIR, "receptor", "receptor.pdbqt")
            if os.path.exists(pdbqt):
                os.remove(pdbqt)
            log(f"PDB dosyasi yuklendi: {os.path.basename(pdb_file)}")
        else:
            config.RECEPTOR_PDB_PATH = os.path.join(BASE_DIR, "receptor", "protein.pdb")

        # Config guncelle
        config.BINDING_SITE = {
            "center_x": float(cx), "center_y": float(cy), "center_z": float(cz),
            "size_x": float(box_size), "size_y": float(box_size), "size_z": float(box_size),
        }
        config.USE_RL_MODE = rl_mode
        config.RL_MAX_STEPS = int(rl_steps)
        config.TOP_N_FOR_LAB = int(top_n)
        config.WEIGHTS = {
            "binding_score": float(b_weight),
            "admet_score": float(a_weight),
            "manufacturability_score": float(m_weight),
            "novelty_score": float(n_weight),
        }

        # ── ADIM 1 ────────────────────────────────────────────────────────────
        run_state["step"] = t("step1")
        run_state["progress"] = 5
        log(f"ADIM 1: SMILES uretimi baslıyor ({'RL modu' if rl_mode else 'Sampling modu'})...")

        import reinvent_runner, reinvent_rl_runner
        if rl_mode:
            candidates = reinvent_rl_runner.generate_candidates_rl()
        else:
            candidates = reinvent_runner.generate_candidates()
        log(f"  {len(candidates)} aday uretildi.")
        run_state["progress"] = 20

        # ── ADIM 2 ────────────────────────────────────────────────────────────
        run_state["step"] = t("step2")
        log("ADIM 2: Receptor hazirlaniyor...")
        import vina_runner
        receptor_pdbqt = vina_runner.prepare_receptor()
        log(f"  Receptor hazir: {receptor_pdbqt}")
        run_state["progress"] = 30

        # ── ADIM 3 ────────────────────────────────────────────────────────────
        run_state["step"] = t("step3")
        log(f"ADIM 3: Docking ({len(candidates)} aday)...")
        docked = []
        total = len(candidates)
        for i, c in enumerate(candidates):
            r = vina_runner.dock_single(c["smiles"], receptor_pdbqt)
            docked.append({**c, **r})
            prog = 30 + int((i + 1) / total * 30)
            run_state["progress"] = prog
            if r["success"]:
                log(f"  {c['id']}: {r['vina_score_kcal']:.1f} kcal/mol")
        success = sum(1 for d in docked if d.get("success"))
        log(f"  Docking: {success}/{total} basarili.")
        run_state["progress"] = 60

        # ── ADIM 4 ────────────────────────────────────────────────────────────
        run_state["step"] = t("step4")
        log("ADIM 4: ADMET degerlendirme...")
        import admet_client
        admet_scored = admet_client.evaluate_all(docked)
        log("  ADMET tamamlandi.")
        run_state["progress"] = 75

        # ── ADIM 5 ────────────────────────────────────────────────────────────
        run_state["step"] = t("step5")
        log("ADIM 5: Skorlama ve rapor olusturuluyor...")
        import scorer, reporter
        all_scored = scorer.compute_scores(admet_scored)
        lab_candidates = scorer.select_lab_candidates(all_scored)
        stats = scorer.get_statistics(all_scored)
        log(f"  {len(all_scored)} aday sıralandı | Ort: {stats.get('avg',0):.2f} | Lab: {len(lab_candidates)}")
        run_state["progress"] = 85

        # ── VALİDASYON ────────────────────────────────────────────────────────
        run_state["step"] = t("stepval")
        log("Validasyon: referans ilaclar test ediliyor...")
        from main import REFERENCE_DRUGS
        validation_results = []
        for drug in REFERENCE_DRUGS:
            r = vina_runner.dock_single(drug["smiles"], receptor_pdbqt)
            validation_results.append({**drug, **r})
            log(f"  {drug['name'][:25]}: {r['vina_score_kcal']:.2f} kcal/mol")
        run_state["progress"] = 92

        # ── RAPORLAR ──────────────────────────────────────────────────────────
        output_dir = os.path.join(BASE_DIR, "results")
        os.makedirs(output_dir, exist_ok=True)

        pdf_path = reporter.save_pdf(protein_name, all_scored, lab_candidates,
                                     output_dir, validation_results)
        csv_path = reporter.save_csv(all_scored, lab_candidates, output_dir, protein_name)
        log(f"  PDF: {os.path.basename(pdf_path)}")
        log(f"  CSV: {os.path.basename(csv_path)}")

        run_state["pdf"] = pdf_path
        run_state["csv"] = csv_path
        run_state["results"] = {
            "total": len(all_scored),
            "lab": len(lab_candidates),
            "avg": stats.get("avg", 0),
            "best": lab_candidates[0] if lab_candidates else None,
            "lab_candidates": lab_candidates,
        }
        run_state["progress"] = 100
        run_state["step"] = t("status_done")
        log("=" * 50)
        log(f"TAMAMLANDI — {len(lab_candidates)} aday laboratuvara hazir.")

    except Exception as e:
        log(f"HATA: {e}")
        run_state["step"] = t("status_error")
        run_state["progress"] = 0
        import traceback
        log(traceback.format_exc()[-500:])
    finally:
        run_state["running"] = False


def start_run(protein_name, pdb_file, cx, cy, cz, box_size, rl_mode, rl_steps, top_n,
              b_weight, a_weight, m_weight, n_weight):
    if run_state["running"]:
        return "Zaten calisıyor!", 0, "", gr.update(value=None), gr.update(value=None)

    if not protein_name:
        return "HATA: Protein adi giriniz.", 0, "", gr.update(value=None), gr.update(value=None)

    t_args = (protein_name, pdb_file, cx, cy, cz, box_size, rl_mode, rl_steps, top_n,
              b_weight, a_weight, m_weight, n_weight)
    thread = threading.Thread(target=run_pipeline, args=t_args, daemon=True)
    thread.start()
    return "Basladi...", 0, "", gr.update(value=None), gr.update(value=None)


def poll_status():
    """Durum guncelleme - her saniye cagirilir."""
    log_text = "\n".join(run_state["log"][-50:])
    progress = run_state["progress"]
    step = run_state["step"]

    results_html = ""
    if run_state["results"]:
        r = run_state["results"]
        best = r.get("best")
        results_html = f"""
<div style='font-family:monospace; padding:12px; background:#1a1a2e; color:#e0e0e0; border-radius:8px;'>
<h3 style='color:#00d4aa;'>Tarama Tamamlandi</h3>
<p>Taranan: <b>{r['total']}</b> &nbsp;|&nbsp; Lab secimi: <b>{r['lab']}</b> &nbsp;|&nbsp; Ort. skor: <b>{r['avg']:.2f}/10</b></p>
"""
        if best:
            results_html += f"""
<hr style='border-color:#333;'>
<h4 style='color:#ffd700;'>En Iyi Aday: {best['id']}</h4>
<p>SMILES: <code style='font-size:11px;'>{best['smiles'][:60]}...</code></p>
<p>Vina: <b>{best.get('vina_score_kcal',0):.2f} kcal/mol</b> &nbsp;|&nbsp;
   ADMET: <b>{best.get('admet_score',0):.1f}/10</b> &nbsp;|&nbsp;
   Kompozit: <b>{best.get('composite_score',0):.2f}/10</b></p>
"""
        # Top 10 table
        lab = r.get("lab_candidates", [])
        if lab:
            results_html += "<hr style='border-color:#333;'><h4>Lab Adaylari</h4>"
            results_html += "<table style='width:100%; font-size:12px; border-collapse:collapse;'>"
            results_html += "<tr style='color:#00d4aa;'><th>#</th><th>ID</th><th>Vina</th><th>ADMET</th><th>Kompozit</th></tr>"
            for i, c in enumerate(lab, 1):
                results_html += f"<tr><td>{i}</td><td>{c['id']}</td><td>{c.get('vina_score_kcal',0):.1f}</td><td>{c.get('admet_score',0):.1f}</td><td>{c.get('composite_score',0):.2f}</td></tr>"
            results_html += "</table>"
        results_html += "</div>"
    else:
        results_html = f"<p style='color:#888;'>{t('no_results')}</p>"

    pdf_update = gr.update(value=run_state["pdf"], visible=run_state["pdf"] is not None)
    csv_update = gr.update(value=run_state["csv"], visible=run_state["csv"] is not None)

    return log_text, progress, step, results_html, pdf_update, csv_update


def change_lang(lang):
    global current_lang
    current_lang = "en" if lang == "English" else "tr"
    return gr.update(label=t("protein_name"))


# ── CSS ───────────────────────────────────────────────────────────────────────
CSS = """
body { font-family: 'Segoe UI', sans-serif; }
.gradio-container { max-width: 1100px !important; margin: 0 auto; }
#title { text-align: center; padding: 20px 0 5px; }
#title h1 { font-size: 2em; font-weight: 700; color: #1a3a5c; margin:0; }
#title p { color: #666; margin: 4px 0 0; }
#start_btn { background: #1a3a5c !important; color: white !important; font-size: 1.1em !important; }
#start_btn:hover { background: #2c5282 !important; }
.progress-bar { height: 8px; border-radius: 4px; }
#log_box textarea { font-family: monospace; font-size: 12px; background: #0d1117; color: #c9d1d9; }
"""

# ── ARAYUZ ────────────────────────────────────────────────────────────────────
with gr.Blocks(css=CSS, title="Drug Discovery Platform v2", theme=gr.themes.Soft()) as app:

    # Baslik
    with gr.Row(elem_id="title"):
        gr.HTML("""
        <div style='text-align:center; padding:20px 0 10px;'>
            <h1 style='font-size:2em; font-weight:700; color:#1a3a5c; margin:0;'>
                Ilac Kesif Platformu v2
            </h1>
            <p style='color:#666; margin:4px 0 0;'>
                REINVENT 4 &nbsp;|&nbsp; AutoDock Vina &nbsp;|&nbsp; admet-ai &nbsp;|&nbsp; RDKit
            </p>
        </div>
        """)

    # Dil secici
    with gr.Row():
        gr.HTML("<div></div>")
        lang_radio = gr.Radio(["Turkce", "English"], value="Turkce", label="Dil / Language", scale=0)

    with gr.Tabs():

        # ── TARAMA SEKMESI ────────────────────────────────────────────────────
        with gr.TabItem("Tarama / Run"):
            with gr.Row():
                with gr.Column(scale=1):
                    gr.Markdown("### Protein Ayarlari")
                    protein_name = gr.Textbox(label="Hedef Protein Adi", placeholder="ornek: EGFR, ACE2, BRAF")
                    pdb_file = gr.File(label="Protein PDB Dosyasi (.pdb)", file_types=[".pdb"])

                    gr.Markdown("### Binding Site Koordinatlari")
                    with gr.Row():
                        cx = gr.Number(label="X", value=26.5)
                        cy = gr.Number(label="Y", value=11.5)
                        cz = gr.Number(label="Z", value=-1.0)
                    box_size = gr.Slider(10, 40, value=25, step=1, label="Kutu Boyutu (A)")

                    gr.Markdown("### Uretim Ayarlari")
                    rl_mode = gr.Checkbox(label="RL Modu (EGFR odakli)", value=True)
                    rl_steps = gr.Slider(10, 500, value=50, step=10, label="RL Adim Sayisi")
                    top_n = gr.Slider(5, 20, value=10, step=1, label="Lab icin En Iyi N Aday")

                    start_btn = gr.Button("Taramayi Baslat", variant="primary", elem_id="start_btn")

                with gr.Column(scale=2):
                    gr.Markdown("### Ilerleme")
                    step_label = gr.Textbox(label="Mevcut Adim", interactive=False, value="Hazir")
                    progress_bar = gr.Slider(0, 100, value=0, label="Ilerleme (%)", interactive=False)

                    gr.Markdown("### Islem Logu")
                    log_box = gr.Textbox(
                        label="Log",
                        lines=18,
                        max_lines=18,
                        interactive=False,
                        elem_id="log_box",
                        placeholder="Tarama baslayinca buraya log yazilacak..."
                    )

        # ── SONUCLAR SEKMESI ──────────────────────────────────────────────────
        with gr.TabItem("Sonuclar / Results"):
            results_html = gr.HTML("<p style='color:#888; padding:20px;'>Henuz sonuc yok. Taramayi baslatin.</p>")
            with gr.Row():
                pdf_download = gr.File(label="PDF Raporu Indir", visible=False)
                csv_download = gr.File(label="CSV Indir", visible=False)

        # ── AYARLAR SEKMESI ───────────────────────────────────────────────────
        with gr.TabItem("Ayarlar / Settings"):
            gr.Markdown("### Agirlik Ayarlari (Toplam = 1.0)")
            with gr.Row():
                b_weight = gr.Slider(0.0, 1.0, value=0.40, step=0.05, label="Baglanma (Vina)")
                a_weight = gr.Slider(0.0, 1.0, value=0.35, step=0.05, label="ADMET")
            with gr.Row():
                m_weight = gr.Slider(0.0, 1.0, value=0.15, step=0.05, label="Uretilebilirlik")
                n_weight = gr.Slider(0.0, 1.0, value=0.10, step=0.05, label="Yenilik")

            gr.Markdown("""
            > **Not:** Agirlik toplami 1.0 olmalidir.
            > Kanser arastirmasi icin ADMET agirligini artirmaniz onerilir.
            """)

        # ── HAKKINDA SEKMESI ──────────────────────────────────────────────────
        with gr.TabItem("Hakkinda / About"):
            gr.Markdown("""
## Ilac Kesif Platformu v2

Bu platform, hedef proteine yonelik yeni ilac adaylarini otomatik olarak uretir,
degerlendirir ve en iyi adaylari secer.

### Araclar
- **REINVENT 4** — AstraZeneca'nin generatif molekul tasarim araci (ChEMBL'de egitilmis)
- **AutoDock Vina** — Scripps Research'in molekuler docking araci
- **admet-ai** — 50+ ADMET ozellik tahmini (hERG, DILI, AMES, biyoyararlilik...)
- **RDKit** — SA Score, Tanimoto benzerlik, Lipinski

### Is Akisi
1. REINVENT 4 ile 50-60 SMILES uretimi (Sampling veya RL modu)
2. Receptor PDB hazirlama (pdbfixer + ADFR)
3. AutoDock Vina ile molekuler docking
4. admet-ai ile ADMET degerlendirme
5. Agirlikli kompozit skorlama
6. PDF + CSV raporu

### Skorlama Formulu
```
Kompozit = Baglanma×0.40 + ADMET×0.35 + Uretim×0.15 + Yenilik×0.10
```

### Referanslar
- Eberhardt et al., J. Chem. Inf. Model. 61, 3891 (2021)
- Loeffler et al., J. Cheminform. 16, 20 (2024)
- Dong et al., J. Chem. Inf. Model. 64, 3 (2024) - admet-ai
            """)

    # ── EVENT HANDLER'LAR ─────────────────────────────────────────────────────
    start_btn.click(
        fn=start_run,
        inputs=[protein_name, pdb_file, cx, cy, cz, box_size,
                rl_mode, rl_steps, top_n, b_weight, a_weight, m_weight, n_weight],
        outputs=[step_label, progress_bar, log_box, pdf_download, csv_download],
    )

    # Her 2 saniyede otomatik guncelle - Gradio 5.x uyumlu
    timer = gr.Timer(value=2)
    timer.tick(
        fn=poll_status,
        outputs=[log_box, progress_bar, step_label, results_html, pdf_download, csv_download],
    )


if __name__ == "__main__":
    print("\n" + "="*60)
    print("  Ilac Kesif Platformu v2 - Masaustu Arayuzu")
    print("  Tarayicinizda acilacak: http://localhost:7860")
    print("="*60 + "\n")
    app.launch(
        server_name="0.0.0.0",
        server_port=7860,
        inbrowser=True,
        share=False,
    )
