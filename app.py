"""
Drug Discovery Platform v2 - Gradio Desktop Interface
Usage: python app.py
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

# Add platform directory to path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
os.chdir(BASE_DIR)

# Auto-setup for Hugging Face / Linux
try:
    import setup_hf
    setup_hf.setup()
except Exception as e:
    print(f"Setup skipped: {e}")

# Global state
current_lang = "en"
run_state = {"running": False, "log": [], "progress": 0, "step": "", "results": None, "pdf": None, "csv": None}


def run_pipeline(protein_name, pdb_file, cx, cy, cz, box_size, rl_mode, rl_steps, top_n,
                 b_weight, a_weight, m_weight, n_weight):
    """Runs the pipeline in a separate thread, updates log and progress."""
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

        # Copy PDB file
        if pdb_file is not None:
            dest = os.path.join(BASE_DIR, "receptor", "protein.pdb")
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copy2(pdb_file, dest)
            config.RECEPTOR_PDB_PATH = dest
            # Remove old PDBQT
            pdbqt = os.path.join(BASE_DIR, "receptor", "receptor.pdbqt")
            if os.path.exists(pdbqt):
                os.remove(pdbqt)
            log(f"PDB file loaded: {os.path.basename(pdb_file)}")
        else:
            config.RECEPTOR_PDB_PATH = os.path.join(BASE_DIR, "receptor", "protein.pdb")

        # Update config
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

        # ── Step 1 ────────────────────────────────────────────────────────────
        run_state["step"] = "1/5 SMILES Generation"
        run_state["progress"] = 5
        log(f"STEP 1: SMILES generation starting ({'RL mode' if rl_mode else 'Sampling mode'})...")

        import reinvent_runner, reinvent_rl_runner
        if rl_mode:
            candidates = reinvent_rl_runner.generate_candidates_rl()
        else:
            candidates = reinvent_runner.generate_candidates()
        log(f"  {len(candidates)} candidates generated.")
        run_state["progress"] = 20

        # ── Step 2 ────────────────────────────────────────────────────────────
        run_state["step"] = "2/5 Receptor Preparation"
        log("STEP 2: Preparing receptor...")
        import vina_runner
        receptor_pdbqt = vina_runner.prepare_receptor()
        log(f"  Receptor ready: {receptor_pdbqt}")
        run_state["progress"] = 30

        # ── Step 3 ────────────────────────────────────────────────────────────
        run_state["step"] = "3/5 Molecular Docking"
        log(f"STEP 3: Docking ({len(candidates)} candidates)...")
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
        log(f"  Docking: {success}/{total} successful.")
        run_state["progress"] = 60

        # ── Step 4 ────────────────────────────────────────────────────────────
        run_state["step"] = "4/5 ADMET Evaluation"
        log("STEP 4: ADMET evaluation...")
        import admet_client
        admet_scored = admet_client.evaluate_all(docked)
        log("  ADMET complete.")
        run_state["progress"] = 75

        # ── Step 5 ────────────────────────────────────────────────────────────
        run_state["step"] = "5/5 Scoring and Report"
        log("STEP 5: Scoring and generating report...")
        import scorer, reporter
        all_scored = scorer.compute_scores(admet_scored)
        lab_candidates = scorer.select_lab_candidates(all_scored)
        stats = scorer.get_statistics(all_scored)
        log(f"  {len(all_scored)} candidates ranked | Avg: {stats.get('avg',0):.2f} | Lab: {len(lab_candidates)}")
        run_state["progress"] = 85

        # ── Validation ────────────────────────────────────────────────────────
        run_state["step"] = "Validation"
        log("Validation: testing reference drugs...")
        from main import REFERENCE_DRUGS
        validation_results = []
        for drug in REFERENCE_DRUGS:
            r = vina_runner.dock_single(drug["smiles"], receptor_pdbqt)
            validation_results.append({**drug, **r})
            log(f"  {drug['name'][:25]}: {r['vina_score_kcal']:.2f} kcal/mol")
        run_state["progress"] = 92

        # ── Reports ───────────────────────────────────────────────────────────
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
        run_state["step"] = "Completed"
        log("=" * 50)
        log(f"COMPLETE — {len(lab_candidates)} candidates ready for laboratory.")

    except Exception as e:
        log(f"ERROR: {e}")
        run_state["step"] = "Error"
        run_state["progress"] = 0
        import traceback
        log(traceback.format_exc()[-500:])
    finally:
        run_state["running"] = False


def start_run(protein_name, pdb_file, cx, cy, cz, box_size, rl_mode, rl_steps, top_n,
              b_weight, a_weight, m_weight, n_weight):
    if run_state["running"]:
        return "Already running!", 0, "", gr.update(value=None), gr.update(value=None)

    if not protein_name:
        return "ERROR: Please enter a protein name.", 0, "", gr.update(value=None), gr.update(value=None)

    t_args = (protein_name, pdb_file, cx, cy, cz, box_size, rl_mode, rl_steps, top_n,
              b_weight, a_weight, m_weight, n_weight)
    thread = threading.Thread(target=run_pipeline, args=t_args, daemon=True)
    thread.start()
    return "Started...", 0, "", gr.update(value=None), gr.update(value=None)


def poll_status():
    """Status update — called every 2 seconds."""
    log_text = "\n".join(run_state["log"][-50:])
    progress = run_state["progress"]
    step = run_state["step"]

    results_html = ""
    if run_state["results"]:
        r = run_state["results"]
        best = r.get("best")
        results_html = f"""
<div style='font-family:monospace; padding:12px; background:#1a1a2e; color:#e0e0e0; border-radius:8px;'>
<h3 style='color:#00d4aa;'>Screening Complete</h3>
<p>Screened: <b>{r['total']}</b> &nbsp;|&nbsp; Lab selection: <b>{r['lab']}</b> &nbsp;|&nbsp; Avg score: <b>{r['avg']:.2f}/10</b></p>
"""
        if best:
            results_html += f"""
<hr style='border-color:#333;'>
<h4 style='color:#ffd700;'>Best Candidate: {best['id']}</h4>
<p>SMILES: <code style='font-size:11px;'>{best['smiles'][:60]}...</code></p>
<p>Vina: <b>{best.get('vina_score_kcal',0):.2f} kcal/mol</b> &nbsp;|&nbsp;
   ADMET: <b>{best.get('admet_score',0):.1f}/10</b> &nbsp;|&nbsp;
   Composite: <b>{best.get('composite_score',0):.2f}/10</b></p>
"""
        lab = r.get("lab_candidates", [])
        if lab:
            results_html += "<hr style='border-color:#333;'><h4>Lab Candidates</h4>"
            results_html += "<table style='width:100%; font-size:12px; border-collapse:collapse;'>"
            results_html += "<tr style='color:#00d4aa;'><th>#</th><th>ID</th><th>Vina</th><th>ADMET</th><th>Composite</th></tr>"
            for i, c in enumerate(lab, 1):
                results_html += f"<tr><td>{i}</td><td>{c['id']}</td><td>{c.get('vina_score_kcal',0):.1f}</td><td>{c.get('admet_score',0):.1f}</td><td>{c.get('composite_score',0):.2f}</td></tr>"
            results_html += "</table>"
        results_html += "</div>"
    else:
        results_html = "<p style='color:#888;'>No results yet. Start a screening run.</p>"

    pdf_update = gr.update(value=run_state["pdf"], visible=run_state["pdf"] is not None)
    csv_update = gr.update(value=run_state["csv"], visible=run_state["csv"] is not None)

    return log_text, progress, step, results_html, pdf_update, csv_update


# ── CSS ───────────────────────────────────────────────────────────────────────
CSS = """
body { font-family: 'Segoe UI', sans-serif; }
.gradio-container { max-width: 1100px !important; margin: 0 auto; }
#start_btn { background: #1a3a5c !important; color: white !important; font-size: 1.1em !important; }
#start_btn:hover { background: #2c5282 !important; }
#log_box textarea { font-family: monospace; font-size: 12px; background: #0d1117; color: #c9d1d9; }
"""

# ── Interface ─────────────────────────────────────────────────────────────────
with gr.Blocks(css=CSS, title="Drug Discovery Platform v2", theme=gr.themes.Soft()) as app:

    with gr.Row():
        gr.HTML("""
        <div style='text-align:center; padding:20px 0 10px;'>
            <h1 style='font-size:2em; font-weight:700; color:#1a3a5c; margin:0;'>
                Drug Discovery Platform v2
            </h1>
            <p style='color:#666; margin:4px 0 0;'>
                REINVENT 4 &nbsp;|&nbsp; AutoDock Vina &nbsp;|&nbsp; admet-ai &nbsp;|&nbsp; RDKit
            </p>
        </div>
        """)

    with gr.Tabs():

        # ── Run Tab ───────────────────────────────────────────────────────────
        with gr.TabItem("Run"):
            with gr.Row():
                with gr.Column(scale=1):
                    gr.Markdown("### Protein Settings")
                    protein_name = gr.Textbox(label="Target Protein Name", placeholder="e.g. EGFR, ACE2, BRAF")
                    pdb_file = gr.File(label="Protein PDB File (.pdb)", file_types=[".pdb"])

                    gr.Markdown("### Binding Site Coordinates")
                    with gr.Row():
                        cx = gr.Number(label="X", value=26.5)
                        cy = gr.Number(label="Y", value=11.5)
                        cz = gr.Number(label="Z", value=-1.0)
                    box_size = gr.Slider(10, 40, value=25, step=1, label="Box Size (Angstrom)")

                    gr.Markdown("### Generation Settings")
                    rl_mode = gr.Checkbox(label="RL Mode (EGFR-focused)", value=True)
                    rl_steps = gr.Slider(10, 500, value=50, step=10, label="RL Steps")
                    top_n = gr.Slider(5, 20, value=10, step=1, label="Top N Candidates for Lab")

                    start_btn = gr.Button("Start Screening", variant="primary", elem_id="start_btn")

                with gr.Column(scale=2):
                    gr.Markdown("### Progress")
                    step_label = gr.Textbox(label="Current Step", interactive=False, value="Ready")
                    progress_bar = gr.Slider(0, 100, value=0, label="Progress (%)", interactive=False)

                    gr.Markdown("### Process Log")
                    log_box = gr.Textbox(
                        label="Log",
                        lines=18,
                        max_lines=18,
                        interactive=False,
                        elem_id="log_box",
                        placeholder="Log will appear here when screening starts..."
                    )

        # ── Results Tab ───────────────────────────────────────────────────────
        with gr.TabItem("Results"):
            results_html = gr.HTML("<p style='color:#888; padding:20px;'>No results yet. Start a screening run.</p>")
            with gr.Row():
                pdf_download = gr.File(label="Download PDF Report", visible=False)
                csv_download = gr.File(label="Download CSV", visible=False)

        # ── Settings Tab ──────────────────────────────────────────────────────
        with gr.TabItem("Settings"):
            gr.Markdown("### Weight Settings (must sum to 1.0)")
            with gr.Row():
                b_weight = gr.Slider(0.0, 1.0, value=0.40, step=0.05, label="Binding (Vina)")
                a_weight = gr.Slider(0.0, 1.0, value=0.35, step=0.05, label="ADMET")
            with gr.Row():
                m_weight = gr.Slider(0.0, 1.0, value=0.15, step=0.05, label="Manufacturability")
                n_weight = gr.Slider(0.0, 1.0, value=0.10, step=0.05, label="Novelty")

            gr.Markdown("""
            > **Note:** Weights must sum to 1.0.
            > For cancer research, increasing the ADMET weight is recommended.
            """)

        # ── About Tab ─────────────────────────────────────────────────────────
        with gr.TabItem("About"):
            gr.Markdown("""
## Drug Discovery Platform v2

This platform automatically generates, evaluates, and selects the best
drug candidates targeting a specific protein.

### Tools
- **REINVENT 4** — AstraZeneca's generative molecular design tool (trained on ChEMBL)
- **AutoDock Vina** — Scripps Research molecular docking engine
- **admet-ai** — 50+ ADMET property prediction (hERG, DILI, AMES, bioavailability...)
- **RDKit** — SA Score, Tanimoto similarity, Lipinski

### Workflow
1. SMILES generation with REINVENT 4 (Sampling or RL mode)
2. Receptor PDB preparation (pdbfixer + ADFR/obabel)
3. Molecular docking with AutoDock Vina
4. ADMET evaluation with admet-ai
5. Weighted composite scoring
6. PDF + CSV report

### Scoring Formula
```
Composite = Binding×0.40 + ADMET×0.35 + Manufacturability×0.15 + Novelty×0.10
```

### References
- Eberhardt et al., J. Chem. Inf. Model. 61, 3891 (2021) — AutoDock Vina
- Loeffler et al., J. Cheminform. 16, 20 (2024) — REINVENT 4
- Swanson et al., J. Chem. Inf. Model. 64, 2 (2024) — admet-ai
            """)

    # ── Event Handlers ────────────────────────────────────────────────────────
    start_btn.click(
        fn=start_run,
        inputs=[protein_name, pdb_file, cx, cy, cz, box_size,
                rl_mode, rl_steps, top_n, b_weight, a_weight, m_weight, n_weight],
        outputs=[step_label, progress_bar, log_box, pdf_download, csv_download],
    )

    # Auto-update every 2 seconds — Gradio 5.x compatible
    timer = gr.Timer(value=2)
    timer.tick(
        fn=poll_status,
        outputs=[log_box, progress_bar, step_label, results_html, pdf_download, csv_download],
    )


if __name__ == "__main__":
    print("\n" + "="*60)
    print("  Drug Discovery Platform v2 - Desktop Interface")
    print("  Opening in browser: http://localhost:7860")
    print("="*60 + "\n")
    app.launch(
        server_name="0.0.0.0",
        server_port=7860,
        inbrowser=True,
        share=False,
    )
