"""
Drug Discovery Platform v2 - REINVENT 4 RL Mode
EGFR-focused SMILES generation using reinforcement learning.
"""
import os
import csv
import glob
import subprocess
import tempfile
import sys
import config

EGFR_REFERENCE_SMILES = [
    "COCCOc1cc2ncnc(Nc3cccc(Cl)c3F)c2cc1OCCOC",    # Erlotinib
    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",  # Gefitinib
    "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCC1CCN(C)CC1",  # Afatinib
]


def _write_rl_toml(toml_path, output_csv, agent_path):
    prior  = os.path.abspath(config.REINVENT_PRIOR_PATH).replace("\\", "/")
    agent  = os.path.abspath(agent_path).replace("\\", "/")
    output = os.path.abspath(output_csv).replace("\\", "/")
    chkpt  = agent.replace(".prior", "_checkpoint.chkpt")
    ref_smiles_lines = ",\n            ".join(f'"{s}"' for s in EGFR_REFERENCE_SMILES)

    toml_text = f'''run_type = "staged_learning"
device   = "cpu"

[parameters]
prior_file   = "{prior}"
agent_file   = "{agent}"
summary_csv_prefix = "{output.replace('.csv', '')}"
batch_size   = 64
randomize_smiles = true
unique_sequences = true

[learning_strategy]
type = "dap"
sigma = 128
rate = 0.0001

[diversity_filter]
type = "IdenticalMurckoScaffold"
bucket_size = 25
minscore = 0.3

[[stage]]
max_steps = {config.RL_MAX_STEPS}
max_score = 0.8
min_steps = 10
chkpt_file = "{chkpt}"

    [stage.scoring]
    type = "geometric_mean"

        [[stage.scoring.component]]
        [stage.scoring.component.TanimotoSimilarity]
            [[stage.scoring.component.TanimotoSimilarity.endpoint]]
            name = "EGFR_similarity"
            weight = 0.4
            params.smiles = [
            {ref_smiles_lines}
            ]
            params.radius = 2
            params.use_counts = true
            params.use_features = true
            transform.type = "sigmoid"
            transform.high = 1.0
            transform.low = 0.0
            transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.QED]
            [[stage.scoring.component.QED.endpoint]]
            name = "drug_likeness"
            weight = 0.3
            transform.type = "sigmoid"
            transform.high = 1.0
            transform.low = 0.0
            transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.SAScore]
            [[stage.scoring.component.SAScore.endpoint]]
            name = "synthesizability"
            weight = 0.2
            transform.type = "reverse_sigmoid"
            transform.high = 1.0
            transform.low = 0.0
            transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.MolecularWeight]
            [[stage.scoring.component.MolecularWeight.endpoint]]
            name = "mol_weight"
            weight = 0.1
            transform.type = "double_sigmoid"
            transform.high = 1.0
            transform.low = 0.0
            transform.coef_div = 500.0
            transform.coef_si  = 20.0
            transform.coef_se  = 20.0
'''
    with open(toml_path, "w", encoding="utf-8") as f:
        f.write(toml_text)


def _run_rl(toml_path):
    cmd = [
        sys.executable, "-c",
        f"from reinvent.reinvent_main import main_script; "
        f"import sys; sys.argv = ['reinvent', r'{toml_path}']; "
        f"main_script()"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=config.RL_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(
            f"REINVENT RL error (code {result.returncode}):\n"
            f"STDERR: {result.stderr[-3000:]}\n"
            f"STDOUT: {result.stdout[-500:]}"
        )
    return result.returncode


def _parse_rl_output(output_dir):
    candidates = []
    seen = set()

    # Use set to avoid duplicate files
    csv_files = list(set(
        glob.glob(os.path.join(output_dir, "*.csv")) +
        glob.glob(os.path.join(output_dir, "**", "*.csv"), recursive=True)
    ))

    for csv_file in sorted(csv_files):
        try:
            with open(csv_file, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    smiles = (row.get("SMILES") or row.get("smiles") or
                              row.get("Smiles") or "").strip()
                    try:
                        score = float(row.get("Score") or row.get("score") or 0)
                    except (ValueError, TypeError):
                        score = 0.0
                    if smiles and smiles not in seen and score > 0.0:
                        seen.add(smiles)
                        candidates.append({"smiles": smiles, "rl_score": score})
        except Exception as e:
            print(f"  CSV read error: {e}")
            continue

    candidates.sort(key=lambda x: x["rl_score"], reverse=True)
    top = candidates[:config.REINVENT_NUM_SMILES]

    result = []
    for i, c in enumerate(top):
        result.append({
            "id":        f"RL{i+1:03d}",
            "smiles":    c["smiles"],
            "name":      f"RL-Candidate-{i+1:03d}",
            "mechanism": f"REINVENT RL EGFR-focused (score={c['rl_score']:.3f})",
        })
    return result


def generate_candidates_rl():
    if not os.path.exists(config.REINVENT_PRIOR_PATH):
        raise FileNotFoundError(f"Prior model not found: {config.REINVENT_PRIOR_PATH}")

    agent_path = config.REINVENT_AGENT_PATH
    if not os.path.exists(agent_path):
        import shutil
        print(f"  Agent file not found, copying from prior: {agent_path}")
        os.makedirs(os.path.dirname(os.path.abspath(agent_path)), exist_ok=True)
        shutil.copy2(config.REINVENT_PRIOR_PATH, agent_path)
    else:
        print(f"  Using existing agent: {agent_path}")

    with tempfile.TemporaryDirectory() as tmpdir:
        toml_path  = os.path.join(tmpdir, "rl_config.toml")
        output_csv = os.path.join(tmpdir, "rl_output.csv")

        _write_rl_toml(toml_path, output_csv, agent_path)

        print(f"  Starting RL training ({config.RL_MAX_STEPS} steps)...")
        _run_rl(toml_path)

        candidates = _parse_rl_output(tmpdir)

    if not candidates:
        print("  RL output empty, falling back to sampling mode...")
        from reinvent_runner import generate_candidates
        return generate_candidates()

    print(f"  RL mode: {len(candidates)} EGFR-focused SMILES generated.")
    return candidates
