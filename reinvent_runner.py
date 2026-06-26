"""
İlaç Keşif Platformu v2 — REINVENT 4 SMILES Üreticisi
"""
import os
import csv
import subprocess
import tempfile
import sys
import config


def _write_toml(toml_path: str, output_csv: str):
    prior  = os.path.abspath(config.REINVENT_PRIOR_PATH).replace("\\", "/")
    output = os.path.abspath(output_csv).replace("\\", "/")
    content = f"""run_type = "sampling"
device   = "cpu"

[parameters]
model_file       = "{prior}"
output_file      = "{output}"
num_smiles       = {config.REINVENT_NUM_SMILES}
unique_molecules = true
randomize_smiles = true
"""
    with open(toml_path, "w", encoding="utf-8") as f:
        f.write(content)


def _run_reinvent(toml_path: str) -> int:
    cmd = [
        sys.executable, "-c",
        f"from reinvent.reinvent_main import main_script; "
        f"import sys; sys.argv = ['reinvent', r'{toml_path}']; "
        f"main_script()"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise RuntimeError(
            f"REINVENT 4 hatası (kod {result.returncode}):\n"
            f"STDERR: {result.stderr[-3000:]}\n"
            f"STDOUT: {result.stdout[-1000:]}"
        )
    return result.returncode


def _parse_output(output_csv: str) -> list:
    candidates = []
    with open(output_csv, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        smiles_col = 0
        if header:
            for i, h in enumerate(header):
                if "smiles" in h.lower():
                    smiles_col = i
                    break
        for i, row in enumerate(reader):
            if not row:
                continue
            smiles = row[smiles_col].strip()
            if smiles:
                candidates.append({
                    "id":        f"DC{i+1:03d}",
                    "smiles":    smiles,
                    "name":      f"Candidate-{i+1:03d}",
                    "mechanism": "REINVENT 4 de novo generation",
                })
    return candidates


def generate_candidates() -> list:
    if not os.path.exists(config.REINVENT_PRIOR_PATH):
        raise FileNotFoundError(f"Prior model bulunamadı: {config.REINVENT_PRIOR_PATH}")

    with tempfile.TemporaryDirectory() as tmpdir:
        toml_path  = os.path.join(tmpdir, "reinvent_sample.toml")
        output_csv = os.path.join(tmpdir, "sampled.csv")

        _write_toml(toml_path, output_csv)
        _run_reinvent(toml_path)

        # output_file belirtilmişse oraya yazar, yoksa samples.csv varsayılan
        if not os.path.exists(output_csv):
            fallback = os.path.join(tmpdir, "samples.csv")
            if os.path.exists(fallback):
                output_csv = fallback
            else:
                raise FileNotFoundError("REINVENT 4 çıktı CSV dosyası bulunamadı.")

        candidates = _parse_output(output_csv)

    print(f"  ✓ REINVENT 4: {len(candidates)} SMILES üretildi.")
    return candidates