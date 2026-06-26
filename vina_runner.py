"""
Drug Discovery Platform v2 - AutoDock Vina Docking Pipeline
"""
import os
import subprocess
import tempfile
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
import meeko
import config

# Windows: vina.exe, Linux: vina
_base = os.path.dirname(os.path.abspath(__file__))
if os.name == "nt":
    VINA_EXE = os.path.join(_base, "vina.exe")
else:
    VINA_EXE = os.path.join(_base, "vina")


def _normalize_vina(kcal_mol):
    best    = config.VINA_BEST_SCORE
    worst   = config.VINA_WORST_SCORE
    span    = best - worst
    clamped = max(best, min(worst, kcal_mol))
    normalized = (clamped - worst) / span * 10.0
    return round(max(0.0, min(10.0, normalized)), 2)


def prepare_receptor(force=False):
    pdb_path   = config.RECEPTOR_PDB_PATH
    pdbqt_path = config.RECEPTOR_PDBQT_PATH

    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Receptor PDB not found: {pdb_path}")

    if os.path.exists(pdbqt_path) and not force:
        print(f"  Receptor PDBQT exists: {pdbqt_path}")
        return pdbqt_path

    print(f"  Preparing receptor: {pdb_path} -> {pdbqt_path}")

    # Clean with pdbfixer
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        fixer = PDBFixer(filename=pdb_path)
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        fixed_pdb = pdb_path.replace(".pdb", "_fixed.pdb")
        with open(fixed_pdb, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        print(f"  pdbfixer: fixed PDB saved to {fixed_pdb}")
    except ImportError:
        print("  pdbfixer/openmm not found, using raw PDB.")
        fixed_pdb = pdb_path

    # Convert to PDBQT (Windows=ADFR, Linux=obabel)
    if os.name == "nt":
        adfr_prepare = r"C:\Program Files (x86)\ADFRsuite-1.0\bin\prepare_receptor.bat"
        print(f"  Converting to PDBQT with ADFR prepare_receptor...")
        result = subprocess.run(
            [adfr_prepare, "-r", fixed_pdb, "-o", pdbqt_path,
             "-A", "checkhydrogens", "-U", "nphs_lps_waters_nonstdres"],
            capture_output=True, text=True
        )
        print(f"  prepare_receptor output: {result.stdout[-300:]}")
    else:
        print(f"  Converting to PDBQT with obabel...")
        result = subprocess.run(
            ["obabel", fixed_pdb, "-O", pdbqt_path,
             "-xr",  # receptor mode — correct Vina atom types
             "--partialcharge", "gasteiger"],
            capture_output=True, text=True
        )
        print(f"  obabel: {result.returncode} | {result.stderr[-300:]}")

    if not os.path.exists(pdbqt_path):
        raise RuntimeError(
            f"Receptor PDBQT could not be created.\n"
            f"Output: {result.stdout}\n{result.stderr[:300]}"
        )

    print(f"  Receptor ready: {pdbqt_path}")
    return pdbqt_path


def _smiles_to_pdbqt(smiles, out_path):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            result = AllChem.EmbedMolecule(mol)
            if result == -1:
                return False
        AllChem.MMFFOptimizeMolecule(mol)
        prep = meeko.MoleculePreparation()
        mol_setups = prep.prepare(mol)
        result = meeko.PDBQTWriterLegacy.write_string(mol_setups[0])
        if isinstance(result, tuple):
            pdbqt_string = result[0]
        else:
            pdbqt_string = result
        with open(out_path, "w") as f:
            f.write(pdbqt_string)
        return True
    except Exception as e:
        print(f"          PDBQT error: {e}")
        return False


def _parse_vina_score(output):
    import re
    for line in output.splitlines():
        m = re.match(r'^\s+1\s+([-\d.]+)', line)
        if m:
            return float(m.group(1))
    for line in output.splitlines():
        m = re.search(r'VINA RESULT:\s+([-\d.]+)', line)
        if m:
            return float(m.group(1))
    return None


def dock_single(smiles, receptor_pdbqt):
    if not os.path.exists(VINA_EXE):
        raise FileNotFoundError(f"Vina executable not found: {VINA_EXE}")

    with tempfile.TemporaryDirectory() as tmpdir:
        ligand_pdbqt = os.path.join(tmpdir, "ligand.pdbqt")
        output_pdbqt = os.path.join(tmpdir, "docked.pdbqt")

        if not _smiles_to_pdbqt(smiles, ligand_pdbqt):
            return {"vina_score_kcal": 0.0, "binding_score": 0.0, "success": False}

        try:
            cmd = [
                VINA_EXE,
                "--receptor",       receptor_pdbqt,
                "--ligand",         ligand_pdbqt,
                "--out",            output_pdbqt,
                "--center_x",       str(config.BINDING_SITE["center_x"]),
                "--center_y",       str(config.BINDING_SITE["center_y"]),
                "--center_z",       str(config.BINDING_SITE["center_z"]),
                "--size_x",         str(config.BINDING_SITE["size_x"]),
                "--size_y",         str(config.BINDING_SITE["size_y"]),
                "--size_z",         str(config.BINDING_SITE["size_z"]),
                "--exhaustiveness", str(config.VINA_EXHAUSTIVENESS),
                "--num_modes",      str(config.VINA_NUM_POSES),
                "--verbosity",      "1",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                print(f"          Vina error: {result.stderr[:200]}")
                return {"vina_score_kcal": 0.0, "binding_score": 0.0, "success": False}
            best_score = _parse_vina_score(result.stdout + result.stderr)
            if best_score is None:
                print(f"          Could not parse Vina score: {(result.stdout+result.stderr)[:200]}")
                return {"vina_score_kcal": 0.0, "binding_score": 0.0, "success": False}
            return {
                "vina_score_kcal": round(best_score, 2),
                "binding_score":   _normalize_vina(best_score),
                "success":         True,
            }
        except Exception:
            return {"vina_score_kcal": 0.0, "binding_score": 0.0, "success": False}


def dock_all(candidates, receptor_pdbqt):
    total   = len(candidates)
    results = []
    for i, c in enumerate(candidates, 1):
        print(f"  [{i:>3}/{total}] Docking: {c['id']} ({c['smiles'][:40]}...)")
        dock_result = dock_single(c["smiles"], receptor_pdbqt)
        results.append({**c, **dock_result})
        if dock_result["success"]:
            print(f"          OK {dock_result['vina_score_kcal']:.1f} kcal/mol -> {dock_result['binding_score']:.1f}/10")
        else:
            print(f"          FAILED")
    success_count = sum(1 for r in results if r["success"])
    print(f"  Docking complete: {success_count}/{total} successful.")
    return results
