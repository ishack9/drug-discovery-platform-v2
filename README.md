# Drug Discovery Platform v2

An automated AI-powered drug discovery pipeline combining generative molecular design, molecular docking, and ADMET prediction.

## Tools

| Tool | Purpose |
|------|---------|
| [REINVENT 4](https://github.com/MolecularAI/REINVENT4) | Generative molecular design (AstraZeneca) |
| [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) | Molecular docking |
| [admet-ai](https://github.com/swansonk14/admet_ai) | ADMET prediction (50+ properties) |
| [RDKit](https://www.rdkit.org/) | Cheminformatics, SA Score, Tanimoto similarity |

## Workflow

```
REINVENT 4 (60 SMILES)
    ↓
AutoDock Vina (molecular docking)
    ↓
admet-ai (ADMET evaluation)
    ↓
Weighted scoring + PDF/CSV report
```

## Scoring Formula

```
Composite = Binding×0.40 + ADMET×0.35 + Manufacturability×0.15 + Novelty×0.10
```

## Features

- **RL Mode**: REINVENT 4 reinforcement learning guided toward EGFR inhibitor scaffolds
- **ADMET profiling**: hERG cardiotoxicity, DILI liver injury, AMES mutagenicity, oral bioavailability, 50+ properties
- **Validation**: Known EGFR inhibitors (erlotinib, gefitinib, osimertinib, afatinib) vs negative control (imatinib)
- **Traffic light system**: Green/Yellow/Red for each ADMET property
- **PDF reports**: Full candidate analysis with ADMET profiles
- **Gradio UI**: Web-based interface with real-time progress

## Installation (Windows)

### Requirements
- Anaconda
- AutoDock Vina 1.2.7 ([download](https://github.com/ccsb-scripps/AutoDock-Vina/releases))
- ADFRsuite ([download](https://ccsb.scripps.edu/adfr/))
- REINVENT 4 prior model ([Zenodo](https://zenodo.org/records/15641297))

### Setup

```bash
# Install dependencies
conda install -c conda-forge rdkit openbabel pdbfixer openmm gemmi
pip install meeko admet-ai reportlab gradio

# Install REINVENT 4
pip install git+https://github.com/MolecularAI/REINVENT4.git

# Run
python main.py EGFR

# Or launch GUI
python app.py
```

## Configuration

Edit `config.py`:

```python
RECEPTOR_PDB_PATH = "receptor/protein.pdb"
BINDING_SITE = {
    "center_x": 26.5,
    "center_y": 11.5,
    "center_z": -1.0,
    "size_x": 25.0,
    "size_y": 25.0,
    "size_z": 25.0,
}
USE_RL_MODE = True
```

## Demo

Live demo available on [Hugging Face Spaces](https://huggingface.co/spaces/ishackT9/drug-discovery-v2)

## Results

Example screening against EGFR (PDB: 4HJO):
- 60 candidates generated and docked
- Average composite score: 7.18/10
- 10 candidates selected for laboratory evaluation
- Best candidate: -11.5 kcal/mol binding affinity

## References

- Loeffler et al., *J. Cheminform.* 16, 20 (2024) — REINVENT 4
- Eberhardt et al., *J. Chem. Inf. Model.* 61, 3891 (2021) — AutoDock Vina
- Swanson et al., *J. Chem. Inf. Model.* 64, 2 (2024) — admet-ai

## License

MIT License