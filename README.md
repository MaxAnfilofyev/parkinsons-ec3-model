# Parkinson’s Energetic Collapse (EC3) Model

This repository contains code and manuscript materials for a computational model of energetic collapse in Parkinson’s disease (the "EC3" model). It includes:

- Python code to run the main simulation and analysis stages
- Saved outputs from recent runs
- A manuscript written in Markdown
- Utilities to build the manuscript PDF

For a quick status view of the project, see **`PROJECT_DASHBOARD.md`**.

---

## Repository Layout

- **`code/`** – Core Python scripts and utilities
  - `run_ec3_all.py` – Main EC3 model runner and helpers
  - `s1_s2_earlier_variants.py` – Stage 1–2 earlier model variants
  - `s3_nonlinear_amplification.py` – Stage 3 nonlinear amplification analyses
  - `s4_s6_robustness_sweeps.py` – Robustness / parameter sweeps
  - `s7_full_continuation.py` – Full continuation / bifurcation-related analyses
  - `stage_utils.py` – Helper utilities for stage directories and README generation
  - `archive/` – Older/experimental model scripts
  - `md_to_pdf.sh` – Convenience script to build the manuscript PDF from Markdown

- **Stage output directories** (each with `runs/` subfolders and a `latest_run/` symlink-like directory):
  - `Stage_1_EC3_Model/`
  - `Stage_2_Bifurcation/`
  - `Stage_3_PhasePlanes/`
  - `Stage_4_TimeCourses/`
  - `Stage_S1_S2_EarlierVariants/`
  - `Stage_S3_NonlinearAmplification/`
  - `Stage_S4_S6_RobustnessSweeps/`
  - `Stage_S7_FullContinuation/`

- **`Manuscript/`** – Manuscript source and supporting files
  - `pd_manuscript.md` – Main manuscript in Markdown
  - `pd_cleaned.pdf` – A compiled PDF version of the manuscript
  - `code_figure_mapping.md` – Mapping between manuscript figures and the code/stages that generate them
  - `significance.md` – Significance statement text

- **`outputs/`** – (Optional) additional figures or derived outputs.

- **Top-level files**
  - `requirements.txt` – Python dependencies for the simulation code
  - `PROJECT_DASHBOARD.md` – High-level project status and open questions

---

## Python Environment

The simulation code is written for **Python ≥ 3.7**.

Install dependencies (preferably in a virtual environment):

```bash
pip install -r requirements.txt
```

`requirements.txt` includes:

- `numpy` – numerical arrays and math
- `scipy` – ODE integration and root-finding
- `matplotlib` – figure/plot generation

---

## How to Run the Code

From the repository root, you can invoke the stage scripts directly. For example:

```bash
# Example: run the main EC3 model / equilibrium scans
python code/run_ec3_all.py

# Stage-specific analysis scripts
python code/s1_s2_earlier_variants.py
python code/s3_nonlinear_amplification.py
python code/s4_s6_robustness_sweeps.py
python code/s7_full_continuation.py
```

Each script writes results into its corresponding `Stage_*` directory, usually under a timestamped folder in `runs/`, and updates the `latest_run/` pointer.

Refer to `PROJECT_DASHBOARD.md` for which stages correspond to which scientific questions and their current completion status.

---

## Manuscript Location and Build Instructions

The manuscript lives under **`Manuscript/`**:

- Main text (Markdown): `Manuscript/pd_manuscript.md`
- Existing compiled PDF: `Manuscript/pd_cleaned.pdf`
- Code-to-figure mapping: `Manuscript/code_figure_mapping.md`
- Significance statement: `Manuscript/significance.md`

To rebuild the manuscript PDF from Markdown, you can use the `md-to-pdf` CLI (Node-based tool) and the helper script in `code/`.

1. Ensure `md-to-pdf` is installed, e.g.:
   ```bash
   npm install -g md-to-pdf
   ```

2. From the repository root, run:
   ```bash
   bash code/md_to_pdf.sh
   ```

This uses the custom config at `pdf-utils/md-to-pdf-config.js` and converts `Manuscript/pd_manuscript.md` to a PDF (by default, saved alongside the Markdown or as configured in that file).

---

## Mapping Figures to Code

To see which scripts and stages generate each figure in the manuscript, consult:

- **`Manuscript/code_figure_mapping.md`**

This file provides the bridge between manuscript figures and the corresponding Python code and stage runs.

---

## Getting Started

1. Create and activate a Python virtual environment.
2. Install dependencies with `pip install -r requirements.txt`.
3. Run one or more stage scripts in `code/` to regenerate the key analyses and figures.
4. Optionally, rebuild the manuscript PDF using `bash code/md_to_pdf.sh`.
5. Use `Manuscript/code_figure_mapping.md` as a guide to connect figures in the paper with the underlying code.
