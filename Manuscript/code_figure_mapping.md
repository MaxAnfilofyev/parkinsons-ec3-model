# Code–Manuscript Mapping for `pd_cleaned.md`

This document maps each Python file in the `code/` directory to the figures and visuals used in `Manuscript/pd_cleaned.md`. It also flags where visuals or supporting code are missing or misaligned with the current text.

**Status legend**

- `OK` – Code exists and matches the manuscript figure/text.
- `CONCEPTUAL` – Figure is a conceptual/mermaid diagram defined directly in the markdown (no Python code).
- `PARTIAL` – Code exists, but semantics or labelling do not fully match the text.
- `MISSING CODE` – Manuscript describes a visual or analysis for which there is no corresponding code in `code/`.
- `LEGACY / UNUSED` – Code produces figures that are not referenced in `pd_cleaned.md`.

---

## 1. Scripts in `code/` – Outputs and Usage

| Script | Main purpose | Key outputs (PNG/CSV) | Used in `pd_cleaned.md`? | Notes |
| ------ | ------------ | --------------------- | ------------------------ | ----- |
| `run_ec3_all.py` | Orchestrates generation of EC3 minimal model figures across three stages (bifurcation, phase planes, time courses) with manuscript-style plotting. | `Stage_2_Bifurcation/<run>/ec3_bifurcation_E_vs_A.png`; `Stage_3_PhasePlanes/<run>/phase_plane_VTA_A_0.40.png`; `Stage_3_PhasePlanes/<run>/phase_plane_SNc_A_1.00.png`; `Stage_4_TimeCourses/<run>/timecourses_VTA_vs_SNc_baseline.png`; `Stage_4_TimeCourses/<run>/timecourses_SNc_vs_VTA_perturbation.png` plus associated CSVs and README files. | **Yes (main text)** | This is the primary driver for all EC3 figures referenced in Sections 2.2–5 of `pd_cleaned.md` via the `Stage_*/*/latest_run/*.png` paths. |
| `ec3_bistable.py` | Standalone EC3 script: scans equilibria vs axonal load `A`, classifies equilibria, and generates a bifurcation plot, a single phase plane at a representative bistable `A`, and matching time courses. | `Stage_2_Bifurcation/runs/<timestamp>/ec3_bifurcation_E_vs_A.png`; `Stage_2_Bifurcation/runs/<timestamp>/ec3_phase_plane_A_<A>.png`; `Stage_2_Bifurcation/runs/<timestamp>/ec3_timecourses_A_<A>.png`; plus `ec3_equilibria_vs_A.csv`, `ec3_bifurcation_summary.txt`. | **Not directly** | Conceptually identical to the Stage 2/3/4 logic used in `run_ec3_all.py`, but `pd_cleaned.md` links to the `latest_run` images produced by `run_ec3_all.py`. Can be considered an alternative/legacy figure generator for EC3. |
| `ec2_stability_scanner.py` | Numerical EC2 scan over `(gamma, beta)` using full ODE integration to detect bistability (via steady states from high/low initial conditions). | `ec2_bistability_scan.csv`; `ec2_deltaE_heatmap.png`; `ec2_bistable_mask.png`; `ec2_example_timecourses.png`. | **No** (for `pd_cleaned.md`) | Supports EC2 robustness-style analyses but EC2 is not explicitly described or cited in the current `pd_cleaned.md`. Potentially useful for future or separate supplementary material. |
| `ec_analytic_bistable.py` | Analytic EC2 fixed-point analysis using cubic roots in `M`; performs fast `(gamma, beta)` sweeps, classifies fixed points and stability, and plots heatmaps. | `ec2_analytic_bistability_maps.npz`; `ec2_analytic_deltaE_heatmap.png`; `ec2_analytic_nstable_heatmap.png`. | **No** (for `pd_cleaned.md`) | The Supplementary text (Section S3) currently describes parameter sweeps over the *final EC3 model* (varying `k_M`, `L_1`, `C`), not EC2 parameter sweeps over `gamma`, `beta`. This script does **not** directly generate the visuals described in S3 of `pd_cleaned.md`. |
| `supp.py` | EC2-based supplementary figure generator; integrates a 2D EC2 system and produces EC2 time courses, bifurcation vs `A`, and a phase plane. | `supplement_ec2_S1_timecourses.png`; `supplement_ec2_S2_bifurcation_A.png`; `supplement_ec2_S3_phase_plane.png`. | **PARTIAL / MISMATCH** | File header states: “Used ONLY for Supplementary Figures S1–S3.” However, the Supplementary S1–S3 text in `pd_cleaned.md` now describes *earlier monostable/feedback-only architectures* and a requirement for nonlinear energy amplification, **not** EC2 per se. The code generates EC2 figures that are not cleanly aligned with the current Supplementary text. |
| `pd_model_figures.py` | 3-variable model (V2) including energy, mitochondrial capacity, and α-syn (`E, M, S`); generates time courses, A×C heatmaps, 1D sweeps, and supplementary trajectories/sensitivity plots. | `figure2_time_courses_v2.png`; `figure3_heatmap_v2_E.png`; `figure3_heatmap_v2_S.png`; `figure4_tipping_point_v2.png`; `supplement_trajectories_v2.png`; `supplement_sensitivity_v2.png`. | **LEGACY / UNUSED** | This earlier 3D model and its figures are not referenced anywhere in `pd_cleaned.md` (which now focuses on the 2-variable EC3 model). These figures could be used in a separate manuscript or an extended supplement but are **not** required for the current narrative. |
| `pd_model_energy_crisis.py` | 3-variable “energetic-crisis” model (EC, V3) similar in structure to V2 but with more detailed α-syn dynamics; generates analogous figure set plus bifurcation and phase-plane plots. | `figure2_time_courses_ec.png`; `figure3_heatmap_ec_E.png`; `figure3_heatmap_ec_S.png`; `figure4_tipping_point_ec.png`; `figure5_bifurcation_A_ec.png`; `figure6_phase_plane_bistable_ec.png`; `supplement_trajectories_ec.png`; `supplement_sensitivity_ec.png`. | **LEGACY / UNUSED** | Like `pd_model_figures.py`, this extended EC model is not cited in `pd_cleaned.md`. The cleaned manuscript treats EC3 as the core minimal model and does not discuss the 3-variable EC systems. |

---

## 2. Manuscript Figures and Visuals → Generating Code

### 2.1 Main-text figures

**Figure 1 – Minimal energetic model of dopaminergic neurons (Section 2.1)**

- **Type:** Mermaid diagram embedded directly in `pd_cleaned.md`.
- **Code:** No Python code; defined as a fenced ```mermaid``` block in the manuscript.
- **Status:** `CONCEPTUAL` (not generated via `code/`).

**Figure 2A – Conceptual energy landscape under load (Section 2.2)**

- **Type:** Mermaid diagram embedded directly in `pd_cleaned.md`.
- **Code:** No Python code; defined as a second ```mermaid``` block.
- **Status:** `CONCEPTUAL`.

**Figure 2B – Phase portrait at SNc-like load (A = 1.00) (Section 2.2)**

- **Path in manuscript:**
  - `![](/Stage_3_PhasePlanes/latest_run/phase_plane_SNc_A_1.00.png)`
- **Generating code:**
  - `run_ec3_all.py` → Stage 3:
    - `make_phase_plane_plot(p_base, A_SNc, os.path.join(run_dir_stage3, "phase_plane_SNc_A_1.00.png"), ...)`
- **Status:** `OK`.
- **Notes:** This is the concrete dynamical realization of the conceptual landscape in Figure 2A. Methods Section 7.7 describes phase-plane computation consistent with this script (up to numerical grid resolution; see Section 3.3 below).

**Figure 3 – Saddle-node bifurcation of energetic reserve vs axonal load (Section 3)**

- **Path in manuscript:**
  - `![](/Stage_2_Bifurcation/latest_run/ec3_bifurcation_E_vs_A.png)`
- **Generating code:**
  - `run_ec3_all.py` → Stage 2:
    - `scan_equilibria_vs_A(p_base, A_min=0.2, A_max=1.4, nA=61, ...)`
    - `plot_bifurcation_E_vs_A(data_rows, os.path.join(run_dir_stage2, "ec3_bifurcation_E_vs_A.png"))`
- **Status:** `OK`.
- **Notes:** The A-range and use of 61 points match Methods Section 7.6. Stability classification via the Jacobian is implemented in `classify_equilibrium`.

**Figure 4A – Phase plane at low axonal load (VTA-like, A = 0.40) (Section 4)**

- **Path in manuscript:**
  - First column of the Figure 4 table: `![](/Stage_3_PhasePlanes/latest_run/phase_plane_VTA_A_0.40.png)`
- **Generating code:**
  - `run_ec3_all.py` → Stage 3:
    - `make_phase_plane_plot(p_base, A_VTA, os.path.join(run_dir_stage3, "phase_plane_VTA_A_0.40.png"), ...)`
- **Status:** `OK`.

**Figure 4B – Phase plane at high axonal load (SNc-like, A = 1.00) (Section 4)**

- **Path in manuscript:**
  - Second column of the Figure 4 table: `![](/Stage_3_PhasePlanes/latest_run/phase_plane_SNc_A_1.00.png)`
- **Generating code:**
  - Same as Figure 2B: `run_ec3_all.py` Stage 3, with `A_SNc = 1.00`.
- **Status:** `OK`.
- **Notes:** The same underlying image is used for both Figure 2B and 4B, but with different narrative emphasis.

**Figure 5A – Baseline trajectories of E(t) under VTA- vs SNc-like load (Section 5)**

- **Path in manuscript:**
  - Table column image: `![](/Stage_4_TimeCourses/latest_run/timecourses_VTA_vs_SNc_baseline.png)`
- **Generating code:**
  - `run_ec3_all.py` → Stage 4:
    - `simulate_timecourse(p_base, A_VTA, y0_common)` and `simulate_timecourse(p_base, A_SNc, y0_common)`
    - Plotted together and saved as `timecourses_VTA_vs_SNc_baseline.png`.
- **Status:** `OK`.

**Figure 5B – Response to a brief energetic insult (collapse vs recovery) (Section 5)**

- **Path in manuscript:**
  - Table column image: `![](/Stage_4_TimeCourses/latest_run/timecourses_SNc_vs_VTA_perturbation.png)`
- **Generating code:**
  - `run_ec3_all.py` → Stage 4:
    - `simulate_with_perturbation(p_base, A_SNc, y_snc_ss, t_pre=50, E_after=0.3, ...)`
    - `simulate_with_perturbation(p_base, A_VTA, y_vta_ss, t_pre=50, E_after=0.3, ...)`
    - Plotted together and saved as `timecourses_SNc_vs_VTA_perturbation.png`.
- **Status:** `OK`.

### 2.2 Supplementary figures described in `pd_cleaned.md`

The Supplementary Results section S1–S5 describes several figures (S1–S7) conceptually, but the markdown does **not** include explicit image paths. Below, we list each described figure and whether there is code in `code/` that directly generates it.

> **Important:** Several supplementary descriptions refer to **earlier, non-bistable model variants**, while the implemented EC2 scripts (`supp.py`, `ec2_stability_scanner.py`, `ec_analytic_bistable.py`) encode different dynamics or parameter-sweep structures. These misalignments are flagged as `MISMATCH` or `MISSING CODE` below.

#### Supplementary Figure S1 – Nullclines & vector field for an earlier, monostable model (Section S1.1)

- **Text description:**
  - *“In a model with linear mitochondrial support and load effects but no nonlinear energy amplification, the energy and mitochondrial nullclines intersect only once … global convergence to a single equilibrium … no bistability.”*
- **Expected content:**
  - Nullclines and vector field for an earlier, strictly monostable architecture (no `k_2 E^2(1-E)` term).
- **Available code in `code/`:**
  - None of the current files (`run_ec3_all.py`, `ec3_bistable.py`, `supp.py`, `ec2_stability_scanner.py`, `ec_analytic_bistable.py`) implement this specific earlier “linear support + linear load + no nonlinear amplification” architecture.
- **Status:** `MISSING CODE` for the exact S1 figure.
- **Related code:**
  - `supp.py` produces EC2 timecourses and a phase plane, but with a different EC2 formulation (nonlinear M-dependence, additional couplings). It is **not** the monostable architecture described in S1.1.

#### Supplementary Figure S2 – Bifurcation scan for a feedback-only model (Section S1.2)

- **Text description:**
  - *“Even with feedback from energy to mitochondrial damage, the system exhibits only a single equilibrium … no saddle-node bifurcation.”*
- **Expected content:**
  - Bifurcation diagram vs `A` for a feedback-only model (without nonlinear energy amplification), showing only a single equilibrium across the scan.
- **Available code in `code/`:**
  - No script explicitly implements this “feedback-only” architecture separate from the final EC3 model.
- **Status:** `MISSING CODE`.

#### Supplementary Figure S3 – Effect of nonlinear energy amplification on nullcline geometry (Section S2)

- **Text description:**
  - Comparison of energy nullclines with and without the `k_2 E^2(1-E)` term, demonstrating how the fold and bistability arise.
- **Expected content:**
  - Plots of energy nullcline(s) and their intersections with the mitochondrial nullcline, with and without nonlinear amplification.
- **Available code in `code/`:**
  - EC3 nullclines and equilibria can be inferred from `ec3_bistable.py` and `run_ec3_all.py` (which compute nullclines numerically), but there is **no dedicated script** that toggles the nonlinear term for a direct S3 comparison.
- **Status:** `MISSING CODE` (for the explicit “with vs without nonlinear term” figure).

#### Supplementary Figures S4–S6 – Parameter sweeps confirming robustness (Section S3)

- **Text description:**
  - Parameter sweeps over **EC3** parameters:
    - S3.1: varying `k_M` (mitochondrial turnover).
    - S3.2: varying `L_1` (load-dependent consumption).
    - S3.3: varying `C` (Ca²⁺ load).
  - Goal: show that the saddle-node structure persists over broad parameter ranges.
- **Expected content:**
  - Heatmaps or bifurcation diagrams showing presence/absence of bistability as parameters vary.
- **Available code in `code/`:**
  - `ec2_stability_scanner.py` and `ec_analytic_bistable.py` perform parameter sweeps over EC2 parameters (`gamma`, `beta`), not over EC3 parameters (`k_M`, `L_1`, `C`).
  - There is no script that sweeps EC3 parameters in the manner described in S3.
- **Status:** `MISMATCH` / `MISSING CODE`.
  - EC2 analytic/numerical sweeps exist but are **not** the EC3 robustness sweeps described in the Supplementary text.

#### Supplementary Figure S7 – Full continuation showing right fold beyond biological range (Section S4)

- **Text description:**
  - *“The complete bifurcation diagram reveals a right-hand fold at axonal loads larger than those observed in dopaminergic neurons …”*
- **Expected content:**
  - An extended bifurcation diagram of EC3 vs `A`, including a right-hand fold at `A > 1.4`, beyond the physiological range.
- **Available code in `code/`:**
  - EC3 bifurcation scan vs `A` is implemented in `ec3_bistable.py` and `run_ec3_all.py` (via `scan_equilibria_vs_A`), but currently uses `A_min=0.2`, `A_max=1.4`.
- **Status:** `PARTIAL`.
  - **What exists:** All infrastructure to perform continuation in `A` and classify equilibria.
  - **What’s missing:** A variant of the scan with `A_max` extended beyond 1.4 and a dedicated plotting routine for the full continuation as described in S4.
- **Suggested update:**
  - Add a small wrapper or parameterized option in `run_ec3_all.py`/`ec3_bistable.py` to extend the `A` range and output a `supplement_S7_full_bifurcation.png` (or similar).

#### EC2-based supplementary figures from `supp.py`

Although `supp.py` is not currently referenced in the markdown, it generates:

- `supplement_ec2_S1_timecourses.png` – EC2 timecourses at a bistable `A`.
- `supplement_ec2_S2_bifurcation_A.png` – EC2 steady-state `E*` and `M*` vs `A` for high/low initial conditions.
- `supplement_ec2_S3_phase_plane.png` – EC2 phase plane at a bistable `A`.

**Status relative to `pd_cleaned.md`:** `PARTIAL / MISMATCH`.

- If the intention is that Supplementary S1–S3 should present EC2 (as an intermediate 2D energetic-catastrophe model distinct from EC3), then the Supplementary text should be adjusted to match the EC2 equations and interpretation currently implemented in `supp.py`.
- If instead S1–S3 are meant to present *earlier, monostable and feedback-only architectures* of EC3 (as currently written), then either:
  - The EC2-specific references in `supp.py` should be relabelled (e.g., as separate EC2 supplement), or
  - New scripts matching those earlier architectures should be added.

---

## 3. Section-by-Section Coverage Check

### 3.1 Abstract & Introduction (Sections 1)

- **Visuals required:** None explicitly; narrative is supported by later figures.
- **Code dependency:** None beyond what is already used for Figures 2–5.
- **Status:** `OK`.

### 3.2 Minimal Energetic Model (Section 2)

- **Visuals:**
  - Figure 1 – Conceptual schematic (`CONCEPTUAL` mermaid; no Python code).
  - Figure 2A – Conceptual energy landscape (`CONCEPTUAL` mermaid).
  - Figure 2B – Phase portrait at SNc-like load (A = 1.00) generated by `run_ec3_all.py`.
- **Code coverage:**
  - `run_ec3_all.py` Stage 3 fully supports Figure 2B.
- **Status:** `OK`.

### 3.3 Bifurcation Analysis (Section 3)

- **Main visual:** Figure 3 (EC3 bifurcation diagram vs `A`), generated by `run_ec3_all.py` Stage 2 (`ec3_bifurcation_E_vs_A.png`).
- **Additional visuals referenced:**
  - *“Supplementary Figures S1–S3”* for extended parameter sweeps / numerical continuation.
- **Code coverage:**
  - Main Figure 3: `OK` (via `run_ec3_all.py`).
  - Supplementary S1–S3 as described in the text: **not fully implemented** (see Section 2.2 above).
- **Status:**
  - Main figure: `OK`.
  - Supplementary references: `MISSING CODE` / `MISMATCH`.

### 3.4 SNc vs VTA comparison (Section 4)

- **Visuals:** Figures 4A and 4B – Phase planes at `A = 0.40` and `A = 1.00`.
- **Code coverage:**
  - Both images produced by `run_ec3_all.py` Stage 3 via `make_phase_plane_plot`.
- **Status:** `OK`.

### 3.5 Perturbation experiments (Section 5)

- **Visuals:** Figures 5A and 5B – baseline trajectories and response to a brief insult.
- **Code coverage:**
  - `run_ec3_all.py` Stage 4 generates both figures using `simulate_timecourse` and `simulate_with_perturbation`.
- **Status:** `OK`.

### 3.6 Discussion (Section 6)

- **Visuals:** None.
- **Status:** `OK` (relies conceptually on Figures 2–5, which are covered).

### 3.7 Methods (Section 7)

Methods outline numerical settings and procedures for simulations, bifurcation scans, and phase-plan computation.

- **7.1–7.3 (Model equations, parameters, initial conditions):**
  - Qualitative correspondence with EC3 equations in `run_ec3_all.py` and `ec3_bistable.py`.
  - Note: The parameter values listed (e.g., `k_1`, `k_2`, `L_0`, `L_1`, `k_M`, `β`) are schematic and do **not** match the calibrated EC3 values hard-coded in `PDParamsEC3`. This is expected given the manuscript’s conceptual presentation.
- **7.4 Numerical integration:**
  - Methods specify RK45 with `rtol=1e-8`, `atol=1e-10`, `max_step=0.1`.
  - Actual code in `run_ec3_all.py` and `ec3_bistable.py` calls `solve_ivp` with default tolerances (i.e., they do not explicitly set `rtol`, `atol`, `max_step`).
  - **Status:** `PARTIAL` – Figures are qualitatively correct, but if strict reproducibility is required, code should be updated to match the stated tolerances.
- **7.5 Steady states and stability:**
  - Implemented in `ec3_bistable.py` and `run_ec3_all.py` via `find_equilibria` and `classify_equilibrium` (using `scipy.optimize.root`).
  - Methods refer to `fsolve`, which is conceptually identical to `root(..., method='hybr')`.
  - **Status:** `OK`.
- **7.6 Bifurcation scan procedure:**
  - Exactly implemented in `scan_equilibria_vs_A` (0.2–1.4, 61 points).
  - **Status:** `OK`.
- **7.7 Phase plane and nullcline computation:**
  - Methods state a 200×200 grid; `run_ec3_all.py`/`ec3_bistable.py` currently use a 25×25 grid for the vector field.
  - **Status:** `PARTIAL` – Qualitative geometry is correctly captured; grid resolution is lower than described. This can be updated if needed.

### 3.8 Supplementary Results (S1–S5)

- **Visuals and code:**
  - S1–S2 describe earlier monostable and feedback-only models; no matching code in `code/`.
  - S2 describes the role of nonlinear energy amplification; no dedicated script to toggle this term.
  - S3 describes EC3 robustness sweeps over `k_M`, `L_1`, `C`; only EC2 sweeps exist in `ec2_stability_scanner.py`/`ec_analytic_bistable.py`.
  - S4 describes a full continuation including a right-hand fold beyond biological range; `scan_equilibria_vs_A` can be extended to support this, but is not currently configured to do so.
  - S5 is textual summary; no additional visuals.
- **Status:**
  - Conceptual narrative: `OK`.
  - Concrete figures S1–S7 as described: `PARTIAL` / `MISSING CODE` as detailed in Section 2.2.

---

## 4. Identified Gaps and Recommended Updates

This section summarizes where visuals or generating code are missing or inconsistent with `pd_cleaned.md`.

### 4.1 Main-text figures

- All main-text figures (Figures 1–5) are either conceptual (mermaid) or backed by `run_ec3_all.py` and the `Stage_*/*/latest_run/*.png` outputs.
- **No missing visuals** for Sections 2–5.

### 4.2 Supplementary figures (S1–S7)

1. **S1 & S2 – Earlier monostable and feedback-only variants**
   - **Issue:** No scripts in `code/` implement these exact architectures.
   - **Options:**
     - (a) Add scripts implementing these earlier models, generating the S1 and S2 figures described in the Supplementary text; or
     - (b) Simplify/remove these references if you do not intend to include the figures.

2. **S3 – Effect of nonlinear energy amplification**
   - **Issue:** No dedicated script toggles the `k_2 E^2(1-E)` term to produce explicit comparisons of nullcline geometry.
   - **Options:**
     - Implement a small EC3 helper (could be integrated into `ec3_bistable.py`) that plots energy nullclines with and without the nonlinear term.

3. **S4–S6 – EC3 robustness sweeps**
   - **Issue:** Text describes sweeps over EC3 parameters (`k_M`, `L_1`, `C`), but only EC2 (`gamma`, `beta`) sweeps are implemented.
   - **Options:**
     - Add an EC3 sweep driver (analogous to `ec_analytic_bistable.py`) that classifies equilibria as a function of `k_M`, `L_1`, `C`, and saves heatmaps/bifurcation-like plots for S4–S6.

4. **S7 – Right-hand fold beyond biological range**
   - **Issue:**
     - `scan_equilibria_vs_A` infrastructure exists, but the A-range is currently limited to [0.2, 1.4].
   - **Options:**
     - Expose `A_min`/`A_max` as arguments in a supplementary wrapper to extend beyond 1.4 and generate the full continuation figure for S7.

5. **EC2-based `supp.py` output vs Supplementary text**
   - **Issue:** `supp.py` generates EC2 figures labelled internally as S1–S3, but the current Supplementary text uses S1–S3 to discuss different (earlier) model variants.
   - **Options:**
     - (a) Update the Supplementary text to explicitly present EC2 as the focus of S1–S3, with equations matching `PDModelParamsEC2` and `pd_ode_ec2`; or
     - (b) Relabel the EC2 figures and reserve S1–S3 for the conceptual/earlier EC3 variants, adding matching code for those.

### 4.3 Methods vs implementation

- **Numerical tolerances and grid resolution:**
  - If strict numerical reproducibility is a priority, update `run_ec3_all.py` and `ec3_bistable.py` to:
    - Pass `rtol`, `atol`, and `max_step` as specified in Section 7.4.
    - Increase vector-field grids to 200×200 in phase-plane plots, as described in Section 7.7.
- These changes are not required for the qualitative narrative but would align the code more precisely with the Methods description.

### 4.4 Legacy 3-variable models

- `pd_model_figures.py` and `pd_model_energy_crisis.py` currently generate a rich set of 3D model figures that are **not referenced** in `pd_cleaned.md`.
- Consider either:
  - (a) Explicitly marking these as legacy/archival in the repository (e.g., moving to `archive/` or noting in a README); or
  - (b) Adding a brief note in the Supplementary section explaining that extended EC models (with α-syn and proteostasis) are available in code but beyond the scope of the present manuscript.

---

## 5. How to Regenerate All Main Figures

For convenience, the main figures supporting `pd_cleaned.md` can be regenerated using:

```bash
python code/run_ec3_all.py
```

This will create timestamped run directories under:

- `Stage_2_Bifurcation/runs/<timestamp>/`
- `Stage_3_PhasePlanes/runs/<timestamp>/`
- `Stage_4_TimeCourses/runs/<timestamp>/`

and update the `Stage_*/*/latest_run` symlinks used by the image paths in `pd_cleaned.md`.
