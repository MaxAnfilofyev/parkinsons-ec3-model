import os
import datetime
from typing import Tuple

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from stage_utils import make_stage_run_dir, write_readme, STAGE_NAMES
from run_ec3_all import PDParamsEC3, scan_equilibria_vs_A

# ============================================================
# GLOBAL PLOTTING STYLE (MATCHES run_ec3_all.py)
# ============================================================

mpl.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

COLOR_WINDOW = "#4b0082"  # dark purple


def new_panel(figsize: Tuple[float, float] = (3.3, 2.2), panel_label: str | None = None):
    """Create a new figure/axes pair with optional panel label.

    Mirrors the helper from `run_ec3_all.py` to keep sizing consistent.
    """
    fig, ax = plt.subplots(figsize=figsize)
    if panel_label is not None:
        ax.text(0.02, 0.95, panel_label, transform=ax.transAxes,
                fontweight="bold", va="top")
    return fig, ax


# ============================================================
# CORE HELPER: COMPUTE BISTABLE WINDOW VS A FOR VARYING PARAMETERS
# ============================================================


def compute_bistable_window_vs_param(
    param_name: str,
    param_values: np.ndarray,
    p_base: PDParamsEC3,
    A_min: float = 0.2,
    A_max: float = 1.4,
    nA: int = 61,
    ngrid: int = 7,
):
    """For each value of a chosen parameter, compute the A-range with 3 equilibria.

    Returns dict with:
      - "param_values": array of parameter values
      - "A_three_min": array of lower bounds (NaN if no bistable window)
      - "A_three_max": array of upper bounds (NaN if no bistable window)
    """
    A_three_min = np.full_like(param_values, np.nan, dtype=float)
    A_three_max = np.full_like(param_values, np.nan, dtype=float)

    for idx, val in enumerate(param_values):
        # Create a modified parameter set with the selected parameter changed.
        kwargs = p_base.__dict__.copy()
        kwargs[param_name] = float(val)
        p_mod = PDParamsEC3(**kwargs)

        _, summary = scan_equilibria_vs_A(
            p_mod,
            A_min=A_min,
            A_max=A_max,
            nA=nA,
            ngrid=ngrid,
        )

        if summary["A_three_min"] is not None:
            A_three_min[idx] = summary["A_three_min"]
            A_three_max[idx] = summary["A_three_max"]

    return {
        "param_values": param_values,
        "A_three_min": A_three_min,
        "A_three_max": A_three_max,
    }


def plot_window_vs_param(result: dict, param_label: str, filename: str, panel_label: str | None = None):
    """Plot bistable A-window (A_three_min / A_three_max) versus a parameter.

    Regions where the window exists are indicated by solid curves and optional
    vertical shading.
    """
    x = result["param_values"]
    A_min = result["A_three_min"]
    A_max = result["A_three_max"]

    fig, ax = new_panel(figsize=(3.3, 2.2), panel_label=panel_label)

    # Plot only points where a window exists
    mask = ~np.isnan(A_min)
    ax.plot(x[mask], A_min[mask], color=COLOR_WINDOW, linestyle="-", label="A_three_min")
    ax.plot(x[mask], A_max[mask], color=COLOR_WINDOW, linestyle="--", label="A_three_max")

    ax.set_xlabel(param_label)
    ax.set_ylabel("Bistable A-window")
    ax.set_ylim(0.0, 1.5)
    ax.legend(frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)


# ============================================================
# MAIN: S4–S6 ROBUSTNESS SWEEPS
# ============================================================


if __name__ == "__main__":
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    stage_name = STAGE_NAMES["S4_S6_RobustnessSweeps"]
    run_dir = make_stage_run_dir(stage_name, timestamp)

    p_base = PDParamsEC3()

    # -------------------------------
    # S3.1: Variation in k_M (mitochondrial turnover)
    # -------------------------------
    kM_values = np.linspace(0.4, 1.2, 25)
    res_kM = compute_bistable_window_vs_param("k_M", kM_values, p_base)
    f_s4 = os.path.join(run_dir, "supp_S4_kM_bistable_window.png")
    plot_window_vs_param(res_kM, param_label="k_M (mitochondrial turnover)", filename=f_s4, panel_label="S4")

    # -------------------------------
    # S3.2: Variation in L1 (load-dependent consumption)
    # -------------------------------
    L1_values = np.linspace(0.4, 1.4, 25)
    res_L1 = compute_bistable_window_vs_param("L1", L1_values, p_base)
    f_s5 = os.path.join(run_dir, "supp_S5_L1_bistable_window.png")
    plot_window_vs_param(res_L1, param_label="L_1 (load-dependent cost)", filename=f_s5, panel_label="S5")

    # -------------------------------
    # S3.3: Variation in C (Ca2+ handling load)
    # -------------------------------
    C_values = np.linspace(0.5, 1.5, 25)
    res_C = compute_bistable_window_vs_param("C", C_values, p_base)
    f_s6 = os.path.join(run_dir, "supp_S6_C_bistable_window.png")
    plot_window_vs_param(res_C, param_label="C (Ca$^{2+}$ load)", filename=f_s6, panel_label="S6")

    # -------------------------------
    # README for this run
    # -------------------------------
    write_readme(
        run_dir,
        header="Stage S4–S6: EC3 Robustness Sweeps",
        body_lines=[
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- supp_S4_kM_bistable_window.png",
            "- supp_S5_L1_bistable_window.png",
            "- supp_S6_C_bistable_window.png",
            "",
            "Description:",
            "Supplementary S4–S6 robustness sweeps for the EC3 model, showing",
            "how the bistable A-window shifts or persists under variation in",
            "(i) mitochondrial turnover k_M, (ii) load-dependent energy cost L1,",
            "and (iii) Ca2+ handling load C.",
        ],
    )

    print("Done.")
    print("Stage S4–S6 run directory:", run_dir)
