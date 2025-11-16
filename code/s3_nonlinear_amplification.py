import os
import datetime
from dataclasses import dataclass
from typing import Tuple

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from stage_utils import make_stage_run_dir, write_readme, STAGE_NAMES
from run_ec3_all import PDParamsEC3

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

COLOR_NULLCLINE = "#4b0082"  # dark purple for nullclines


def new_panel(figsize: Tuple[float, float] = (3.3, 2.0), panel_label: str | None = None):
    """Create a new figure/axes pair with optional panel label.

    Mirrors the helper from `run_ec3_all.py` to keep sizing consistent.
    """
    fig, ax = plt.subplots(figsize=figsize)
    if panel_label is not None:
        ax.text(0.02, 0.95, panel_label, transform=ax.transAxes,
                fontweight="bold", va="top")
    return fig, ax


# ============================================================
# EC3-LIKE DYNAMICS WITH OPTIONAL NONLINEAR AMPLIFICATION TERM
# ============================================================


def dEdt(E: float, M: float, p: PDParamsEC3, k2_override: float | None = None) -> float:
    """Energy dynamics with optional override of the nonlinear term k2.

    When `k2_override` is 0, this reduces to a purely linear mitochondrial
    support model, as discussed in Supplementary S2.
    """
    k2 = p.k2 if k2_override is None else k2_override
    return (
        p.k1 * M * (1.0 - E)
        + k2 * (E**2) * (1.0 - E)
        - (p.L0 + p.L1 * p.A * p.C) * E
    )


def dMdt(E: float, M: float, p: PDParamsEC3) -> float:
    """Mitochondrial dynamics (same as EC3)."""
    return (
        p.k_M * (1.0 - M)
        - p.beta * p.A * p.C * M * (1.0 - E)
    )


# ============================================================
# S3 FIGURE: EFFECT OF NONLINEAR ENERGY AMPLIFICATION ON NULLCLINES
# ============================================================


def plot_S3_nullcline_comparison(p: PDParamsEC3, A_value: float, filename: str) -> None:
    """Generate Supplementary S3-style figure.

    Left panel: k2 = 0 (no nonlinear energy amplification) – energy nullcline
    is monotonic.
    Right panel: k2 = p.k2 (full EC3 nonlinearity) – energy nullcline is
    folded, producing up to three intersections with the mitochondrial
    nullcline.
    """
    # Fix load at an SNc-like value inside the bistable window.
    pA = p.with_A(A_value)

    E_vals = np.linspace(0.0, 1.0, 60)
    M_vals = np.linspace(0.0, 1.0, 60)
    EE, MM = np.meshgrid(E_vals, M_vals)

    # --- Panel A: k2 = 0 ---
    dE_grid_lin = np.zeros_like(EE)
    dM_grid_lin = np.zeros_like(MM)

    for i in range(EE.shape[0]):
        for j in range(EE.shape[1]):
            dE_grid_lin[i, j] = dEdt(EE[i, j], MM[i, j], pA, k2_override=0.0)
            dM_grid_lin[i, j] = dMdt(EE[i, j], MM[i, j], pA)

    # --- Panel B: k2 = p.k2 ---
    dE_grid_nonlin = np.zeros_like(EE)
    dM_grid_nonlin = np.zeros_like(MM)

    for i in range(EE.shape[0]):
        for j in range(EE.shape[1]):
            dE_grid_nonlin[i, j] = dEdt(EE[i, j], MM[i, j], pA, k2_override=None)
            dM_grid_nonlin[i, j] = dMdt(EE[i, j], MM[i, j], pA)

    # --- Create figure with two panels ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.4, 3.0), sharex=True, sharey=True)

    # Panel A: linear model (k2 = 0)
    ax1.set_title("k2 = 0 (no nonlinear amplification)")
    try:
        ax1.contour(
            EE,
            MM,
            dE_grid_lin,
            levels=[0.0],
            linestyles="solid",
            linewidths=1.5,
            colors=COLOR_NULLCLINE,
        )
    except Exception:
        pass
    try:
        ax1.contour(
            EE,
            MM,
            dM_grid_lin,
            levels=[0.0],
            linestyles="dashed",
            linewidths=1.5,
            colors=COLOR_NULLCLINE,
        )
    except Exception:
        pass
    ax1.set_xlabel("E (energy)")
    ax1.set_ylabel("M (mitochondrial capacity)")

    # Panel B: nonlinear model (k2 > 0)
    ax2.set_title(f"k2 = {pA.k2:.2f} (nonlinear amplification)")
    try:
        ax2.contour(
            EE,
            MM,
            dE_grid_nonlin,
            levels=[0.0],
            linestyles="solid",
            linewidths=1.5,
            colors=COLOR_NULLCLINE,
        )
    except Exception:
        pass
    try:
        ax2.contour(
            EE,
            MM,
            dM_grid_nonlin,
            levels=[0.0],
            linestyles="dashed",
            linewidths=1.5,
            colors=COLOR_NULLCLINE,
        )
    except Exception:
        pass
    ax2.set_xlabel("E (energy)")

    for ax in (ax1, ax2):
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)

    fig.suptitle("Effect of nonlinear energy amplification on nullcline geometry")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    fig.savefig(filename, dpi=300)
    plt.close(fig)


# ============================================================
# MAIN: GENERATE S3 FIGURE
# ============================================================


if __name__ == "__main__":
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    stage_name = STAGE_NAMES["S3_NonlinearAmplification"]
    run_dir = make_stage_run_dir(stage_name, timestamp)

    p_base = PDParamsEC3()

    # Use an SNc-like load inside the bistable regime
    A_demo = 1.0
    s3_filename = os.path.join(run_dir, "supp_S3_nonlinear_amplification_nullclines.png")

    plot_S3_nullcline_comparison(p_base, A_demo, s3_filename)

    write_readme(
        run_dir,
        header="Stage S3: Effect of Nonlinear Energy Amplification",
        body_lines=[
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- supp_S3_nonlinear_amplification_nullclines.png",
            "",
            "Description:",
            "Supplementary S3 figure comparing EC3-like nullcline geometry ",
            "with k2 = 0 (no nonlinear energy amplification) versus k2 > 0, ",
            "illustrating how the nonlinear term generates a folded energy ",
            "nullcline and enables multiple equilibria.",
        ],
    )

    print("Done.")
    print("Stage S3 run directory:", run_dir)
