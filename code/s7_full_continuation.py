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

COLOR_NULLCLINE = "#4b0082"  # dark purple for branches


def new_panel(figsize: Tuple[float, float] = (3.3, 2.2), panel_label: str | None = None):
    """Create a new figure/axes pair with optional panel label.

    Mirrors the helper from `run_ec3_all.py` to keep sizing consistent.
    """
    fig, ax = plt.subplots(figsize=figsize)
    if panel_label is not None:
        ax.text(0.02, 0.95, panel_label, transform=ax.transAxes,
                fontweight="bold", va="top")
    return fig, ax


def plot_full_bifurcation_E_vs_A(data_rows, filename: str, A_phys_max: float = 1.4):
    """Plot full continuation of E* vs A over an extended A-range.

    Stable, saddle, and unstable equilibria are shown with consistent
    colors and markers. The physiologically relevant range of A (up to
    `A_phys_max`) is shaded lightly to emphasize that the right-hand
    fold lies beyond it.
    """
    A_stable, E_stable = [], []
    A_saddle, E_saddle = [], []
    A_unstable, E_unstable = [], []

    for row in data_rows:
        if row["label"] == "stable":
            A_stable.append(row["A"])
            E_stable.append(row["E"])
        elif row["label"] == "saddle":
            A_saddle.append(row["A"])
            E_saddle.append(row["E"])
        else:
            A_unstable.append(row["A"])
            E_unstable.append(row["E"])

    fig, ax = new_panel(figsize=(5.0, 2.2), panel_label="S7")

    # Use high-contrast colors and slightly larger markers for clarity
    if A_stable:
        ax.scatter(A_stable, E_stable, s=22, marker="o",
                   color="#1f77b4", label="stable")
    if A_unstable:
        ax.scatter(A_unstable, E_unstable, s=22, marker="o",
                   facecolors="none", edgecolors="#ff7f0e",
                   linewidths=1.0, label="unstable")
    if A_saddle:
        ax.scatter(A_saddle, E_saddle, s=28, marker="x",
                   color="red", label="saddle")

    # Shade physiologically relevant A-range for dopaminergic neurons
    A_min = min(row["A"] for row in data_rows)
    A_max = max(row["A"] for row in data_rows)
    ax.axvspan(A_min, A_phys_max, color="#fff2b2", alpha=0.5,
               label="physiological A-range")

    ax.set_xlabel("Axonal load A")
    ax.set_ylabel("Steady-state energy E*")
    ax.set_xlim(A_min, A_max)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Full continuation of E* vs A (including right fold)")

    # Place legend to the right with markers on the right and right-aligned labels
    legend = ax.legend(
        frameon=False,
        loc="upper right",
        markerfirst=False,
    )
    for txt in legend.get_texts():
        txt.set_ha("right")

    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.22, top=0.88)
    fig.savefig(filename, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    stage_name = STAGE_NAMES["S7_FullContinuation"]
    run_dir = make_stage_run_dir(stage_name, timestamp)

    p_base = PDParamsEC3()

    # Extend A beyond the physiological range to reveal the right-hand fold
    A_min = 0.2
    A_max = 2.0
    nA = 90

    data_rows, summary = scan_equilibria_vs_A(
        p_base,
        A_min=A_min,
        A_max=A_max,
        nA=nA,
        ngrid=7,
    )

    out_png = os.path.join(run_dir, "supp_S7_full_bifurcation_E_vs_A.png")
    plot_full_bifurcation_E_vs_A(data_rows, out_png, A_phys_max=1.4)

    # README
    write_readme(
        run_dir,
        header="Stage S7: Full Continuation of EC3 Bifurcation vs A",
        body_lines=[
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- supp_S7_full_bifurcation_E_vs_A.png",
            "",
            "Description:",
            "Supplementary S7 figure showing the full continuation of E* vs A ",
            "for the EC3 model, including the right-hand fold located at ",
            "axonal loads beyond the physiological range of dopaminergic ",
            "neurons.",
        ],
    )

    print("Done.")
    print("Stage S7 run directory:", run_dir)
