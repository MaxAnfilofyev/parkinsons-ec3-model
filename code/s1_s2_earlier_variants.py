import os
import datetime
from dataclasses import dataclass
from typing import Tuple, List

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib as mpl
import matplotlib.pyplot as plt

from stage_utils import make_stage_run_dir, write_readme, STAGE_NAMES

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
# PARAMETERS AND DYNAMICS FOR EARLIER MODEL VARIANTS (S1/S2)
# ============================================================

@dataclass
class PDParamsEarly:
    """Parameter set for earlier (pre-EC3) energetic models.

    These are simplified relatives of the EC3 parameters, chosen to
    be qualitatively consistent with the Supplementary text:

    - Variant 1: linear mitochondrial support, linear energy drain,
      no nonlinear energy amplification.
    - Variant 2: feedback-only model, where energetic deficit modulates
      the energy equation via a (1-M)*E term, but still lacks the
      nonlinear E^2(1-E) amplification.
    """

    # Energy / mitochondria parameters
    k1: float = 0.5        # mitochondrial support term
    L0: float = 0.1        # baseline energy consumption
    L1: float = 1.0        # load-dependent energy consumption
    k_M: float = 0.5       # mitochondrial repair / turnover
    beta: float = 0.8      # load- and energy-dependent mito damage

    # Feedback-only extra term (variant 2)
    k_fb: float = 0.3      # strength of (1-M)*E feedback in dE/dt

    # Structural / Ca2+ loads
    A: float = 0.8         # axonal load
    C: float = 1.0         # Ca2+ load

    def with_A(self, A_new: float) -> "PDParamsEarly":
        d = self.__dict__.copy()
        d["A"] = A_new
        return PDParamsEarly(**d)


# ---------- Variant 1: linear support, linear load, no nonlinear amplification


def dEdt_variant1(E: float, M: float, p: PDParamsEarly) -> float:
    """Energy dynamics for the first earlier variant (S1).

    - Linear mitochondrial support: k1 * M * (1 - E)
    - Linear energy consumption: (L0 + L1*A*C) * E
    - No nonlinear E^2(1-E) amplification term.
    """
    return p.k1 * M * (1.0 - E) - (p.L0 + p.L1 * p.A * p.C) * E


def dMdt_common(E: float, M: float, p: PDParamsEarly) -> float:
    """Shared mitochondrial dynamics used in both variants.

    Mitochondria recover toward M = 1 and are damaged when load is
    high and energy is low.
    """
    return p.k_M * (1.0 - M) - p.beta * p.A * p.C * M * (1.0 - E)


def ode_variant1(t: float, y: np.ndarray, p: PDParamsEarly) -> np.ndarray:
    E, M = y
    return np.array([dEdt_variant1(E, M, p), dMdt_common(E, M, p)])


# ---------- Variant 2: feedback-only model (adds (1-M)*E term) ----------


def dEdt_variant2(E: float, M: float, p: PDParamsEarly) -> float:
    """Energy dynamics for the second earlier variant (S2).

    Adds a feedback term of the form (1 - M) * E, as described in the
    Supplementary text, while retaining linear mitochondrial support
    and linear load-dependent consumption.
    """
    base = dEdt_variant1(E, M, p)
    feedback = p.k_fb * (1.0 - M) * E
    return base + feedback


def ode_variant2(t: float, y: np.ndarray, p: PDParamsEarly) -> np.ndarray:
    E, M = y
    return np.array([dEdt_variant2(E, M, p), dMdt_common(E, M, p)])


# ============================================================
# EQUILIBRIA AND STABILITY HELPERS (REUSED FOR BOTH VARIANTS)
# ============================================================


def jacobian_variant1(E: float, M: float, p: PDParamsEarly) -> np.ndarray:
    """Jacobian for variant 1 at (E, M)."""
    # dEdt = k1*M*(1-E) - (L0 + L1*A*C)*E
    dEdE = p.k1 * M * (-1.0) - (p.L0 + p.L1 * p.A * p.C)
    dEdM = p.k1 * (1.0 - E)

    # dMdt = k_M*(1-M) - beta*A*C*M*(1-E)
    dMdE = p.beta * p.A * p.C * M
    dMdM = -p.k_M - p.beta * p.A * p.C * (1.0 - E)

    return np.array([[dEdE, dEdM], [dMdE, dMdM]])


def jacobian_variant2(E: float, M: float, p: PDParamsEarly) -> np.ndarray:
    """Jacobian for variant 2 at (E, M)."""
    # dEdt_v2 = dEdt_v1 + k_fb*(1-M)*E
    # For the feedback term: f_fb = k_fb*(1-M)*E
    # df_fb/dE = k_fb*(1-M)
    # df_fb/dM = -k_fb*E
    J1 = jacobian_variant1(E, M, p)
    dEdE_fb = p.k_fb * (1.0 - M)
    dEdM_fb = -p.k_fb * E

    dEdE = J1[0, 0] + dEdE_fb
    dEdM = J1[0, 1] + dEdM_fb

    # dMdt is unchanged
    dMdE = J1[1, 0]
    dMdM = J1[1, 1]

    return np.array([[dEdE, dEdM], [dMdE, dMdM]])


def classify_equilibrium(E: float, M: float, p: PDParamsEarly, variant: int) -> str:
    """Classify equilibrium stability for a given variant.

    Returns 'stable', 'saddle', or 'unstable'.
    """
    if variant == 1:
        J = jacobian_variant1(E, M, p)
    else:
        J = jacobian_variant2(E, M, p)

    eigvals = np.linalg.eigvals(J)
    reals = [ev.real for ev in eigvals]

    if all(r < 0 for r in reals):
        return "stable"
    elif any(r > 0 for r in reals) and any(r < 0 for r in reals):
        return "saddle"
    else:
        return "unstable"


def find_equilibria(p: PDParamsEarly, variant: int, ngrid: int = 7, tol: float = 1e-3) -> List[Tuple[float, float, str]]:
    """Grid-based multi-start root search for equilibria in [0,1]x[0,1].

    Returns a list of (E, M, label) tuples, deduplicated within `tol`.
    """
    seeds_E = np.linspace(0.05, 0.95, ngrid)
    seeds_M = np.linspace(0.05, 0.95, ngrid)

    solutions: List[Tuple[float, float, str]] = []

    for E0 in seeds_E:
        for M0 in seeds_M:
            try:
                if variant == 1:
                    sol = root(
                        lambda x: np.array([
                            dEdt_variant1(x[0], x[1], p),
                            dMdt_common(x[0], x[1], p),
                        ]),
                        np.array([E0, M0]),
                        method="hybr",
                    )
                else:
                    sol = root(
                        lambda x: np.array([
                            dEdt_variant2(x[0], x[1], p),
                            dMdt_common(x[0], x[1], p),
                        ]),
                        np.array([E0, M0]),
                        method="hybr",
                    )
            except Exception:
                continue

            if not sol.success:
                continue

            E_star, M_star = sol.x
            if not (-1e-4 <= E_star <= 1.0 + 1e-4 and -1e-4 <= M_star <= 1.0 + 1e-4):
                continue

            if any(
                (abs(E_star - e) < tol and abs(M_star - m) < tol)
                for (e, m, _) in solutions
            ):
                continue

            label = classify_equilibrium(E_star, M_star, p, variant)
            solutions.append((E_star, M_star, label))

    return solutions


# ============================================================
# S1: NULLCLINES + VECTOR FIELD FOR MONOSTABLE VARIANT
# ============================================================


def plot_S1_nullclines_and_field(p: PDParamsEarly, A_value: float, filename: str) -> None:
    """Generate Supplementary S1-style figure for variant 1.

    Shows:
    - Vector field in (E, M).
    - dE/dt = 0 (solid purple) and dM/dt = 0 (dashed purple) nullclines.
    - Single equilibrium marker, expected to be a globally attracting point.
    """
    pA = p.with_A(A_value)

    E_vals = np.linspace(0.0, 1.0, 50)
    M_vals = np.linspace(0.0, 1.0, 50)
    EE, MM = np.meshgrid(E_vals, M_vals)

    dE_grid = np.zeros_like(EE)
    dM_grid = np.zeros_like(MM)

    for i in range(EE.shape[0]):
        for j in range(EE.shape[1]):
            dE_grid[i, j] = dEdt_variant1(EE[i, j], MM[i, j], pA)
            dM_grid[i, j] = dMdt_common(EE[i, j], MM[i, j], pA)

    speed = np.sqrt(dE_grid**2 + dM_grid**2)
    speed[speed == 0] = 1e-8
    dE_unit = dE_grid / speed
    dM_unit = dM_grid / speed

    fig, ax = new_panel(figsize=(3.0, 3.0), panel_label="S1")

    # Vector field
    step = 2
    ax.quiver(
        EE[::step, ::step],
        MM[::step, ::step],
        dE_unit[::step, ::step],
        dM_unit[::step, ::step],
        color="gray",
        alpha=0.5,
        linewidth=0.3,
    )

    # Nullclines
    try:
        ax.contour(
            EE,
            MM,
            dE_grid,
            levels=[0.0],
            linestyles="solid",
            linewidths=1.5,
            colors=COLOR_NULLCLINE,
        )
    except Exception:
        pass

    try:
        ax.contour(
            EE,
            MM,
            dM_grid,
            levels=[0.0],
            linestyles="dashed",
            linewidths=1.5,
            colors=COLOR_NULLCLINE,
        )
    except Exception:
        pass

    # Equilibrium(s)
    eqs = find_equilibria(pA, variant=1, ngrid=9)
    for (E_star, M_star, label) in eqs:
        if label == "stable":
            ax.scatter(E_star, M_star, c="black", s=40, marker="o")
        elif label == "saddle":
            ax.scatter(E_star, M_star, c="red", s=50, marker="x")
        else:
            ax.scatter(E_star, M_star, c="grey", s=40, marker="^")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("E (energy)")
    ax.set_ylabel("M (mitochondrial capacity)")
    ax.set_title(f"Earlier Variant 1 (monostable), A = {A_value:.2f}")

    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)


# ============================================================
# S2: BIFURCATION-STYLE SCAN FOR FEEDBACK-ONLY VARIANT
# ============================================================


def scan_equilibria_vs_A_variant2(
    p_base: PDParamsEarly,
    A_min: float = 0.1,
    A_max: float = 2.0,
    nA: int = 50,
    ngrid: int = 7,
):
    """Scan equilibria vs A for the feedback-only model (variant 2).

    Returns a list of (A, E, M, label) dicts. For this architecture we
    expect a single equilibrium (monostable) across the scanned range.
    """
    As = np.linspace(A_min, A_max, nA)
    data_rows = []

    for A in As:
        pA = p_base.with_A(A)
        eqs = find_equilibria(pA, variant=2, ngrid=ngrid)
        for (E_star, M_star, label) in eqs:
            data_rows.append({
                "A": float(A),
                "E": float(E_star),
                "M": float(M_star),
                "label": label,
            })

    return data_rows


def plot_S2_bifurcation_variant2(data_rows, filename: str) -> None:
    """Plot E* vs A for the feedback-only model (S2).

    Stable, saddle, and unstable equilibria are shown with different
    markers, but in practice we expect only a single stable branch.
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

    fig, ax = new_panel(figsize=(3.3, 2.0), panel_label="S2")

    if A_stable:
        ax.scatter(A_stable, E_stable, s=18, marker="o",
                   color=COLOR_NULLCLINE, label="stable")
    if A_unstable:
        ax.scatter(A_unstable, E_unstable, s=18, marker="o",
                   facecolors="none", edgecolors=COLOR_NULLCLINE,
                   label="unstable")
    if A_saddle:
        ax.scatter(A_saddle, E_saddle, s=24, marker="x",
                   color="red", label="saddle")

    ax.set_xlabel("Axonal load A")
    ax.set_ylabel("Steady-state energy E*")
    ax.set_title("Earlier Variant 2 (feedback-only): E* vs A")
    ax.set_xlim(min(row["A"] for row in data_rows),
                max(row["A"] for row in data_rows))
    ax.set_ylim(0.0, 1.0)

    ax.legend(frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)


# ============================================================
# MAIN: GENERATE S1 AND S2 FIGURES
# ============================================================


if __name__ == "__main__":
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Create a run directory for S1/S2 using the shared stage utility.
    stage_name = STAGE_NAMES["S1_S2_EarlierVariants"]
    run_dir = make_stage_run_dir(stage_name, timestamp)

    p_base = PDParamsEarly()

    # -------------------------------
    # S1: Nullclines + vector field
    # -------------------------------
    A_demo = 1.0  # SNc-like load level used for illustration
    s1_filename = os.path.join(run_dir, "supp_S1_nullclines_monostable.png")
    plot_S1_nullclines_and_field(p_base, A_demo, s1_filename)

    # -------------------------------
    # S2: Bifurcation-style scan
    # -------------------------------
    data_rows_s2 = scan_equilibria_vs_A_variant2(
        p_base,
        A_min=0.1,
        A_max=2.0,
        nA=60,
        ngrid=7,
    )
    s2_filename = os.path.join(run_dir, "supp_S2_bifurcation_feedback_only.png")
    plot_S2_bifurcation_variant2(data_rows_s2, s2_filename)

    # -------------------------------
    # README for this run
    # -------------------------------
    write_readme(
        run_dir,
        header="Stage S1–S2: Earlier Monostable and Feedback-Only Variants",
        body_lines=[
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- supp_S1_nullclines_monostable.png",
            "- supp_S2_bifurcation_feedback_only.png",
            "",
            "Description:",
            "Supplementary S1: nullclines and vector field for an earlier, ",
            "monostable energetic model without nonlinear energy amplification.",
            "Supplementary S2: bifurcation-style scan for a feedback-only ",
            "model lacking nonlinear energy restoration, showing a single ",
            "equilibrium branch across a wide range of axonal loads A.",
        ],
    )

    print("Done.")
    print("Stage S1–S2 run directory:", run_dir)
