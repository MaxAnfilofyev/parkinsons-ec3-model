import numpy as np
from dataclasses import dataclass
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib.pyplot as plt
import csv
import textwrap

# ============================================================
# EC3 PARAMETERS (CALIBRATED FOR TRUE BISTABILITY)
# ============================================================

@dataclass
class PDParamsEC3:
    # Energy regeneration via mitochondria
    k1: float = 0.323525842503625   # strength of M-dependent energy restoration

    # Nonlinear positive feedback on energy (e.g. activity-dependent pump/repair)
    k2: float = 5.764702398127752   # strength of E^2(1-E) term

    # Baseline + load-dependent energy drain
    L0: float = 0.8329793043013896  # baseline leak
    L1: float = 0.7138379754309159  # extra drain per unit A*C

    # Mitochondrial turnover + load-dependent damage
    k_M: float = 0.7138550347382676 # repair / replacement rate
    beta: float = 1.5444631099298312# damage coupling to A*C*(1-E)

    # Structural / Ca2+ loads
    A: float = 0.8825917266901168   # axonal arbor load
    C: float = 1.0                  # Ca2+ handling load

    def with_A(self, A_new: float):
        """Return a copy of this parameter set with modified A."""
        d = self.__dict__.copy()
        d["A"] = A_new
        return PDParamsEC3(**d)


# ============================================================
# EC3 ODEs
# ============================================================

def dEdt(E: float, M: float, p: PDParamsEC3) -> float:
    """
    dE/dt = k1*M*(1 - E)       # mitochondrial-powered restoration
           + k2*E^2*(1 - E)    # nonlinear positive feedback on energy
           - (L0 + L1*A*C)*E   # baseline + load-dependent drain
    """
    return (
        p.k1 * M * (1.0 - E)
        + p.k2 * (E ** 2) * (1.0 - E)
        - (p.L0 + p.L1 * p.A * p.C) * E
    )


def dMdt(E: float, M: float, p: PDParamsEC3) -> float:
    """
    dM/dt = k_M*(1 - M)             # repair / replacement
           - beta*A*C*M*(1 - E)     # damage when loads high and energy low
    """
    return (
        p.k_M * (1.0 - M)
        - p.beta * p.A * p.C * M * (1.0 - E)
    )


def ec3_ode(t, y, p: PDParamsEC3):
    E, M = y
    return np.array([dEdt(E, M, p), dMdt(E, M, p)])


# ============================================================
# EQUILIBRIA + STABILITY
# ============================================================

def jacobian(E: float, M: float, p: PDParamsEC3) -> np.ndarray:
    """
    Jacobian matrix of the EC3 system at (E, M).
    """
    # dEdt = k1*M*(1-E) + k2*E^2*(1-E) - (L0 + L1*A*C)*E
    # Partial derivatives:
    dEdE = (
        p.k1 * M * (-1.0)
        + p.k2 * (2.0 * E * (1.0 - E) + (E ** 2) * (-1.0))  # d/dE of E^2(1-E)
        - (p.L0 + p.L1 * p.A * p.C)
    )
    dEdM = p.k1 * (1.0 - E)

    # dMdt = k_M*(1-M) - beta*A*C*M*(1-E)
    dMdE = p.beta * p.A * p.C * M   # derivative of -beta*A*C*M*(1-E) wrt E
    dMdM = -p.k_M - p.beta * p.A * p.C * (1.0 - E)

    return np.array([[dEdE, dEdM],
                     [dMdE, dMdM]])


def classify_equilibrium(E: float, M: float, p: PDParamsEC3):
    """
    Classify equilibrium by linear stability:
    - 'stable'  (both eigenvalues have negative real part)
    - 'saddle'  (eigenvalues with opposite sign real parts)
    - 'unstable' (both positive or complex with positive real part)
    """
    J = jacobian(E, M, p)
    eigvals = np.linalg.eigvals(J)
    reals = [ev.real for ev in eigvals]

    if all(r < 0 for r in reals):
        label = "stable"
    elif any(r > 0 for r in reals) and any(r < 0 for r in reals):
        label = "saddle"
    else:
        label = "unstable"

    return label, eigvals


def find_equilibria(p: PDParamsEC3, ngrid: int = 7, tol: float = 1e-3):
    """
    Grid-based multi-start root search for equilibria in [0,1]x[0,1].
    Returns a list of (E, M, label) tuples with duplicates removed.
    """
    seeds_E = np.linspace(0.05, 0.95, ngrid)
    seeds_M = np.linspace(0.05, 0.95, ngrid)

    solutions = []

    for E0 in seeds_E:
        for M0 in seeds_M:
            try:
                sol = root(
                    lambda x: np.array([
                        dEdt(x[0], x[1], p),
                        dMdt(x[0], x[1], p)
                    ]),
                    np.array([E0, M0]),
                    method="hybr",
                )
            except Exception:
                continue

            if not sol.success:
                continue

            E_star, M_star = sol.x
            # Require physical domain
            if not (-1e-4 <= E_star <= 1.0 + 1e-4 and -1e-4 <= M_star <= 1.0 + 1e-4):
                continue

            # Deduplicate
            if any(
                (abs(E_star - e) < tol and abs(M_star - m) < tol)
                for (e, m, _) in solutions
            ):
                continue

            label, _ = classify_equilibrium(E_star, M_star, p)
            solutions.append((E_star, M_star, label))

    return solutions


# ============================================================
# BIFURCATION SCAN VS A
# ============================================================

def scan_equilibria_vs_A(
    p_base: PDParamsEC3,
    A_min: float = 0.2,
    A_max: float = 1.4,
    nA: int = 61,
    ngrid: int = 7,
):
    """
    For a range of A values, compute all equilibria and their stability.
    Returns:
      data_rows: list of dicts with keys ["A","E","M","label"]
      summary: dict with high-level info (A ranges for 3 equilibria, etc.)
    """
    As = np.linspace(A_min, A_max, nA)
    data_rows = []
    counts = []

    for A in As:
        pA = p_base.with_A(A)
        eqs = find_equilibria(pA, ngrid=ngrid)
        counts.append(len(eqs))

        for (E_star, M_star, label) in eqs:
            data_rows.append({
                "A": float(A),
                "E": float(E_star),
                "M": float(M_star),
                "label": label,
            })

    # Identify A range where we see 3 equilibria (bistable window)
    A_three = [A for A, c in zip(As, counts) if c == 3]
    if len(A_three) > 0:
        A_three_min = float(np.min(A_three))
        A_three_max = float(np.max(A_three))
    else:
        A_three_min = None
        A_three_max = None

    summary = {
        "A_min": float(A_min),
        "A_max": float(A_max),
        "nA": int(nA),
        "ngrid": int(ngrid),
        "A_values_with_3_equilibria": [float(a) for a in A_three],
        "A_three_min": A_three_min,
        "A_three_max": A_three_max,
    }

    return data_rows, summary


# ============================================================
# PLOTTING
# ============================================================

def plot_bifurcation_E_vs_A(data_rows, filename: str):
    """
    Plot E* vs A, with markers colored by stability.
    """
    A_stable = []
    E_stable = []
    A_saddle = []
    E_saddle = []
    A_unstable = []
    E_unstable = []

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

    plt.figure(figsize=(6, 4))
    if A_stable:
        plt.scatter(A_stable, E_stable, s=20, marker="o", label="stable")
    if A_saddle:
        plt.scatter(A_saddle, E_saddle, s=20, marker="x", label="saddle")
    if A_unstable:
        plt.scatter(A_unstable, E_unstable, s=20, marker="^", label="unstable")

    plt.xlabel("A (axonal load)")
    plt.ylabel("E* (steady-state energy)")
    plt.title("EC3 Bifurcation Diagram: E* vs A")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def make_phase_plane_plot(p: PDParamsEC3, A_value: float, filename: str):
    """
    Generate a phase-plane plot at a given A:
    - nullclines (dE/dt=0, dM/dt=0)
    - vector field
    - a few trajectories
    - equilibrium markers with stability labels
    """
    pA = p.with_A(A_value)

    # Grid for vector field and nullclines
    E_vals = np.linspace(0, 1, 25)
    M_vals = np.linspace(0, 1, 25)
    EE, MM = np.meshgrid(E_vals, M_vals)

    dE_grid = np.zeros_like(EE)
    dM_grid = np.zeros_like(MM)

    for i in range(EE.shape[0]):
        for j in range(EE.shape[1]):
            dE_grid[i, j] = dEdt(EE[i, j], MM[i, j], pA)
            dM_grid[i, j] = dMdt(EE[i, j], MM[i, j], pA)

    speed = np.sqrt(dE_grid**2 + dM_grid**2)
    speed[speed == 0] = 1e-8  # avoid division by zero
    dE_unit = dE_grid / speed
    dM_unit = dM_grid / speed

    plt.figure(figsize=(6, 5))

    # Vector field (unit arrows for direction only)
    plt.quiver(EE, MM, dE_unit, dM_unit, alpha=0.4, linewidth=0.3)

    # Nullclines
    try:
        cs1 = plt.contour(
            EE, MM, dE_grid,
            levels=[0.0],
            linestyles="solid",
            linewidths=1.5
        )
        cs1.collections[0].set_label("dE/dt = 0")
    except Exception:
        pass

    try:
        cs2 = plt.contour(
            EE, MM, dM_grid,
            levels=[0.0],
            linestyles="dashed",
            linewidths=1.5
        )
        cs2.collections[0].set_label("dM/dt = 0")
    except Exception:
        pass

    # Equilibria at this A
    eqs = find_equilibria(pA, ngrid=9)
    for (E_star, M_star, label) in eqs:
        if label == "stable":
            plt.scatter(E_star, M_star, c="black", s=40, marker="o", label="stable eq")
        elif label == "saddle":
            plt.scatter(E_star, M_star, c="red", s=40, marker="x", label="saddle eq")
        else:
            plt.scatter(E_star, M_star, c="grey", s=40, marker="^", label="unstable eq")

    # A few trajectories from different initial conditions
    initials = [
        (0.05, 0.1),
        (0.2, 0.8),
        (0.8, 0.2),
        (0.95, 0.9),
        (0.4, 0.4),
        (0.7, 0.7),
    ]

    for (E0, M0) in initials:
        sol = solve_ivp(
            lambda t, y: ec3_ode(t, y, pA),
            t_span=(0, 200),
            y0=[E0, M0],
            max_step=0.5,
            dense_output=False,
        )
        plt.plot(sol.y[0, :], sol.y[1, :], linewidth=0.8)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel("E (energy)")
    plt.ylabel("M (mitochondrial capacity)")
    plt.title(f"EC3 Phase Plane at A = {A_value:.2f}")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def make_timecourse_plot(p: PDParamsEC3, A_value: float, filename: str):
    """
    Time courses at a bistable A for two different initial conditions,
    showing convergence to different attractors.
    """
    pA = p.with_A(A_value)

    initials = {
        "IC1_low": [0.05, 0.2],
        "IC2_high": [0.9, 0.9],
    }

    t_eval = np.linspace(0, 300, 2000)

    plt.figure(figsize=(7, 5))

    for label, y0 in initials.items():
        sol = solve_ivp(
            lambda t, y: ec3_ode(t, y, pA),
            t_span=(0, t_eval[-1]),
            y0=y0,
            t_eval=t_eval,
            max_step=0.5,
        )
        plt.plot(
            sol.t, sol.y[0, :],
            label=f"E(t), {label}"
        )

    plt.xlabel("t (arbitrary units)")
    plt.ylabel("E(t)")
    plt.title(f"EC3 Time Courses at A = {A_value:.2f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# ============================================================
# DATA / SUMMARY WRITERS
# ============================================================

def save_equilibria_csv(data_rows, filename: str):
    """
    Save equilibria data to CSV: columns A,E,M,label.
    """
    fieldnames = ["A", "E", "M", "label"]
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in data_rows:
            writer.writerow(row)


def save_summary_txt(summary: dict, filename: str):
    """
    Save a human-readable summary of the bifurcation scan.
    """
    lines = []
    lines.append("EC3 Bifurcation Scan Summary")
    lines.append("=" * 40)
    lines.append(f"A_min = {summary['A_min']}")
    lines.append(f"A_max = {summary['A_max']}")
    lines.append(f"nA    = {summary['nA']}")
    lines.append(f"ngrid = {summary['ngrid']}")
    lines.append("")
    lines.append("A values with 3 equilibria (bistable window):")
    if summary["A_values_with_3_equilibria"]:
        lines.append(
            textwrap.fill(
                ", ".join(f"{a:.3f}" for a in summary["A_values_with_3_equilibria"]),
                width=78,
            )
        )
        lines.append("")
        lines.append(f"A_three_min = {summary['A_three_min']}")
        lines.append(f"A_three_max = {summary['A_three_max']}")
    else:
        lines.append("None (no A with 3 equilibria found).")

    with open(filename, "w") as f:
        f.write("\n".join(lines))


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == "__main__":
    import os, datetime

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = f"Stage_2_Bifurcation/runs/{timestamp}"
    os.makedirs(run_dir, exist_ok=True)

    def out(path):
        return os.path.join(run_dir, path)

    
    # Base parameter set that exhibits true bistability vs A
    p_base = PDParamsEC3()

    # 1) Scan equilibria vs A and classify stability
    data_rows, summary = scan_equilibria_vs_A(
        p_base,
        A_min=0.2,
        A_max=1.4,
        nA=61,
        ngrid=7,
    )

    # 2) Save data and summary
    save_equilibria_csv(data_rows, out("ec3_equilibria_vs_A.csv"))
    save_summary_txt(summary, out("ec3_bifurcation_summary.txt"))

    # 3) Bifurcation-style plot E*(A) with stability coding
    plot_bifurcation_E_vs_A(data_rows, out("ec3_bifurcation_E_vs_A.png"))

    # 4) Phase-plane plot at a representative bistable A
    #    If the scan found a bistable window, choose its midpoint; else fallback.
    if summary["A_three_min"] is not None:
        A_demo = 0.5 * (summary["A_three_min"] + summary["A_three_max"])
    else:
        A_demo = p_base.A  # fallback (but for this calibrated set, we DO have a window)

    make_phase_plane_plot(p_base, A_demo, out(f"ec3_phase_plane_A_{A_demo:.2f}.png"))

    # 5) Time courses at the same A_demo to show different attractors
    make_timecourse_plot(p_base, A_demo, out(f"ec3_timecourses_A_{A_demo:.2f}.png"))

    print("Done.")
    print("Generated files:")
    print("  - ec3_equilibria_vs_A.csv")
    print("  - ec3_bifurcation_summary.txt")
    print("  - ec3_bifurcation_E_vs_A.png")
    print(f"  - ec3_phase_plane_A_{A_demo:.2f}.png")
    print(f"  - ec3_timecourses_A_{A_demo:.2f}.png")
    with open(out("README.txt"), "w") as f:
        f.write(f"""
            EC3 Run
            Timestamp: {timestamp}

            Files:
            - bifurcation_E_vs_A.png
            - ec3_equilibria_vs_A.csv
            - ec3_bifurcation_summary.txt

            Parameters:
            {p_base}

            Interpretation notes:
            (Add short notes here manually later)
            """)
    
    # Update or create a latest_run symlink for this stage directory
    stage_dir = os.path.dirname(os.path.dirname(run_dir))
    latest = os.path.join(stage_dir, "latest_run")

    if os.path.islink(latest) or os.path.exists(latest):
        os.remove(latest)

    os.symlink(run_dir, latest, target_is_directory=True)
