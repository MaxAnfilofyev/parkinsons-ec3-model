import os
import csv
import textwrap
import datetime
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib.pyplot as plt


# ============================================================
# UTILITIES: FOLDERS, RUN DIRS, READMEs
# ============================================================

def make_stage_run_dir(stage_name: str, timestamp: str) -> str:
    """
    Create Stage_<...>/runs/<timestamp>/ and update latest_run link/text.
    Returns absolute path to the run directory.
    """
    root = os.getcwd()
    stage_dir = os.path.join(root, stage_name)
    runs_dir = os.path.join(stage_dir, "runs")

    os.makedirs(runs_dir, exist_ok=True)

    run_dir = os.path.join(runs_dir, timestamp)
    os.makedirs(run_dir, exist_ok=True)

    latest_link = os.path.join(stage_dir, "latest_run")
    # Remove existing symlink or file
    if os.path.islink(latest_link) or os.path.exists(latest_link):
        try:
            os.remove(latest_link)
        except OSError:
            pass

    # Try to create symlink; if not allowed, write a text file instead
    try:
        os.symlink(run_dir, latest_link, target_is_directory=True)
    except OSError:
        with open(os.path.join(stage_dir, "latest_run.txt"), "w") as f:
            f.write(run_dir)

    return run_dir


def write_readme(run_dir: str, header: str, body_lines):
    """
    Create a simple README.txt summarizing this run.
    """
    path = os.path.join(run_dir, "README.txt")
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("=" * len(header) + "\n\n")
        for line in body_lines:
            f.write(line + "\n")


# ============================================================
# EC3 PARAMETERS AND DYNAMICS
# ============================================================

@dataclass
class PDParamsEC3:
    # Energy regeneration via mitochondria
    k1: float = 0.323525842503625

    # Nonlinear positive feedback on energy
    k2: float = 5.764702398127752

    # Baseline + load-dependent energy drain
    L0: float = 0.8329793043013896
    L1: float = 0.7138379754309159

    # Mitochondrial turnover + load-dependent damage
    k_M: float = 0.7138550347382676
    beta: float = 1.5444631099298312

    # Structural / Ca2+ loads
    A: float = 0.8825917266901168
    C: float = 1.0

    def with_A(self, A_new: float):
        d = self.__dict__.copy()
        d["A"] = A_new
        return PDParamsEC3(**d)


def dEdt(E: float, M: float, p: PDParamsEC3) -> float:
    return (
        p.k1 * M * (1.0 - E)
        + p.k2 * (E**2) * (1.0 - E)
        - (p.L0 + p.L1 * p.A * p.C) * E
    )


def dMdt(E: float, M: float, p: PDParamsEC3) -> float:
    return (
        p.k_M * (1.0 - M)
        - p.beta * p.A * p.C * M * (1.0 - E)
    )


def ec3_ode(t, y, p: PDParamsEC3):
    E, M = y
    return np.array([dEdt(E, M, p), dMdt(E, M, p)])


# ============================================================
# EQUILIBRIA AND STABILITY
# ============================================================

def jacobian(E: float, M: float, p: PDParamsEC3) -> np.ndarray:
    dEdE = (
        p.k1 * M * (-1.0)
        + p.k2 * (2.0 * E * (1.0 - E) + (E**2) * (-1.0))
        - (p.L0 + p.L1 * p.A * p.C)
    )
    dEdM = p.k1 * (1.0 - E)

    dMdE = p.beta * p.A * p.C * M
    dMdM = -p.k_M - p.beta * p.A * p.C * (1.0 - E)

    return np.array([[dEdE, dEdM],
                     [dMdE, dMdM]])


def classify_equilibrium(E: float, M: float, p: PDParamsEC3):
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
    seeds_E = np.linspace(0.05, 0.95, ngrid)
    seeds_M = np.linspace(0.05, 0.95, ngrid)

    solutions = []

    for E0 in seeds_E:
        for M0 in seeds_M:
            try:
                sol = root(
                    lambda x: np.array([dEdt(x[0], x[1], p),
                                        dMdt(x[0], x[1], p)]),
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

            label, _ = classify_equilibrium(E_star, M_star, p)
            solutions.append((E_star, M_star, label))

    return solutions


def scan_equilibria_vs_A(
    p_base: PDParamsEC3,
    A_min: float = 0.2,
    A_max: float = 1.4,
    nA: int = 61,
    ngrid: int = 7,
):
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
# SAVERS
# ============================================================

def save_equilibria_csv(data_rows, filename: str):
    fieldnames = ["A", "E", "M", "label"]
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in data_rows:
            writer.writerow(row)


def save_summary_txt(summary: dict, filename: str):
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
# PLOTS: BIFURCATION, PHASE PLANES, TIMECOURSES
# ============================================================

def plot_bifurcation_E_vs_A(data_rows, filename: str):
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

    plt.figure(figsize=(7, 4))
    if A_stable:
        plt.scatter(A_stable, E_stable, s=20, marker="o", label="stable")
    if A_saddle:
        plt.scatter(A_saddle, E_saddle, s=30, marker="x", label="saddle")
    if A_unstable:
        plt.scatter(A_unstable, E_unstable, s=20, marker="^", label="unstable")

    plt.xlabel("A (axonal load)")
    plt.ylabel("E* (steady-state energy)")
    plt.title("EC3 Bifurcation Diagram: E* vs A")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def make_phase_plane_plot(p: PDParamsEC3, A_value: float, filename: str, title_suffix: str):
    pA = p.with_A(A_value)

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
    speed[speed == 0] = 1e-8
    dE_unit = dE_grid / speed
    dM_unit = dM_grid / speed

    plt.figure(figsize=(7, 6))
    plt.quiver(EE, MM, dE_unit, dM_unit, alpha=0.4, linewidth=0.3)

    try:
        cs1 = plt.contour(EE, MM, dE_grid, levels=[0.0],
                          linestyles="solid", linewidths=1.5)
        cs1.collections[0].set_label("dE/dt=0")
    except Exception:
        pass

    try:
        cs2 = plt.contour(EE, MM, dM_grid, levels=[0.0],
                          linestyles="dashed", linewidths=1.5)
        cs2.collections[0].set_label("dM/dt=0")
    except Exception:
        pass

    eqs = find_equilibria(pA, ngrid=9)
    for (E_star, M_star, label) in eqs:
        if label == "stable":
            plt.scatter(E_star, M_star, c="black", s=45, marker="o")
        elif label == "saddle":
            plt.scatter(E_star, M_star, c="red", s=50, marker="x")
        else:
            plt.scatter(E_star, M_star, c="grey", s=40, marker="^")

    initials = [
        (0.05, 0.2),
        (0.2, 0.8),
        (0.8, 0.2),
        (0.9, 0.9),
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
    plt.title(f"EC3 Phase Plane at A = {A_value:.2f} ({title_suffix})")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def simulate_timecourse(p: PDParamsEC3, A_value: float, y0, t_end=300, n_points=2000):
    pA = p.with_A(A_value)
    t_eval = np.linspace(0, t_end, n_points)
    sol = solve_ivp(
        lambda t, y: ec3_ode(t, y, pA),
        t_span=(0, t_end),
        y0=y0,
        t_eval=t_eval,
        max_step=0.5,
    )
    return sol.t, sol.y


def simulate_with_perturbation(p: PDParamsEC3,
                               A_value: float,
                               y0,
                               t_pre: float = 50.0,
                               t_post: float = 300.0,
                               E_after: float = 0.3,
                               n_points_pre: int = 500,
                               n_points_post: int = 1500):
    """
    Simulate from 0..t_pre, then set E to E_after while keeping M continuous,
    then simulate t_pre..t_post.
    """
    pA = p.with_A(A_value)

    # Pre-perturb
    t_eval1 = np.linspace(0, t_pre, n_points_pre)
    sol1 = solve_ivp(
        lambda t, y: ec3_ode(t, y, pA),
        t_span=(0, t_pre),
        y0=y0,
        t_eval=t_eval1,
        max_step=0.5,
    )

    E_last, M_last = sol1.y[0, -1], sol1.y[1, -1]
    y_pert = np.array([E_after, M_last])

    # Post-perturb
    t_eval2 = np.linspace(t_pre, t_post, n_points_post)
    sol2 = solve_ivp(
        lambda t, y: ec3_ode(t, y, pA),
        t_span=(t_pre, t_post),
        y0=y_pert,
        t_eval=t_eval2,
        max_step=0.5,
    )

    # Concatenate, avoiding duplicate time at t_pre
    t_combined = np.concatenate([sol1.t, sol2.t[1:]])
    y_combined = np.concatenate([sol1.y, sol2.y[:, 1:]], axis=1)

    return t_combined, y_combined


def save_timecourse_csv(filename: str, t, E, M, meta: dict = None):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        if meta:
            for k, v in meta.items():
                writer.writerow([f"# {k}", v])
        writer.writerow(["t", "E", "M"])
        for ti, Ei, Mi in zip(t, E, M):
            writer.writerow([ti, Ei, Mi])


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Create run dirs for stages we use here
    run_dir_stage2 = make_stage_run_dir("Stage_2_Bifurcation", timestamp)
    run_dir_stage3 = make_stage_run_dir("Stage_3_PhasePlanes", timestamp)
    run_dir_stage4 = make_stage_run_dir("Stage_4_TimeCourses", timestamp)

    p_base = PDParamsEC3()

    # --------------------------------------------------------
    # STAGE 2: BIFURCATION SCAN
    # --------------------------------------------------------
    data_rows, summary = scan_equilibria_vs_A(
        p_base, A_min=0.2, A_max=1.4, nA=61, ngrid=7
    )

    save_equilibria_csv(
        data_rows,
        os.path.join(run_dir_stage2, "ec3_equilibria_vs_A.csv"),
    )
    save_summary_txt(
        summary,
        os.path.join(run_dir_stage2, "ec3_bifurcation_summary.txt"),
    )
    plot_bifurcation_E_vs_A(
        data_rows,
        os.path.join(run_dir_stage2, "ec3_bifurcation_E_vs_A.png"),
    )

    write_readme(
        run_dir_stage2,
        "Stage 2: EC3 Bifurcation Scan",
        [
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- ec3_equilibria_vs_A.csv",
            "- ec3_bifurcation_summary.txt",
            "- ec3_bifurcation_E_vs_A.png",
            "",
            "Description:",
            "Bifurcation diagram E* vs A with stability coding.",
        ],
    )

    # --------------------------------------------------------
    # STAGE 3: PHASE PLANES (SNc and VTA)
    # --------------------------------------------------------
    A_VTA = 0.40
    A_SNc = 1.00

    make_phase_plane_plot(
        p_base,
        A_VTA,
        os.path.join(run_dir_stage3, "phase_plane_VTA_A_0.40.png"),
        title_suffix="VTA-like (low load, monostable)",
    )

    make_phase_plane_plot(
        p_base,
        A_SNc,
        os.path.join(run_dir_stage3, "phase_plane_SNc_A_1.00.png"),
        title_suffix="SNc-like (high load, bistable)",
    )

    write_readme(
        run_dir_stage3,
        "Stage 3: Phase Planes for VTA and SNc",
        [
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- phase_plane_VTA_A_0.40.png",
            "- phase_plane_SNc_A_1.00.png",
            "",
            "Description:",
            "Phase planes (E,M) for low-load VTA-like and high-load SNc-like neurons,",
            "showing vector fields, nullclines, and equilibria.",
        ],
    )

    # --------------------------------------------------------
    # STAGE 4: TIME COURSES & PERTURBATIONS
    # --------------------------------------------------------
    # Baseline timecourses (no perturbation)
    y0_common = [0.9, 0.9] # Start near high-energy state

    t_vta, y_vta = simulate_timecourse(p_base, A_VTA, y0_common)
    t_snc, y_snc = simulate_timecourse(p_base, A_SNc, y0_common)

    plt.figure(figsize=(7, 4))
    # Baseline curves
    plt.plot(t_vta, y_vta[0, :], label="VTA baseline (A=0.40)", linewidth=2)
    plt.plot(t_snc, y_snc[0, :], label="SNc baseline (A=1.00)", linewidth=2)
    
    # Compute (or manually set) equilibria for dashed reference lines
    # If you already know them: 
    # E_vta_eq = 0.767   # example
    # E_snc_eq = 0.565   # example
    E_vta_eq = y_vta[0, -1]
    E_snc_eq = y_snc[0, -1]
    # Faint dashed horizontal lines at equilibria
    plt.axhline(E_vta_eq, color="blue", linestyle="--", alpha=0.3)
    plt.axhline(E_snc_eq, color="orange", linestyle="--", alpha=0.3)

    plt.xlabel("t (a.u.)")
    plt.ylabel("E(t)")
    plt.title("Baseline Time Courses: VTA vs SNc")
    plt.legend()
    plt.tight_layout()
    plt.ylim(0.5, 0.8)
    plt.savefig(os.path.join(run_dir_stage4, "timecourses_VTA_vs_SNc_baseline.png"),
                dpi=300)
    plt.close()

    save_timecourse_csv(
        os.path.join(run_dir_stage4, "timecourse_VTA_baseline.csv"),
        t_vta, y_vta[0, :], y_vta[1, :],
        meta={"A": A_VTA, "scenario": "VTA baseline"},
    )
    save_timecourse_csv(
        os.path.join(run_dir_stage4, "timecourse_SNc_baseline.csv"),
        t_snc, y_snc[0, :], y_snc[1, :],
        meta={"A": A_SNc, "scenario": "SNc baseline"},
    )

    # Perturbation tests
    # Start perturbation from near-steady-state upper branch
    y_snc_ss = y_snc[:, -1]
    y_vta_ss = y_vta[:, -1]

    # SNc collapse: drop E below saddle (e.g., to 0.3)
    t_snc_pert, y_snc_pert = simulate_with_perturbation(
        p_base, A_SNc, y_snc_ss, t_pre=50.0, t_post=300.0, E_after=0.3
    )

    # VTA recovery: same perturbation, but system should recover
    t_vta_pert, y_vta_pert = simulate_with_perturbation(
        p_base, A_VTA, y_vta_ss, t_pre=50.0, t_post=300.0, E_after=0.3
    )

    # Plot SNc collapse vs VTA recovery in one panel
    plt.figure(figsize=(7, 4))
    plt.plot(t_vta_pert, y_vta_pert[0, :], label="VTA (A=0.40) with perturbation")
    plt.plot(t_snc_pert, y_snc_pert[0, :], label="SNc (A=1.00) with perturbation")
    plt.axvline(50.0, color="k", linestyle="--", linewidth=0.8, label="perturbation")
    plt.xlabel("t (a.u.)")
    plt.ylabel("E(t)")
    plt.title("Perturbation: SNc Collapse vs VTA Recovery")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir_stage4, "timecourses_SNc_vs_VTA_perturbation.png"),
                dpi=300)
    plt.close()

    save_timecourse_csv(
        os.path.join(run_dir_stage4, "timecourse_SNc_perturbation.csv"),
        t_snc_pert, y_snc_pert[0, :], y_snc_pert[1, :],
        meta={"A": A_SNc, "scenario": "SNc perturbation (E->0.3 at t=50)"},
    )
    save_timecourse_csv(
        os.path.join(run_dir_stage4, "timecourse_VTA_perturbation.csv"),
        t_vta_pert, y_vta_pert[0, :], y_vta_pert[1, :],
        meta={"A": A_VTA, "scenario": "VTA perturbation (E->0.3 at t=50)"},
    )

    write_readme(
        run_dir_stage4,
        "Stage 4: Time Courses and Perturbation Experiments",
        [
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- timecourses_VTA_vs_SNc_baseline.png",
            "- timecourse_VTA_baseline.csv",
            "- timecourse_SNc_baseline.csv",
            "- timecourses_SNc_vs_VTA_perturbation.png",
            "- timecourse_SNc_perturbation.csv",
            "- timecourse_VTA_perturbation.csv",
            "",
            "Description:",
            "Baseline E(t) trajectories for VTA-like (A=0.40) and SNc-like (A=1.00) neurons,",
            "plus perturbation experiments where E is transiently driven to 0.3 at t=50.",
            "SNc collapses to low-energy attractor; VTA recovers to high-energy state.",
        ],
    )

    print("Done.")
    print("Run directories:")
    print("  Stage_2_Bifurcation:", run_dir_stage2)
    print("  Stage_3_PhasePlanes:", run_dir_stage3)
    print("  Stage_4_TimeCourses:", run_dir_stage4)
