import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple, Dict
from scipy.integrate import solve_ivp


# ======================================================
# ENERGETIC-CRISIS MODEL (V3): PARAMETERS + ODE
# ======================================================

@dataclass
class PDModelParamsEC:
    # neuron-type parameters
    A: float              # axonal arbor factor
    C: float              # Ca2+ load factor

    # energy & mitochondrial parameters
    E_max: float = 1.0
    k_E: float   = 1.0    # rate E relaxes to its setpoint
    k_reg: float = 0.2    # mitochondrial repair / biogenesis

    d0: float  = 0.005    # baseline mito damage
    d_L: float = 0.02     # load-dependent damage
    d_S: float = 0.02     # alpha-syn–dependent damage

    # alpha-syn parameters
    beta_basal:  float = 0.0005
    beta_stress: float = 0.02
    k_clear:     float = 0.08
    K_m:         float = 0.05   # MM constant for clearance
    k_agg:       float = 0.0    # optional S^2 aggregation (off by default)

    # thresholds & nonlinearity for crisis
    E_crit:   float = 0.6       # energetic crisis threshold
    E_Scrit:  float = 0.7       # proteostasis threshold
    sigma_E:  float = 0.02
    sigma_S:  float = 0.03
    n_E:      int   = 6         # steepness of energetic crisis
    k_crisis: float = 0.3       # extra drain when below E_crit


def smooth_heaviside(x: float, sigma: float) -> float:
    """Smooth step; ~0 for x<<0, ~1 for x>>0."""
    return 1.0 / (1.0 + np.exp(-x / sigma))


def pd_ode_ec(t: float, y: np.ndarray, p: PDModelParamsEC) -> np.ndarray:
    """
    Energetic-crisis model.
    y = [E, M, S]
    """
    E, M, S = y
    L = p.A * p.C  # chronic load

    # 1) Energy: relaxes toward setpoint, plus crisis drain when E<E_crit
    E_base = p.E_max * M / (1.0 + L)
    dE_dt  = p.k_E * (E_base - E)

    H_E = smooth_heaviside(p.E_crit - E, p.sigma_E)  # 0 when E >> E_crit
    crisis_term = p.k_crisis * (H_E ** p.n_E) * E    # only strong when E below E_crit
    dE_dt -= crisis_term

    # 2) Mitochondria: repair vs damage, gated by energetic crisis and alpha-syn
    dM_dt = (
        p.k_reg * (1.0 - M)                            # repair to M=1
        - (p.d0 + p.d_L * L) * (H_E ** p.n_E) * M      # crisis damage
        - p.d_S * S * M                                # alpha-syn–mediated damage
    )

    # 3) Alpha-syn: basal + stress - nonlinear clearance (+ optional aggregation)
    H_S = smooth_heaviside(p.E_Scrit - E, p.sigma_S)  # turns on when E below E_Scrit

    # Michaelis-Menten-like clearance, strengthened by M^2
    clearance = p.k_clear * (M ** 2) * S / (p.K_m + S)

    dS_dt = (
        p.beta_basal
        + p.beta_stress * H_S
        + p.k_agg * S * S
        - clearance
    )

    return np.array([dE_dt, dM_dt, dS_dt])


# ======================================================
# GENERIC SIMULATION HELPERS
# ======================================================

def run_simulation(
    params: PDModelParamsEC,
    y0: Tuple[float, float, float] = (1.0, 1.0, 0.0),
    t_span: Tuple[float, float] = (0.0, 500.0),
    t_eval: np.ndarray = None,
    rtol: float = 1e-6,
    atol: float = 1e-8,
) -> Dict[str, np.ndarray]:
    """
    Run a single simulation for given parameters.
    Returns dict with 't', 'E', 'M', 'S'.
    """
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 2000)

    sol = solve_ivp(
        fun=lambda t, y: pd_ode_ec(t, y, params),
        t_span=t_span,
        y0=np.array(y0, dtype=float),
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )

    if not sol.success:
        raise RuntimeError("ODE solver failed: " + sol.message)

    E, M, S = sol.y
    return {"t": sol.t, "E": E, "M": M, "S": S}


def classify_outcome(E: np.ndarray, M: np.ndarray, S: np.ndarray) -> Dict[str, float]:
    """Simple summary: final E, M, S."""
    return {
        "E_final": float(E[-1]),
        "M_final": float(M[-1]),
        "S_final": float(S[-1]),
    }


# ======================================================
# FIGURE 2 — TIME COURSES (RESILIENT VS VULNERABLE)
# ======================================================

def example_stress_test():
    base = PDModelParamsEC(A=1.0, C=1.0)

    # Resilient neuron: low load
    params_resilient = PDModelParamsEC(
        A=0.5,
        C=0.5,
        E_max=base.E_max,
        k_E=base.k_E,
        k_reg=base.k_reg,
        d0=base.d0,
        d_L=base.d_L,
        d_S=base.d_S,
        beta_basal=base.beta_basal,
        beta_stress=base.beta_stress,
        k_clear=base.k_clear,
        K_m=base.K_m,
        k_agg=base.k_agg,
        E_crit=base.E_crit,
        E_Scrit=base.E_Scrit,
        sigma_E=base.sigma_E,
        sigma_S=base.sigma_S,
        n_E=base.n_E,
        k_crisis=base.k_crisis,
    )

    # Vulnerable SNc-like neuron: high load
    params_vulnerable = PDModelParamsEC(
        A=2.5,
        C=2.5,
        E_max=base.E_max,
        k_E=base.k_E,
        k_reg=base.k_reg,
        d0=base.d0,
        d_L=base.d_L,
        d_S=base.d_S,
        beta_basal=base.beta_basal,
        beta_stress=base.beta_stress,
        k_clear=base.k_clear,
        K_m=base.K_m,
        k_agg=base.k_agg,
        E_crit=base.E_crit,
        E_Scrit=base.E_Scrit,
        sigma_E=base.sigma_E,
        sigma_S=base.sigma_S,
        n_E=base.n_E,
        k_crisis=base.k_crisis,
    )

    t_eval = np.linspace(0.0, 500.0, 2000)
    sim_res = run_simulation(params_resilient, t_eval=t_eval)
    sim_vul = run_simulation(params_vulnerable, t_eval=t_eval)
    return sim_res, sim_vul


def plot_figure2_time_courses(output_path: str = "figure2_time_courses_ec.png"):
    sim_res, sim_vul = example_stress_test()

    t_res, E_res, M_res, S_res = sim_res["t"], sim_res["E"], sim_res["M"], sim_res["S"]
    t_vul, E_vul, M_vul, S_vul = sim_vul["t"], sim_vul["E"], sim_vul["M"], sim_vul["S"]

    fig, axes = plt.subplots(3, 1, figsize=(6, 8), sharex=True)

    axes[0].plot(t_res, E_res, label="Resilient")
    axes[0].plot(t_vul, E_vul, linestyle="--", label="Vulnerable")
    axes[0].set_ylabel("Energetic reserve E")
    axes[0].legend(loc="best")

    axes[1].plot(t_res, M_res, label="Resilient")
    axes[1].plot(t_vul, M_vul, linestyle="--", label="Vulnerable")
    axes[1].set_ylabel("Mitochondrial capacity M")

    axes[2].plot(t_res, S_res, label="Resilient")
    axes[2].plot(t_vul, S_vul, linestyle="--", label="Vulnerable")
    axes[2].set_ylabel("α-syn burden S")
    axes[2].set_xlabel("Time (arbitrary units)")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ======================================================
# FIGURE 3 — A×C HEATMAPS (FINAL E AND S)
# ======================================================

def sweep_A_C(
    A_values: np.ndarray,
    C_values: np.ndarray,
    base_params: PDModelParamsEC,
    t_span: Tuple[float, float] = (0.0, 500.0),
    t_eval: np.ndarray = None,
) -> Dict[str, np.ndarray]:
    nA = len(A_values)
    nC = len(C_values)

    E_grid = np.zeros((nA, nC))
    M_grid = np.zeros((nA, nC))
    S_grid = np.zeros((nA, nC))

    for i, A in enumerate(A_values):
        for j, C in enumerate(C_values):
            params = PDModelParamsEC(
                A=A,
                C=C,
                E_max=base_params.E_max,
                k_E=base_params.k_E,
                k_reg=base_params.k_reg,
                d0=base_params.d0,
                d_L=base_params.d_L,
                d_S=base_params.d_S,
                beta_basal=base_params.beta_basal,
                beta_stress=base_params.beta_stress,
                k_clear=base_params.k_clear,
                K_m=base_params.K_m,
                k_agg=base_params.k_agg,
                E_crit=base_params.E_crit,
                E_Scrit=base_params.E_Scrit,
                sigma_E=base_params.sigma_E,
                sigma_S=base_params.sigma_S,
                n_E=base_params.n_E,
                k_crisis=base_params.k_crisis,
            )
            sim = run_simulation(params, t_span=t_span, t_eval=t_eval)
            summary = classify_outcome(sim["E"], sim["M"], sim["S"])
            E_grid[i, j] = summary["E_final"]
            M_grid[i, j] = summary["M_final"]
            S_grid[i, j] = summary["S_final"]

    return {
        "A_values": A_values,
        "C_values": C_values,
        "E_final": E_grid,
        "M_final": M_grid,
        "S_final": S_grid,
    }


def plot_figure3_heatmaps(output_prefix: str = "figure3_heatmap_ec"):
    base = PDModelParamsEC(A=1.0, C=1.0)

    A_values = np.linspace(0.3, 2.5, 25)
    C_values = np.linspace(0.3, 2.5, 25)

    sweep = sweep_A_C(A_values, C_values, base)

    # Final E
    fig1, ax1 = plt.subplots(figsize=(5, 4))
    im1 = ax1.imshow(
        sweep["E_final"],
        origin="lower",
        extent=[C_values[0], C_values[-1], A_values[0], A_values[-1]],
        aspect="auto",
    )
    ax1.set_xlabel("Ca²⁺ load C")
    ax1.set_ylabel("Axonal arbor A")
    ax1.set_title("Final energetic reserve E")
    fig1.colorbar(im1, ax=ax1, label="E_final")
    fig1.tight_layout()
    fig1.savefig(f"{output_prefix}_E.png", dpi=300)
    plt.close(fig1)

    # Final S
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    im2 = ax2.imshow(
        sweep["S_final"],
        origin="lower",
        extent=[C_values[0], C_values[-1], A_values[0], A_values[-1]],
        aspect="auto",
    )
    ax2.set_xlabel("Ca²⁺ load C")
    ax2.set_ylabel("Axonal arbor A")
    ax2.set_title("Final α-syn burden S")
    fig2.colorbar(im2, ax=ax2, label="S_final")
    fig2.tight_layout()
    fig2.savefig(f"{output_prefix}_S.png", dpi=300)
    plt.close(fig2)


# ======================================================
# FIGURE 4 — 1D SWEEP OVER A (TIPPING CURVE)
# ======================================================

def sweep_A_1d(
    A_values: np.ndarray,
    base_params: PDModelParamsEC,
    t_span: Tuple[float, float] = (0.0, 500.0),
    t_eval: np.ndarray = None,
) -> Dict[str, np.ndarray]:
    E_final = np.zeros_like(A_values)
    M_final = np.zeros_like(A_values)
    S_final = np.zeros_like(A_values)

    for k, A in enumerate(A_values):
        params = PDModelParamsEC(
            A=A,
            C=base_params.C,
            E_max=base_params.E_max,
            k_E=base_params.k_E,
            k_reg=base_params.k_reg,
            d0=base_params.d0,
            d_L=base_params.d_L,
            d_S=base_params.d_S,
            beta_basal=base_params.beta_basal,
            beta_stress=base_params.beta_stress,
            k_clear=base_params.k_clear,
            K_m=base_params.K_m,
            k_agg=base_params.k_agg,
            E_crit=base_params.E_crit,
            E_Scrit=base_params.E_Scrit,
            sigma_E=base_params.sigma_E,
            sigma_S=base_params.sigma_S,
            n_E=base_params.n_E,
            k_crisis=base_params.k_crisis,
        )
        sim = run_simulation(params, t_span=t_span, t_eval=t_eval)
        summary = classify_outcome(sim["E"], sim["M"], sim["S"])
        E_final[k] = summary["E_final"]
        M_final[k] = summary["M_final"]
        S_final[k] = summary["S_final"]

    return {
        "A_values": A_values,
        "E_final": E_final,
        "M_final": M_final,
        "S_final": S_final,
    }


def plot_figure4_tipping(output_path: str = "figure4_tipping_point_ec.png"):
    base = PDModelParamsEC(A=1.0, C=1.0)

    A_values = np.linspace(0.3, 3.0, 30)
    sweep = sweep_A_1d(A_values, base_params=base)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(A_values, sweep["E_final"], label="E_final")
    ax.plot(A_values, sweep["S_final"], label="S_final")
    ax.set_xlabel("Axonal arbor A")
    ax.set_ylabel("Final state")
    ax.set_title("Energetic fragility vs arbor size A")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ======================================================
# SUPPLEMENTARY FIGURES
# ======================================================

def plot_supplement_trajectories(output_path: str = "supplement_trajectories_ec.png"):
    base = PDModelParamsEC(A=1.0, C=1.0)
    t_eval = np.linspace(0.0, 500.0, 2000)

    configs = [
        ("Low load", 0.5, 0.5),
        ("Medium load", 1.0, 1.0),
        ("High load", 2.5, 2.5),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(6, 8), sharex=True)

    for label, A, C in configs:
        params = PDModelParamsEC(
            A=A,
            C=C,
            E_max=base.E_max,
            k_E=base.k_E,
            k_reg=base.k_reg,
            d0=base.d0,
            d_L=base.d_L,
            d_S=base.d_S,
            beta_basal=base.beta_basal,
            beta_stress=base.beta_stress,
            k_clear=base.k_clear,
            K_m=base.K_m,
            k_agg=base.k_agg,
            E_crit=base.E_crit,
            E_Scrit=base.E_Scrit,
            sigma_E=base.sigma_E,
            sigma_S=base.sigma_S,
            n_E=base.n_E,
            k_crisis=base.k_crisis,
        )
        sim = run_simulation(params, t_eval=t_eval)
        axes[0].plot(sim["t"], sim["E"], label=label)
        axes[1].plot(sim["t"], sim["M"], label=label)
        axes[2].plot(sim["t"], sim["S"], label=label)

    axes[0].set_ylabel("E")
    axes[1].set_ylabel("M")
    axes[2].set_ylabel("S")
    axes[2].set_xlabel("Time")

    for ax in axes:
        ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_supplement_sensitivity(output_path: str = "supplement_sensitivity_ec.png"):
    """
    Sensitivity to energy relaxation rate k_E for a vulnerable neuron.
    """
    base = PDModelParamsEC(A=2.5, C=2.5)
    t_eval = np.linspace(0.0, 500.0, 2000)

    k_E_values = np.linspace(0.5, 2.0, 20)
    E_final = np.zeros_like(k_E_values)
    S_final = np.zeros_like(k_E_values)

    for i, kE in enumerate(k_E_values):
        params = PDModelParamsEC(
            A=base.A,
            C=base.C,
            E_max=base.E_max,
            k_E=kE,
            k_reg=base.k_reg,
            d0=base.d0,
            d_L=base.d_L,
            d_S=base.d_S,
            beta_basal=base.beta_basal,
            beta_stress=base.beta_stress,
            k_clear=base.k_clear,
            K_m=base.K_m,
            k_agg=base.k_agg,
            E_crit=base.E_crit,
            E_Scrit=base.E_Scrit,
            sigma_E=base.sigma_E,
            sigma_S=base.sigma_S,
            n_E=base.n_E,
            k_crisis=base.k_crisis,
        )
        sim = run_simulation(params, t_eval=t_eval)
        summary = classify_outcome(sim["E"], sim["M"], sim["S"])
        E_final[i] = summary["E_final"]
        S_final[i] = summary["S_final"]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(k_E_values, E_final, label="E_final")
    ax.plot(k_E_values, S_final, label="S_final")
    ax.set_xlabel("k_E (energy relaxation rate)")
    ax.set_ylabel("Final state")
    ax.set_title("Sensitivity to energy dynamics (vulnerable neuron)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ======================================================
# BIFURCATION ANALYSIS — BRANCHES VS ARBOR SIZE A
# ======================================================

def bifurcation_branches_A(
    A_values: np.ndarray,
    C: float = 1.0,
    t_end: float = 800.0,
) -> Dict[str, np.ndarray]:
    """
    For each A, integrate from two widely separated initial conditions:

      - 'high' IC: near healthy state (E=1, M=1, S=0)
      - 'low'  IC: pre-collapsed (E=0.1, M=0.6, S=0.1)

    Returns arrays for E* and S* on each branch.
    Where the two branches differ, the system is bistable.
    """
    E_hi, S_hi = [], []
    E_lo, S_lo = [], []

    for A in A_values:
        # High-energy initial condition
        p_hi = PDModelParamsEC(A=A, C=C)
        sim_hi = run_simulation(p_hi, y0=(1.0, 1.0, 0.0), t_span=(0.0, t_end))
        summ_hi = classify_outcome(sim_hi["E"], sim_hi["M"], sim_hi["S"])
        E_hi.append(summ_hi["E_final"])
        S_hi.append(summ_hi["S_final"])

        # Low-energy initial condition
        p_lo = PDModelParamsEC(A=A, C=C)
        sim_lo = run_simulation(p_lo, y0=(0.1, 0.6, 0.1), t_span=(0.0, t_end))
        summ_lo = classify_outcome(sim_lo["E"], sim_lo["M"], sim_lo["S"])
        E_lo.append(summ_lo["E_final"])
        S_lo.append(summ_lo["S_final"])

    return {
        "A_values": A_values,
        "E_hi": np.array(E_hi),
        "S_hi": np.array(S_hi),
        "E_lo": np.array(E_lo),
        "S_lo": np.array(S_lo),
    }
def plot_figure5_bifurcation_A(
    output_path: str = "figure5_bifurcation_A_ec.png",
    C: float = 1.0,
):
    """
    Bifurcation-style plot:

      Top panel: E* vs A for two branches (high-IC, low-IC)
      Bottom panel: S* vs A for two branches

    The region where |E_hi - E_lo| > eps is the bistable regime.
    """
    # Choose a range of A that spans healthy → vulnerable
    A_values = np.linspace(0.4, 1.0, 40)
    branches = bifurcation_branches_A(A_values, C=C)

    A_vals = branches["A_values"]
    E_hi, S_hi = branches["E_hi"], branches["S_hi"]
    E_lo, S_lo = branches["E_lo"], branches["S_lo"]

    # Identify bistable window
    eps = 1e-3
    diff = np.abs(E_hi - E_lo)
    bistable_mask = diff > eps

    fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

    # --- Panel 1: Energetic reserve E* ---
    ax = axes[0]
    ax.plot(A_vals, E_hi, label="High-E branch (high IC)")
    ax.plot(A_vals, E_lo, linestyle="--", label="Low-E branch (low IC)")
    ax.set_ylabel("E* (steady-state energetic reserve)")
    ax.set_title("Bifurcation vs arbor size A")

    # Shade bistable region if present
    if np.any(bistable_mask):
        A_bi = A_vals[bistable_mask]
        ax.axvspan(A_bi[0], A_bi[-1], alpha=0.15)
        # Optionally annotate
        ax.text(
            0.5 * (A_bi[0] + A_bi[-1]),
            0.5 * (np.max(E_hi) + np.min(E_lo)),
            "bistable\nregion",
            ha="center",
            va="center",
        )

    ax.legend(loc="best")

    # --- Panel 2: α-syn burden S* ---
    ax2 = axes[1]
    ax2.plot(A_vals, S_hi, label="High-E branch (high IC)")
    ax2.plot(A_vals, S_lo, linestyle="--", label="Low-E branch (low IC)")
    ax2.set_xlabel("Axonal arbor A")
    ax2.set_ylabel("S* (steady-state α-syn)")

    if np.any(bistable_mask):
        A_bi = A_vals[bistable_mask]
        ax2.axvspan(A_bi[0], A_bi[-1], alpha=0.15)

    ax2.legend(loc="best")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

# ======================================================
# STEADY STATE HELPER
# ======================================================

def steady_state_for_A(
    A: float,
    C: float = 1.0,
    y0: Tuple[float, float, float] = (1.0, 1.0, 0.0),
    t_end: float = 1000.0,
) -> Dict[str, float]:
    """
    Integrate the full 3D system for given A, C and return approximate steady state.
    """
    params = PDModelParamsEC(A=A, C=C)
    sim = run_simulation(params, y0=y0, t_span=(0.0, t_end))
    return classify_outcome(sim["E"], sim["M"], sim["S"])

# ======================================================
# PHASE PLANE / NULLCLINES AT FIXED S*
# ======================================================

def compute_phase_plane_data(
    A: float,
    C: float = 1.0,
    E_range: Tuple[float, float] = (0.0, 1.0),
    M_range: Tuple[float, float] = (0.5, 1.0),
    n_grid: int = 100,
    t_steady: float = 1000.0,
) -> Dict[str, np.ndarray]:
    """
    Compute dE/dt and dM/dt on an (E, M) grid, holding S fixed at the
    steady-state value S* for a high-energy initial condition.

    This effectively gives a 2D slice of the 3D dynamics near the high-E branch.
    """
    # Get S* at high-energy IC
    ss_high = steady_state_for_A(A, C, y0=(1.0, 1.0, 0.0), t_end=t_steady)
    S_star = ss_high["S_final"]

    params = PDModelParamsEC(A=A, C=C)

    E_vals = np.linspace(E_range[0], E_range[1], n_grid)
    M_vals = np.linspace(M_range[0], M_range[1], n_grid)
    EE, MM = np.meshgrid(E_vals, M_vals)

    dE = np.zeros_like(EE)
    dM = np.zeros_like(EE)

    for i in range(n_grid):
        for j in range(n_grid):
            y = np.array([EE[i, j], MM[i, j], S_star])
            d = pd_ode_ec(0.0, y, params)
            dE[i, j] = d[0]
            dM[i, j] = d[1]

    return {
        "A": A,
        "C": C,
        "S_star": S_star,
        "E_vals": E_vals,
        "M_vals": M_vals,
        "EE": EE,
        "MM": MM,
        "dE": dE,
        "dM": dM,
    }

def plot_figure6_phase_plane_bistable(
    A_bistable: float = 0.55,
    C: float = 1.0,
    output_path: str = "figure6_phase_plane_bistable_ec.png",
):
    """
    Phase-plane plot in (E, M) for a value of A inside the bistable window.

    - dE/dt = 0 and dM/dt = 0 nullclines
    - Vector field (dE/dt, dM/dt)
    - Two trajectories: high-energy IC and low-energy IC
    """
    # ---------- phase-plane data ----------
    plane = compute_phase_plane_data(
        A=A_bistable,
        C=C,
        E_range=(0.0, 1.0),
        M_range=(0.55, 1.5),
        n_grid=60,
    )

    EE, MM = plane["EE"], plane["MM"]
    dE, dM = plane["dE"], plane["dM"]

    # Normalize vector field for quiver
    mag = np.sqrt(dE**2 + dM**2) + 1e-9
    u = dE / mag
    v = dM / mag

    # ---------- trajectories ----------
    params = PDModelParamsEC(A=A_bistable, C=C)
    t_eval = np.linspace(0.0, 500.0, 2000)

    # High-energy trajectory
    sim_hi = run_simulation(params, y0=(1.0, 1.0, 0.0), t_eval=t_eval)
    # Low-energy trajectory
    sim_lo = run_simulation(params, y0=(0.1, 0.6, 0.1), t_eval=t_eval)

    # ---------- plot ----------
    fig, ax = plt.subplots(figsize=(6, 5))

    # Vector field (down-sampled so it doesn't clutter)
    step = 4
    ax.quiver(
        EE[::step, ::step],
        MM[::step, ::step],
        u[::step, ::step],
        v[::step, ::step],
        angles="xy",
        scale_units="xy",
        scale=20.0,
        alpha=0.5,
    )

    # Nullclines: dE/dt = 0 and dM/dt = 0
    ax.contour(EE, MM, dE, levels=[0.0], linewidths=2.0, linestyles="solid")
    ax.contour(EE, MM, dM, levels=[0.0], linewidths=2.0, linestyles="dashed")

    # Trajectories in (E, M)
    ax.plot(sim_hi["E"], sim_hi["M"], label="High-energy IC")
    ax.plot(sim_lo["E"], sim_lo["M"], label="Low-energy IC")

    # Mark approximate steady states (final points)
    ax.scatter(sim_hi["E"][-1], sim_hi["M"][-1])
    ax.scatter(sim_lo["E"][-1], sim_lo["M"][-1])

    ax.set_xlabel("Energetic reserve E")
    ax.set_ylabel("Mitochondrial capacity M")
    ax.set_title(f"Phase plane at A = {A_bistable:.2f}, C = {C:.1f}")
    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

def plot_figure6_phase_plane_bistable_analytic(
    A_bistable: float = 0.55,
    C: float = 1.0,
    output_path: str = "figure6_phase_plane_bistable_ec.png",
):
    """
    Phase-plane plot in (E, M) for a value of A inside the bistable window,
    using *analytic* nullclines:

        dE/dt = 0 -> E = E_max - (A*C/k_E) * M^2
        dM/dt = 0 -> M = M_max - (alpha*A*C/k_M) * E

    - Gray vector field computed from the full ODE with S fixed at S*
      (as in compute_phase_plane_data).
    - Blue/orange trajectories from high vs low initial conditions.
    - Smooth purple nullclines overlaid.

    This removes numerical artifacts from contour-based nullclines.
    """
    # ----------------- parameters -----------------
    params = PDModelParamsEC(A=A_bistable, C=C)

    # ----------------- vector field -----------------
    plane = compute_phase_plane_data(
        A=A_bistable,
        C=C,
        E_range=(0.0, 1.0),
        M_range=(0.55, 1.5),
        n_grid=60,
    )
    EE, MM = plane["EE"], plane["MM"]
    dE, dM = plane["dE"], plane["dM"]

    mag = np.sqrt(dE**2 + dM**2) + 1e-9
    u = dE / mag
    v = dM / mag

    # ----------------- trajectories -----------------
    t_eval = np.linspace(0.0, 500.0, 2000)

    sim_hi = run_simulation(params, y0=(1.0, 1.0, 0.0), t_eval=t_eval)
    sim_lo = run_simulation(params, y0=(0.1, 0.6, 0.1), t_eval=t_eval)

    # ----------------- analytic nullclines -----------------
    E_min, E_max_plot = EE.min(), EE.max()
    M_min, M_max_plot = MM.min(), MM.max()

    # E-nullcline: E = E_max - (A*C/k_E) * M^2
    M_nc = np.linspace(M_min, M_max_plot, 400)
    E_nc = (
        params.E_max
        - (A_bistable * C / params.k_E) * (M_nc**2)
    )

    # M-nullcline: M = 1.0 (horizontal line)
    E_nc2 = np.linspace(E_min, E_max_plot, 400)
    M_nc2 = np.ones_like(E_nc2)

    # Filter both to plotting window
    mask_E = (E_nc >= E_min) & (E_nc <= E_max_plot)
    E_nc = E_nc[mask_E]
    M_nc = M_nc[mask_E]

    mask_M = (M_nc2 >= M_min) & (M_nc2 <= M_max_plot)
    E_nc2 = E_nc2[mask_M]
    M_nc2 = M_nc2[mask_M]

    # ----------------- plot -----------------
    fig, ax = plt.subplots(figsize=(6, 5))

    # Vector field
    step = 3
    ax.quiver(
        EE[::step, ::step],
        MM[::step, ::step],
        u[::step, ::step],
        v[::step, ::step],
        angles="xy",
        scale_units="xy",
        scale=20.0,
        alpha=0.4,
        color="gray",
    )

    # Analytic nullclines
    ax.plot(E_nc, M_nc, color="purple", linewidth=2, label="dE/dt = 0")
    ax.plot(E_nc2, M_nc2, color="purple", linestyle="--", linewidth=2, label="dM/dt = 0")

    # Trajectories
    ax.plot(sim_hi["E"], sim_hi["M"], label="High-energy IC", color="C0")
    ax.plot(sim_lo["E"], sim_lo["M"], label="Low-energy IC", color="C1")

    # Mark final steady states
    ax.scatter(sim_hi["E"][-1], sim_hi["M"][-1], color="C0")
    ax.scatter(sim_lo["E"][-1], sim_lo["M"][-1], color="C1")

    ax.set_xlabel("Energetic reserve E")
    ax.set_ylabel("Mitochondrial capacity M")
    ax.set_title(f"Phase plane at A = {A_bistable:.2f}, C = {C:.1f}")
    ax.set_xlim(E_min, E_max_plot)
    ax.set_ylim(M_min, M_max_plot)

    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ======================================================
# MAIN
# ======================================================

if __name__ == "__main__":
    #plot_figure2_time_courses()
    #plot_figure3_heatmaps()
    #plot_figure4_tipping()
    #plot_supplement_trajectories()
    #plot_supplement_sensitivity()
    #plot_figure5_bifurcation_A()
    #plot_figure6_phase_plane_bistable()
    plot_figure6_phase_plane_bistable_analytic()

    print("Energetic-crisis figures generated (including bifurcation + phase plane).")

