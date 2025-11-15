import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple, Dict
from scipy.integrate import solve_ivp


# ============================================
# MODEL V2: PARAMETERS + ODE
# ============================================

@dataclass
class PDModelParamsV2:
    # neuron-type parameters
    A: float  # axonal arbor factor
    C: float  # Ca2+ load factor

    # energy & mitochondrial parameters
    E_max: float = 1.0
    k_E: float   = 1.0      # rate E relaxes to its setpoint
    k_reg: float = 0.2      # mitochondrial repair / biogenesis

    d0: float    = 0.01     # baseline mito damage
    d_L: float   = 0.02     # load-dependent damage
    d_S: float   = 0.02     # alpha-syn–dependent damage

    # alpha-syn parameters
    beta_basal:  float = 0.0005
    beta_stress: float = 0.02
    k_clear:     float = 0.10

    # thresholds
    E_crit:   float = 0.6   # energetic threshold for mito damage
    E_Scrit:  float = 0.7   # energetic threshold for proteostasis failure
    sigma_E:  float = 0.03
    sigma_S:  float = 0.03


def smooth_heaviside(x: float, sigma: float) -> float:
    """Smooth step function ~0 for x<<0, ~1 for x>>0."""
    return 1.0 / (1.0 + np.exp(-x / sigma))


def pd_ode_v2(t: float, y: np.ndarray, p: PDModelParamsV2) -> np.ndarray:
    """
    y = [E, M, S]
    """
    E, M, S = y
    L = p.A * p.C  # chronic load

    # Energy setpoint: high when M is high and load L is small
    E_base = p.E_max * M / (1.0 + L)

    # 1) Energy relaxes to setpoint
    dE_dt = p.k_E * (E_base - E)

    # 2) Mitochondria: repair vs damage
    H_E = smooth_heaviside(p.E_crit - E, p.sigma_E)   # turns ON when E < E_crit
    dM_dt = (
        p.k_reg * (1.0 - M)                                  # repair
        - (p.d0 + p.d_L * L) * H_E * M                       # energy-gated damage
        - p.d_S * S * M                                      # alpha-syn damage
    )

    # 3) Alpha-syn: basal + stress - clearance
    H_S = smooth_heaviside(p.E_Scrit - E, p.sigma_S)  # turns ON when E < E_Scrit
    dS_dt = (
        p.beta_basal
        + p.beta_stress * H_S
        - p.k_clear * M * S
    )

    return np.array([dE_dt, dM_dt, dS_dt])


# ============================================
# GENERIC SIMULATION HELPERS
# ============================================

def run_simulation(
    params: PDModelParamsV2,
    y0: Tuple[float, float, float] = (1.0, 1.0, 0.0),
    t_span: Tuple[float, float] = (0.0, 500.0),
    t_eval: np.ndarray = None,
    rtol: float = 1e-6,
    atol: float = 1e-8,
) -> Dict[str, np.ndarray]:
    """
    Run a single simulation for given parameters.
    Returns a dict with 't', 'E', 'M', 'S'.
    """
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 2000)

    sol = solve_ivp(
        fun=lambda t, y: pd_ode_v2(t, y, params),
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
    """
    Extract simple summary stats from a trajectory.
    Here we take final values.
    """
    return {
        "E_final": E[-1],
        "M_final": M[-1],
        "S_final": S[-1],
    }


# ============================================
# FIGURE 2: TIME COURSES (RESILIENT VS VULNERABLE)
# ============================================

def example_stress_test():
    base = PDModelParamsV2(A=1.0, C=1.0)

    # Resilient neuron: small arbor + low Ca2+ load
    params_resilient = PDModelParamsV2(
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
        E_crit=base.E_crit,
        E_Scrit=base.E_Scrit,
        sigma_E=base.sigma_E,
        sigma_S=base.sigma_S,
    )

    # Vulnerable SNc-like neuron: large arbor + high Ca2+ load
    params_vulnerable = PDModelParamsV2(
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
        E_crit=base.E_crit,
        E_Scrit=base.E_Scrit,
        sigma_E=base.sigma_E,
        sigma_S=base.sigma_S,
    )

    t_eval = np.linspace(0.0, 500.0, 2000)
    sim_res = run_simulation(params_resilient, t_eval=t_eval)
    sim_vul = run_simulation(params_vulnerable, t_eval=t_eval)
    return sim_res, sim_vul


def plot_figure2_time_courses(output_path: str = "figure2_time_courses_v2.png"):
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


# ============================================
# FIGURE 3: A × C HEATMAPS OF FINAL E AND S
# ============================================

def sweep_A_C(
    A_values: np.ndarray,
    C_values: np.ndarray,
    base_params: PDModelParamsV2,
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
            params = PDModelParamsV2(
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
                E_crit=base_params.E_crit,
                E_Scrit=base_params.E_Scrit,
                sigma_E=base_params.sigma_E,
                sigma_S=base_params.sigma_S,
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


def plot_figure3_heatmaps(output_prefix: str = "figure3_heatmap_v2"):
    base = PDModelParamsV2(A=1.0, C=1.0)

    A_values = np.linspace(0.3, 2.5, 25)
    C_values = np.linspace(0.3, 2.5, 25)

    sweep = sweep_A_C(A_values, C_values, base)

    # Final E heatmap
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

    # Final S heatmap
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


# ============================================
# FIGURE 4: 1D SWEEP OVER A (TIPPING BEHAVIOR)
# ============================================

def sweep_A_1d(
    A_values: np.ndarray,
    base_params: PDModelParamsV2,
    t_span: Tuple[float, float] = (0.0, 500.0),
    t_eval: np.ndarray = None,
) -> Dict[str, np.ndarray]:
    E_final = np.zeros_like(A_values)
    M_final = np.zeros_like(A_values)
    S_final = np.zeros_like(A_values)

    for k, A in enumerate(A_values):
        params = PDModelParamsV2(
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
            E_crit=base_params.E_crit,
            E_Scrit=base_params.E_Scrit,
            sigma_E=base_params.sigma_E,
            sigma_S=base_params.sigma_S,
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


def plot_figure4_tipping(output_path: str = "figure4_tipping_point_v2.png"):
    base = PDModelParamsV2(A=1.0, C=1.0)

    A_values = np.linspace(0.3, 3.0, 40)
    sweep = sweep_A_1d(A_values, base_params=base)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(A_values, sweep["E_final"], label="E_final")
    ax.plot(A_values, sweep["S_final"], label="S_final")
    ax.set_xlabel("Axonal arbor A")
    ax.set_ylabel("Final state")
    ax.set_title("Tipping behavior vs arbor size A")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ============================================
# SUPPLEMENTARY: TRAJECTORIES & SENSITIVITY
# ============================================

def plot_supplement_trajectories(output_path: str = "supplement_trajectories_v2.png"):
    base = PDModelParamsV2(A=1.0, C=1.0)
    t_eval = np.linspace(0.0, 500.0, 2000)

    configs = [
        ("Low load", 0.5, 0.5),
        ("Medium load", 1.0, 1.0),
        ("High load", 2.5, 2.5),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(6, 8), sharex=True)

    for label, A, C in configs:
        params = PDModelParamsV2(
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
            E_crit=base.E_crit,
            E_Scrit=base.E_Scrit,
            sigma_E=base.sigma_E,
            sigma_S=base.sigma_S,
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


def plot_supplement_sensitivity(output_path: str = "supplement_sensitivity_v2.png"):
    """
    Sensitivity to ATP production/energy relaxation: vary k_E
    for a vulnerable neuron (high load).
    """
    base = PDModelParamsV2(A=2.5, C=2.5)
    t_eval = np.linspace(0.0, 500.0, 2000)

    k_E_values = np.linspace(0.5, 2.0, 20)
    E_final = np.zeros_like(k_E_values)
    S_final = np.zeros_like(k_E_values)

    for i, kE in enumerate(k_E_values):
        params = PDModelParamsV2(
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
            E_crit=base.E_crit,
            E_Scrit=base.E_Scrit,
            sigma_E=base.sigma_E,
            sigma_S=base.sigma_S,
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


# ============================================
# MAIN: GENERATE ALL FIGURES
# ============================================

if __name__ == "__main__":
    plot_figure2_time_courses()
    plot_figure3_heatmaps()
    plot_figure4_tipping()
    plot_supplement_trajectories()
    plot_supplement_sensitivity()
    print("Figures (v2) generated.")
