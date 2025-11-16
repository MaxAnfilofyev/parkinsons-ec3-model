# ================================================================
# Supplementary 2D energetic-catastrophe model (EC2)
# ================================================================
# This extends the simple energy-crisis model by adding a mild
# feedback from energetic stress to mitochondrial damage:
#
#   dE/dt = k_E (E_max - E) - A M^2 C
#   dM/dt = k_M (1 - M)     - beta A C (1 - E)
#
# Used ONLY for Supplementary Figures S1–S3.
# ================================================================

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# -----------------------------
# Parameters for EC2 model
# -----------------------------
@dataclass
class PDModelParamsEC2:
    E_max: float = 1.0     # baseline energetic reserve
    k_E: float = 1.0       # energy relaxation rate
    k_M: float = 0.5       # mitochondrial recovery rate
    A: float = 0.55        # axonal arbor size (load)
    C: float = 1.0         # Ca2+ load
    beta: float = 0.6      # coupling: stress-driven mitochondrial damage
    gamma: float = 1.0     # coupling: mitochondrial damage-driven energy loss


# -----------------------------
# ODE system: EC2 model
# -----------------------------
def pd_ode_ec2(t: float, y: np.ndarray, params: PDModelParamsEC2) -> np.ndarray:
    """
    2D energetic-catastrophe model with E<->M coupling.
    State: y = [E, M]
    """
    E, M = y

    
    #dEdt = params.k_E * (params.E_max - E) - params.A * (M ** 2) * params.C
    dEdt = params.k_E*(params.E_max - E) - params.A*M**2*params.C - params.gamma*(1 - M)*E

    #dMdt = params.k_M * (1.0 - M) - params.beta * params.A * params.C * (1.0 - E)
    #dMdt = params.k_M * (1.0 - M) - params.beta * params.A * params.C * M * (1.0 - E)
    #dMdt = params.k_M * (1 - M) - params.beta * params.A * params.C * M * (1 - E)**2
    dMdt = params.k_M * (1 - M) - params.beta * params.A * params.C * M**2 * (1 - E)





    return np.array([dEdt, dMdt])


# -----------------------------
# Simulation helper for EC2
# -----------------------------
def run_simulation_ec2(
    params: PDModelParamsEC2,
    y0=(1.0, 1.0),
    t_span=(0.0, 500.0),
    t_eval=None,
    rtol=1e-6,
    atol=1e-9,
):
    """
    Integrate the EC2 system.
    Returns dict with t, E, M.
    """
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 2000)

    sol = solve_ivp(
        fun=lambda t, y: pd_ode_ec2(t, y, params),
        t_span=t_span,
        y0=np.array(y0, dtype=float),
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )

    E, M = sol.y
    return {"t": sol.t, "E": E, "M": M}


# ================================================================
# Supplementary Figure S1:
# Time courses for resilient vs vulnerable trajectories (EC2)
# ================================================================
def plot_supp_fig_S1_timecourses_ec2(
    A_bistable: float = 0.55,
    C: float = 1.0,
    output_path: str = "supplement_ec2_S1_timecourses.png",
):
    """
    Supplementary Figure S1:
    Time courses of energetic reserve E and mitochondrial capacity M
    for two initial conditions at the same (A, C), illustrating
    resilient vs vulnerable trajectories in the EC2 model.
    """
    params = PDModelParamsEC2(A=A_bistable, C=C)

    t_eval = np.linspace(0.0, 500.0, 2000)

    # High-energy initial condition
    sim_hi = run_simulation_ec2(params, y0=(1.0, 1.0), t_eval=t_eval)
    # Low-energy initial condition (or slightly damaged mitochondria)
    sim_lo = run_simulation_ec2(params, y0=(0.3, 0.8), t_eval=t_eval)

    fig, axes = plt.subplots(2, 1, figsize=(6, 7), sharex=True)

    # Energetic reserve E
    axes[0].plot(sim_hi["t"], sim_hi["E"], label="High-energy IC", color="C0")
    axes[0].plot(sim_lo["t"], sim_lo["E"], label="Low-energy IC", color="C1", linestyle="--")
    axes[0].set_ylabel("Energetic reserve E")
    axes[0].legend(loc="best")

    # Mitochondrial capacity M
    axes[1].plot(sim_hi["t"], sim_hi["M"], label="High-energy IC", color="C0")
    axes[1].plot(sim_lo["t"], sim_lo["M"], label="Low-energy IC", color="C1", linestyle="--")
    axes[1].set_xlabel("Time (arbitrary units)")
    axes[1].set_ylabel("Mitochondrial capacity M")

    fig.suptitle(f"Supplementary Fig. S1 — EC2 time courses at A={A_bistable:.2f}, C={C:.1f}")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ================================================================
# Supplementary Figure S2:
# Bifurcation vs arbor size A (steady-state E*, M*)
# ================================================================
def plot_supp_fig_S2_bifurcation_A_ec2(
    A_values=None,
    C: float = 1.0,
    output_path: str = "supplement_ec2_S2_bifurcation_A.png",
):
    """
    Supplementary Figure S2:
    Sweep over axonal arbor size A, and for each A integrate from
    high- and low-energy initial conditions. Plot steady-state E*
    (and optionally M*) vs A to reveal a bistable window.
    """
    if A_values is None:
        # Tunable: adjust range/spacing as needed
        A_values = np.linspace(0.3, 1.2, 40)

    E_hi_final = []
    E_lo_final = []
    M_hi_final = []
    M_lo_final = []

    t_eval = np.linspace(0.0, 500.0, 2000)

    for A in A_values:
        params = PDModelParamsEC2(A=A, C=C)

        sim_hi = run_simulation_ec2(params, y0=(1.0, 1.0), t_eval=t_eval)
        sim_lo = run_simulation_ec2(params, y0=(0.3, 0.8), t_eval=t_eval)

        E_hi_final.append(sim_hi["E"][-1])
        E_lo_final.append(sim_lo["E"][-1])
        M_hi_final.append(sim_hi["M"][-1])
        M_lo_final.append(sim_lo["M"][-1])

    E_hi_final = np.array(E_hi_final)
    E_lo_final = np.array(E_lo_final)
    M_hi_final = np.array(M_hi_final)
    M_lo_final = np.array(M_lo_final)

    fig, axes = plt.subplots(2, 1, figsize=(6, 7), sharex=True)

    # E* vs A
    axes[0].plot(A_values, E_hi_final, label="High-energy branch", color="C0")
    axes[0].plot(A_values, E_lo_final, label="Low-energy branch", color="C1", linestyle="--")
    axes[0].set_ylabel("Steady-state E*")
    axes[0].legend(loc="best")

    # M* vs A (optional)
    axes[1].plot(A_values, M_hi_final, label="High-energy branch", color="C0")
    axes[1].plot(A_values, M_lo_final, label="Low-energy branch", color="C1", linestyle="--")
    axes[1].set_xlabel("Axonal arbor size A")
    axes[1].set_ylabel("Steady-state M*")

    fig.suptitle("Supplementary Fig. S2 — EC2 bifurcation vs A")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ================================================================
# Supplementary Figure S3:
# Phase plane with proper saddle (EC2)
# ================================================================
def plot_supp_fig_S3_phase_plane_ec2(
    A_bistable: float = 0.55,
    C: float = 1.0,
    output_path: str = "supplement_ec2_S3_phase_plane.png",
):
    """
    Supplementary Figure S3:
    Phase-plane plot for EC2 model at a value of A inside the bistable window.
    Shows:
      - Vector field
      - Analytic nullclines (dE/dt=0, dM/dt=0)
      - Trajectories from high- and low-energy initial conditions.
    """
    params = PDModelParamsEC2(A=A_bistable, C=C)

    # Grid for vector field
    E_min, E_max = 0.0, 1.0
    M_min, M_max = 0.4, 1.5
    n_grid = 40

    E_vals = np.linspace(E_min, E_max, n_grid)
    M_vals = np.linspace(M_min, M_max, n_grid)
    EE, MM = np.meshgrid(E_vals, M_vals)

    dE = np.zeros_like(EE)
    dM = np.zeros_like(MM)

    for i in range(n_grid):
        for j in range(n_grid):
            dE_ij, dM_ij = pd_ode_ec2(0.0, np.array([EE[i, j], MM[i, j]]), params)
            dE[i, j] = dE_ij
            dM[i, j] = dM_ij

    mag = np.sqrt(dE**2 + dM**2) + 1e-12
    u = dE / mag
    v = dM / mag

    # Trajectories
    t_eval = np.linspace(0.0, 500.0, 2000)
    sim_hi = run_simulation_ec2(params, y0=(1.0, 1.0), t_eval=t_eval)
    sim_lo = run_simulation_ec2(params, y0=(0.3, 0.8), t_eval=t_eval)

    # Analytic nullclines
    # dE/dt = 0 -> E = E_max - (A*C/k_E) * M^2
    M_nc = np.linspace(M_min, M_max, 400)
    E_nc = params.E_max - (params.A * params.C / params.k_E) * (M_nc**2)

    # dM/dt = 0 -> M = 1 - (beta*A*C/k_M)*(1 - E)
    E_nc2 = np.linspace(E_min, E_max, 400)
    M_nc2 = 1.0 - (params.beta * params.A * params.C / params.k_M) * (1.0 - E_nc2)

    # Filter to plot window
    mask_E = (E_nc >= E_min) & (E_nc <= E_max)
    E_nc = E_nc[mask_E]
    M_nc = M_nc[mask_E]

    mask_M = (M_nc2 >= M_min) & (M_nc2 <= M_max)
    E_nc2 = E_nc2[mask_M]
    M_nc2 = M_nc2[mask_M]

    # Plot
    fig, ax = plt.subplots(figsize=(6, 5))

    step = 2  # subsample for arrows
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

    # Nullclines
    ax.plot(E_nc, M_nc, color="purple", linewidth=2, label="dE/dt = 0")
    ax.plot(E_nc2, M_nc2, color="purple", linestyle="--", linewidth=2, label="dM/dt = 0")

    # Trajectories
    ax.plot(sim_hi["E"], sim_hi["M"], color="C0", label="High-energy IC")
    ax.plot(sim_lo["E"], sim_lo["M"], color="C1", linestyle="--", label="Low-energy IC")

    # Mark endpoints
    ax.scatter(sim_hi["E"][-1], sim_hi["M"][-1], color="C0")
    ax.scatter(sim_lo["E"][-1], sim_lo["M"][-1], color="C1")

    ax.set_xlim(E_min, E_max)
    ax.set_ylim(M_min, M_max)
    ax.set_xlabel("Energetic reserve E")
    ax.set_ylabel("Mitochondrial capacity M")
    ax.set_title(f"Supplementary Fig. S3 — EC2 phase plane at A={A_bistable:.2f}, C={C:.1f}")
    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


# ================================================================
# Optional: convenience runner to generate all three supplementary
# figures in one call.
# ================================================================
def generate_all_supplementary_ec2_figures(
    A_bistable: float = 0.55,
    C: float = 1.0,
):
    plot_supp_fig_S1_timecourses_ec2(A_bistable=A_bistable, C=C)
    plot_supp_fig_S2_bifurcation_A_ec2(C=C)
    plot_supp_fig_S3_phase_plane_ec2(A_bistable=A_bistable, C=C)


generate_all_supplementary_ec2_figures()