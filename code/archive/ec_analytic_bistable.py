import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# EC2 MODEL: ANALYTIC STEADY-STATES VIA CUBIC IN M
# ============================================================

"""
EC2 dynamics:
    dE/dt = kE (Emax − E) − A M^2 C − gamma (1 − M) E
    dM/dt = kM (1 − M) − beta A C M (1 − E)

At steady state, from dM/dt = 0 we get:
    E(M) = 1 + kM/(A C beta) − kM/(A C beta M)

Plugging this into dE/dt = 0 yields a cubic in M:

    P(M) = a3 M^3 + a2 M^2 + a1 M + a0 = 0

with coefficients (derived symbolically):

    a3 = -A**2 * C**2 * beta
    a2 = A * C * beta * gamma + gamma * kM
    a1 = (A * C * Emax * beta * kE
          - A * C * beta * gamma
          - A * C * beta * kE
          - 2 * gamma * kM
          - kE * kM)
    a0 = gamma * kM + kE * kM

Bistability corresponds to this cubic having three real roots in (0, 1),
with at least two of them stable fixed points of the full 2D system.
"""


# ============================================================
# ANALYTIC UTILITIES
# ============================================================

def cubic_coeffs_M(A, C, beta, gamma, kE=1.0, Emax=1.0, kM=0.5):
    """Return coefficients (a3, a2, a1, a0) of P(M)=0."""
    a3 = -A**2 * C**2 * beta
    a2 = A * C * beta * gamma + gamma * kM
    a1 = (A * C * Emax * beta * kE
          - A * C * beta * gamma
          - A * C * beta * kE
          - 2 * gamma * kM
          - kE * kM)
    a0 = gamma * kM + kE * kM
    return np.array([a3, a2, a1, a0], dtype=float)


def E_from_M(M, A, C, beta, kM=0.5):
    """
    Steady-state E(M) from dM/dt = 0:

        E = 1 + kM/(A C beta) − kM/(A C beta M)
    """
    return 1.0 + kM/(A * C * beta) - kM/(A * C * beta * M)


def jacobian_entries(E, M, A, C, beta, gamma, kE=1.0, kM=0.5):
    """
    Jacobian of the EC2 system at (E, M):

        f(E,M) = dE/dt
        g(E,M) = dM/dt

    f_E = ∂f/∂E = -kE - gamma (1 - M)
    f_M = ∂f/∂M = -2 A C M + gamma E
    g_E = ∂g/∂E = beta A C M
    g_M = ∂g/∂M = -kM - beta A C (1 - E)
    """
    f_E = -kE - gamma * (1.0 - M)
    f_M = -2.0 * A * C * M + gamma * E
    g_E = beta * A * C * M
    g_M = -kM - beta * A * C * (1.0 - E)
    return f_E, f_M, g_E, g_M


def is_stable(E, M, A, C, beta, gamma, kE=1.0, kM=0.5):
    """
    Classify stability via Jacobian trace/determinant.
    Stable if tr < 0 and det > 0.
    """
    f_E, f_M, g_E, g_M = jacobian_entries(E, M, A, C, beta, gamma, kE=kE, kM=kM)
    tr = f_E + g_M
    det = f_E * g_M - f_M * g_E
    return (tr < 0.0) and (det > 0.0)


def fixed_points_for_params(A, C, beta, gamma,
                            kE=1.0, Emax=1.0, kM=0.5,
                            tol_imag=1e-6):
    """
    For given parameters, solve the cubic P(M)=0 analytically (via np.roots),
    then compute the corresponding E(M), and keep only physical fixed points:

        0 < M < 1 and 0 < E < 1

    Returns:
        fixed_points: list of dicts with keys:
            'E', 'M', 'stable' (bool)
    """
    coeffs = cubic_coeffs_M(A, C, beta, gamma, kE=kE, Emax=Emax, kM=kM)
    roots = np.roots(coeffs)

    fixed_points = []
    for r in roots:
        if abs(r.imag) < tol_imag:
            M = r.real
            if 0.0 < M < 1.0:
                E = E_from_M(M, A, C, beta, kM=kM)
                if 0.0 < E < 1.0:
                    stable = is_stable(E, M, A, C, beta, gamma, kE=kE, kM=kM)
                    fixed_points.append({"E": E, "M": M, "stable": stable})
    return fixed_points


def classify_bistability(A, C, beta, gamma,
                         kE=1.0, Emax=1.0, kM=0.5):
    """
    Return:
        n_fixed (int)   - total physical fixed points in (0,1)^2
        n_stable (int)  - number of stable fixed points
        dE_stable (float) - separation in E between stable fixed points
                            (0 if fewer than 2 stables)
        fixed_points (list) - detailed fixed point list
    """
    fps = fixed_points_for_params(A, C, beta, gamma,
                                  kE=kE, Emax=Emax, kM=kM)
    n_fixed = len(fps)
    stable_Es = [fp["E"] for fp in fps if fp["stable"]]
    n_stable = len(stable_Es)
    if n_stable >= 2:
        dE = max(stable_Es) - min(stable_Es)
    else:
        dE = 0.0
    return n_fixed, n_stable, dE, fps


# ============================================================
# PARAMETER SWEEP USING ANALYTIC FIXED POINTS (FAST)
# ============================================================

if __name__ == "__main__":
    # Core biological parameters (same as before)
    A_value = 0.55
    C_value = 1.0
    kE = 1.0
    Emax = 1.0
    kM = 0.5

    # Parameter grid for gamma and beta
    gamma_min, gamma_max, gamma_points = 1.0, 4.0, 40
    beta_min, beta_max, beta_points = 1.0, 3.0, 40

    gamma_vals = np.linspace(gamma_min, gamma_max, gamma_points)
    beta_vals = np.linspace(beta_min, beta_max, beta_points)

    # Arrays to store results
    n_fixed_map   = np.zeros((gamma_points, beta_points))
    n_stable_map  = np.zeros((gamma_points, beta_points))
    dE_stable_map = np.zeros((gamma_points, beta_points))

    best_dE = 0.0
    best_idx = None  # (i, j)

    print("Starting analytic fixed-point scan (no ODE integration)...")
    for i, gamma in enumerate(gamma_vals):
        for j, beta in enumerate(beta_vals):
            n_fixed, n_stable, dE, fps = classify_bistability(
                A=A_value, C=C_value, beta=beta, gamma=gamma,
                kE=kE, Emax=Emax, kM=kM
            )

            n_fixed_map[i, j]   = n_fixed
            n_stable_map[i, j]  = n_stable
            dE_stable_map[i, j] = dE

            if dE > best_dE:
                best_dE = dE
                best_idx = (i, j)

        print(f"gamma={gamma:.3f}: best ΔE so far = {best_dE:.4f}")

    print("Scan complete.")
    if best_idx is not None:
        best_gamma = gamma_vals[best_idx[0]]
        best_beta  = beta_vals[best_idx[1]]
        print(f"Max separation between stable fixed points: ΔE = {best_dE:.4f}")
        print(f"at gamma = {best_gamma:.3f}, beta = {best_beta:.3f}")
    else:
        print("No multi-stable region with ≥2 stable fixed points found on this grid.")

    # --------------------------------------------------------
    # SAVE MAPS TO DISK FOR LATER USE
    # --------------------------------------------------------
    np.savez(
        "ec2_analytic_bistability_maps.npz",
        gamma_vals=gamma_vals,
        beta_vals=beta_vals,
        n_fixed_map=n_fixed_map,
        n_stable_map=n_stable_map,
        dE_stable_map=dE_stable_map,
    )
    print("Saved analytic maps to ec2_analytic_bistability_maps.npz")

    # --------------------------------------------------------
    # PLOT HEATMAPS
    # --------------------------------------------------------
    extent = [beta_min, beta_max, gamma_min, gamma_max]

    # ΔE between stable fixed points
    plt.figure(figsize=(8, 6))
    im = plt.imshow(dE_stable_map,
                    origin="lower",
                    extent=extent,
                    aspect="auto")
    plt.colorbar(im, label="ΔE between stable fixed points")
    plt.xlabel("beta")
    plt.ylabel("gamma")
    plt.title(f"EC2: ΔE between stable fixed points (A={A_value}, C={C_value})")
    plt.tight_layout()
    plt.savefig("ec2_analytic_deltaE_heatmap.png", dpi=300)
    plt.close()
    print("Saved ec2_analytic_deltaE_heatmap.png")

    # Number of stable fixed points (0,1,2,…)
    plt.figure(figsize=(8, 6))
    im2 = plt.imshow(n_stable_map,
                     origin="lower",
                     extent=extent,
                     aspect="auto")
    plt.colorbar(im2, label="# stable fixed points")
    plt.xlabel("beta")
    plt.ylabel("gamma")
    plt.title(f"EC2: number of stable fixed points (A={A_value}, C={C_value})")
    plt.tight_layout()
    plt.savefig("ec2_analytic_nstable_heatmap.png", dpi=300)
    plt.close()
    print("Saved ec2_analytic_nstable_heatmap.png")
