import numpy as np
from dataclasses import dataclass
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import csv

# ============================================================
# PARAMETER CLASS FOR EC2 MODEL
# ============================================================

@dataclass
class PDParamsEC2:
    k_E: float = 1.0
    E_max: float = 1.0
    k_M: float = 0.5
    beta: float = 1.0
    A: float = 0.55
    C: float = 1.0
    gamma: float = 1.0


# ============================================================
# EC2 ODE SYSTEM
# ============================================================

def pd_ode_ec2(t, y, p: PDParamsEC2):
    """
    EC2 energetic crisis model with feedback:
        dE/dt = kE (Emax − E) − A M^2 C − gamma (1 − M) E
        dM/dt = kM (1 − M) − beta A C M (1 − E)
    """
    E, M = y

    dE = (
        p.k_E * (p.E_max - E)
        - p.A * (M ** 2) * p.C
        - p.gamma * (1.0 - M) * E
    )

    dM = (
        p.k_M * (1.0 - M)
        - p.beta * p.A * p.C * M * (1.0 - E)
    )

    return np.array([dE, dM])


# ============================================================
# SIMULATION UTILITIES
# ============================================================

def simulate_ec2(params: PDParamsEC2, y0, tmax=500.0, n_points=4000):
    """
    Simulate EC2 from initial conditions y0 = [E0, M0].
    """
    sol = solve_ivp(
        fun=lambda t, y: pd_ode_ec2(t, y, params),
        t_span=(0.0, tmax),
        y0=y0,
        dense_output=True,
        max_step=0.1,
    )
    t = np.linspace(0.0, tmax, n_points)
    y = sol.sol(t)
    return t, y


def steady_state_ec2(params: PDParamsEC2, y0, tmax=1000.0):
    """
    Run a long simulation and return the final state as the steady state.
    """
    t, y = simulate_ec2(params, y0, tmax=tmax, n_points=4000)
    # y is shape (2, n_points)
    E_final = y[0, -1]
    M_final = y[1, -1]
    return np.array([E_final, M_final])


# ============================================================
# BISTABILITY CHECK
# ============================================================

def check_bistability_for_params(params: PDParamsEC2,
                                 y0_high=[1.0, 1.0],
                                 y0_low=[0.2, 0.5],
                                 tmax=1000.0,
                                 deltaE_threshold=0.02):
    """
    For a given parameter set, compute steady states from high and low
    initial conditions and decide if they are distinct (bistable).

    Returns:
        (E_hi, M_hi, E_lo, M_lo, deltaE, is_bistable)
    """
    ss_hi = steady_state_ec2(params, y0=y0_high, tmax=tmax)
    ss_lo = steady_state_ec2(params, y0=y0_low, tmax=tmax)

    E_hi, M_hi = ss_hi
    E_lo, M_lo = ss_lo

    deltaE = abs(E_hi - E_lo)
    is_bistable = deltaE > deltaE_threshold

    return E_hi, M_hi, E_lo, M_lo, deltaE, is_bistable


# ============================================================
# MAIN SCAN + FIGURE GENERATION
# ============================================================

if __name__ == "__main__":

    # --------------------------------------------------------
    # CONFIGURATION
    # --------------------------------------------------------
    # Core load parameters (can be tuned)
    A_value = 0.55
    C_value = 1.0

    # Gamma and beta sweep ranges
    gamma_min, gamma_max, gamma_points = 1.0, 4.0, 16  # inclusive range
    beta_min, beta_max, beta_points = 1.0, 3.0, 11     # inclusive range

    # Simulation and bistability thresholds
    tmax_ss = 1000.0        # time horizon for reaching steady state
    deltaE_threshold = 0.02 # minimum E separation to call it bistable

    # Initial conditions for high/low branches
    y0_high = [1.0, 1.0]
    y0_low  = [0.2, 0.5]

    # Output file names
    csv_filename = "ec2_bistability_scan.csv"
    heatmap_delta_filename = "ec2_deltaE_heatmap.png"
    heatmap_mask_filename  = "ec2_bistable_mask.png"
    example_tc_filename    = "ec2_example_timecourses.png"

    # --------------------------------------------------------
    # PREPARE GRIDS
    # --------------------------------------------------------
    gamma_vals = np.linspace(gamma_min, gamma_max, gamma_points)
    beta_vals  = np.linspace(beta_min,  beta_max,  beta_points)

    # Allocate result arrays
    E_hi_map = np.zeros((gamma_points, beta_points))
    M_hi_map = np.zeros((gamma_points, beta_points))
    E_lo_map = np.zeros((gamma_points, beta_points))
    M_lo_map = np.zeros((gamma_points, beta_points))
    deltaE_map = np.zeros((gamma_points, beta_points))
    bistable_map = np.zeros((gamma_points, beta_points), dtype=bool)

    # Track the most bistable point found
    best_deltaE = -np.inf
    best_indices = (None, None)

    # --------------------------------------------------------
    # SWEEP OVER GAMMA & BETA
    # --------------------------------------------------------
    print("Starting EC2 bistability scan...")
    for i, gamma in enumerate(gamma_vals):
        for j, beta in enumerate(beta_vals):
            # Create params for this point
            params = PDParamsEC2(
                k_E=1.0,
                E_max=1.0,
                k_M=0.5,
                beta=beta,
                A=A_value,
                C=C_value,
                gamma=gamma
            )

            E_hi, M_hi, E_lo, M_lo, deltaE, is_bistable = check_bistability_for_params(
                params,
                y0_high=y0_high,
                y0_low=y0_low,
                tmax=tmax_ss,
                deltaE_threshold=deltaE_threshold
            )

            E_hi_map[i, j] = E_hi
            M_hi_map[i, j] = M_hi
            E_lo_map[i, j] = E_lo
            M_lo_map[i, j] = M_lo
            deltaE_map[i, j] = deltaE
            bistable_map[i, j] = is_bistable

            if deltaE > best_deltaE:
                best_deltaE = deltaE
                best_indices = (i, j)

            print(f"gamma={gamma:.3f}, beta={beta:.3f} "
                  f"=> E_hi={E_hi:.3f}, E_lo={E_lo:.3f}, ΔE={deltaE:.4f}, "
                  f"bistable={is_bistable}")

    print("Scan complete.")
    print(f"Largest ΔE found: {best_deltaE:.4f} "
          f"at gamma={gamma_vals[best_indices[0]]:.3f}, "
          f"beta={beta_vals[best_indices[1]]:.3f}")

    # --------------------------------------------------------
    # SAVE RESULTS TO CSV
    # --------------------------------------------------------
    print(f"Saving CSV results to {csv_filename} ...")
    with open(csv_filename, mode="w", newline="") as f:
        writer = csv.writer(f)
        header = [
            "gamma", "beta",
            "E_hi", "M_hi",
            "E_lo", "M_lo",
            "deltaE",
            "bistable"
        ]
        writer.writerow(header)
        for i, gamma in enumerate(gamma_vals):
            for j, beta in enumerate(beta_vals):
                writer.writerow([
                    float(gamma),
                    float(beta),
                    float(E_hi_map[i, j]),
                    float(M_hi_map[i, j]),
                    float(E_lo_map[i, j]),
                    float(M_lo_map[i, j]),
                    float(deltaE_map[i, j]),
                    int(bistable_map[i, j])
                ])
    print("CSV saved.")

    # --------------------------------------------------------
    # PLOT HEATMAP OF ΔE
    # --------------------------------------------------------
    print(f"Saving ΔE heatmap to {heatmap_delta_filename} ...")
    plt.figure(figsize=(8, 6))
    # imshow expects [row, col] => gamma is y-axis, beta is x-axis
    extent = [beta_min, beta_max, gamma_min, gamma_max]
    im = plt.imshow(deltaE_map,
                    origin="lower",
                    extent=extent,
                    aspect="auto")
    plt.colorbar(im, label="ΔE = |E_hi − E_lo|")
    plt.xlabel("beta")
    plt.ylabel("gamma")
    plt.title(f"EC2 ΔE Heatmap (A={A_value}, C={C_value})")
    plt.tight_layout()
    plt.savefig(heatmap_delta_filename, dpi=300)
    plt.close()
    print("ΔE heatmap saved.")

    # --------------------------------------------------------
    # PLOT BISTABLE MASK (0/1)
    # --------------------------------------------------------
    print(f"Saving bistable mask heatmap to {heatmap_mask_filename} ...")
    plt.figure(figsize=(8, 6))
    im2 = plt.imshow(bistable_map.astype(float),
                     origin="lower",
                     extent=extent,
                     aspect="auto")
    plt.colorbar(im2, label="Bistable (1) / Not (0)")
    plt.xlabel("beta")
    plt.ylabel("gamma")
    plt.title(f"EC2 Bistability Mask (ΔE > {deltaE_threshold}, A={A_value}, C={C_value})")
    plt.tight_layout()
    plt.savefig(heatmap_mask_filename, dpi=300)
    plt.close()
    print("Bistable mask heatmap saved.")

    # --------------------------------------------------------
    # EXAMPLE TIME COURSES AT MOST BISTABLE POINT
    # --------------------------------------------------------
    if best_indices[0] is not None:
        best_gamma = gamma_vals[best_indices[0]]
        best_beta  = beta_vals[best_indices[1]]
        print("Simulating example time courses at most bistable point...")
        best_params = PDParamsEC2(
            k_E=1.0,
            E_max=1.0,
            k_M=0.5,
            beta=best_beta,
            A=A_value,
            C=C_value,
            gamma=best_gamma
        )

        t_hi, y_hi = simulate_ec2(best_params, y0_high, tmax=500.0, n_points=4000)
        t_lo, y_lo = simulate_ec2(best_params, y0_low,  tmax=500.0, n_points=4000)

        E_hi_tc, M_hi_tc = y_hi
        E_lo_tc, M_lo_tc = y_lo

        plt.figure(figsize=(10, 6))

        plt.subplot(2, 1, 1)
        plt.plot(t_hi, E_hi_tc, label="E (high IC)")
        plt.plot(t_lo, E_lo_tc, label="E (low IC)", linestyle="--")
        plt.ylabel("Energy E")
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(t_hi, M_hi_tc, label="M (high IC)")
        plt.plot(t_lo, M_lo_tc, label="M (low IC)", linestyle="--")
        plt.xlabel("Time")
        plt.ylabel("Mitochondria M")
        plt.legend()

        plt.suptitle(
            f"EC2 Example Time Courses at Most Bistable Point\n"
            f"gamma={best_gamma:.3f}, beta={best_beta:.3f}, A={A_value}, C={C_value}"
        )
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(example_tc_filename, dpi=300)
        plt.close()
        print(f"Example timecourses saved to {example_tc_filename}.")
    else:
        print("No valid best point recorded; skipping example timecourses.")

    print("All done.")
