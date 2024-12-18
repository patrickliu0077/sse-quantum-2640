"""
Finite-Size Analysis for the 1D Heisenberg Model using SSE
==========================================================

Enhanced version with publication-quality system sizes.
This script systematically investigates finite-size effects by running
the SSE simulation at multiple system sizes L over a temperature range,
collecting data for each (L, T) in a single output CSV.
You can then plot or analyze how the observables (energy, magnetization, etc.) 
vary with L to discuss finite-size scaling.

Usage:
  python analysis/finite_size_analysis.py

Configuration:
  - sizes: Publication-quality system sizes (L=16,24,32,48,64)
  - T_min, T_max, n_temps: Temperature sweep parameters
  - total_sweeps, thermal_sweeps: Monte Carlo parameters
  - The data is stored to "data/finite_size_results.csv"
"""

import os
import sys
import numpy as np
import csv

# Add the src directory to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from sse import run_sse_simulation

# Publication-quality system sizes
sizes = [16, 24, 32, 48, 64]  # Standard sizes for 1D quantum spin chains
T_min = 0.5
T_max = 3.0
n_temps = 20  # Fine temperature resolution
temperatures = np.linspace(T_min, T_max, n_temps)

# Monte Carlo parameters - increased for larger systems
J = 1.0
total_sweeps = 30000   # More sweeps for better statistics
thermal_sweeps = 6000  # Longer thermalization
periodic = True

# Output CSV
os.makedirs("data", exist_ok=True)
output_path = os.path.join("data", "finite_size_results_publication.csv")

def compute_observables(energies_sse, energies_sum, mags, beta):
    """
    Compute observables from measurement time series.
    """
    e_array = np.array(energies_sum, dtype=float)
    m_array = np.array(mags, dtype=float)
    E_avg = np.mean(e_array)
    E2_avg = np.mean(e_array**2)
    M_avg = np.mean(m_array)
    M2_avg = np.mean(m_array**2)

    T = 1.0 / beta
    C = (E2_avg - E_avg**2) / (T**2)
    Chi = (M2_avg - M_avg**2) / T
    return E_avg, C, M_avg, Chi

def main():
    results = []
    for L in sizes:
        print(f"\nStarting calculations for L={L}")
        M = 2 * L  # expansion order
        for T in temperatures:
            beta = 1.0 / T
            print(f"  Temperature T={T:.3f}...", end=" ", flush=True)

            # Run SSE
            sim_data = run_sse_simulation(
                num_sites=L,
                beta=beta,
                M=M,
                num_sweeps=total_sweeps,
                J=J,
                periodic=periodic
            )

            # Discard thermal sweeps
            energies_sse = sim_data["energy_sse"][thermal_sweeps:]
            energies_sum = sim_data["energy_sum"][thermal_sweeps:]
            mags = sim_data["magnetization"][thermal_sweeps:]

            E_avg, spec_heat, M_avg, chi = compute_observables(
                energies_sse, energies_sum, mags, beta
            )

            results.append({
                "L": L,
                "T": T,
                "E": E_avg,
                "C": spec_heat,
                "M": M_avg,
                "Chi": chi
            })
            print(f"E={E_avg:.4f}, C={spec_heat:.4f}, M={M_avg:.4f}, χ={chi:.4f}")

    # Write to CSV
    fieldnames = ["L", "T", "E", "C", "M", "Chi"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nPublication-quality finite-size analysis results saved to {output_path}")
    print("\nSystem sizes used:")
    print(f"L = {sizes}")
    print("\nNext steps:")
    print("1. Run 'python analysis/plot_finite_size.py' to create publication plots")
    print("2. Check convergence of observables with system size")
    print("3. Analyze low-temperature behavior and finite-size effects")

if __name__ == "__main__":
    main()
