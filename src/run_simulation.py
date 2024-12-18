"""
Run SSE Simulations Over a Temperature Range & Validate with ED
===============================================================

This script implements:
1. A temperature sweep with SSE for larger system sizes and more sweeps. 
2. An optional validation routine for small systems using exact diagonalization (ED) 
   from QuSpin for a few selected temperatures.

Usage:
  python run_simulation.py

Configuration variables at the top define:
  - L, M, total_sweeps, thermal_sweeps, etc. for the SSE runs.
  - A separate small_system_ED() function runs ED for L=4 or L=6 to 
    cross-check energies against SSE results at a few chosen temperatures.

Output:
  - For SSE runs, results are stored in data/results_temperature_sweep.csv
  - ED validation prints a comparison table to the console.

Dependencies:
  - sse.py, hamiltonian_quspin.py must be in the same folder or python path.
  - QuSpin installed for ED.
  - data directory for saving CSV results.
"""

import os
import numpy as np
import csv

from sse import run_sse_simulation
from hamiltonian_quspin import build_heisenberg_hamiltonian


#####################
# Configuration for SSE
#####################
# Increase system size and sweeps for more representative physics
L = 12                  # chain length
J = 1.0                 # coupling
M = 2 * L               # SSE expansion order
total_sweeps = 10000    # total sweeps
thermal_sweeps = 2000   # sweeps to discard as thermalization
periodic = True         # boundary condition

# Temperature sweep configuration
T_min = 0.5
T_max = 3.0
n_temps = 6  # number of temperatures to sweep
temperatures = np.linspace(T_min, T_max, n_temps)

# Ensure data directory exists
os.makedirs("data", exist_ok=True)

#####################
# SSE Observables Helper
#####################
def compute_observables(energies_sse, energies_sum, mags, beta):
    """
    Compute thermodynamic observables from SSE data:
      - E = <E>
      - E^2 = <E^2>
      - C = ( <E^2> - <E>^2 ) / (T^2 )   [k_B=1]
      - M = <M>
      - Chi = ( <M^2> - <M>^2 ) / T

    We'll use energies_sum as the measure of the system's "true" diagonal energy.
    This is a simplification. A refined SSE approach might rely on expansion-based
    estimators, but for demonstration we'll use local diagonal energies.

    For magnetization, we directly average spin magnetizations. This might be
    pinned if the system gets stuck in certain spin states, requiring more
    advanced loop updates or acceptance rules to reduce metastability.
    """
    # Convert to arrays
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


#####################
# ED Validation for Small Systems
#####################
def exact_diagonalization_energy(L, T, J, bc="periodic"):
    """
    Use hamiltonian_quspin.py to build a Heisenberg Hamiltonian for small L, 
    then compute the exact partition function at temperature T and return 
    the exact thermal average of energy and magnetization (Sz).
    """
    import numpy as np
    H, basis = build_heisenberg_hamiltonian(L=L, J=J, bc=bc)
    E_vals, V = H.eigh()

    beta = 1.0 / T
    Boltz_factors = np.exp(-beta * E_vals)
    Z = np.sum(Boltz_factors)

    # Exact average energy
    E_exact = np.sum(E_vals * Boltz_factors) / Z

    # For magnetization, we compute Sz average with the same Boltzmann weighting
    # We can define a Sz_total operator or sum site by site. For brevity, we'll sum site by site.
    # Each basis state has a definite Sz. We'll do a basis-based approach:
    # basis states -> we can use the fact that basis states have well-defined Sz.

    # We'll do an 'Sz_total' operator in QuSpin:
    from quspin.operators import hamiltonian

    static = [["z", [[1.0, i] for i in range(L)]]]  # sum of Sz_i
    Sz_total = hamiltonian(static, [], basis=basis, dtype=np.float64)
    # Evaluate <Sz_total> exactly
    # Because we have E_vecs, we can do a diagonalization approach or a direct approach.
    # We'll do a direct approach: for each eigenvalue/eigenvector, compute <psi|Sz_total|psi>.
    # Then weight by Boltzmann factor.

    # We can do (psi^dagger Sz_total psi) for each eigenstate
    # E_vecs has shape (dim, dim). columns are eigenvectors
    # We'll use 'Sz_total.matrix_ele' or do the matrix multiplication.

    # We'll do a simpler approach: build the diagonal representation or just use the built-in: 
    # 'Sz_total.expt_value(V, i)' -> but that only does for a single state index i if I recall.
    # We'll do it manually:

    n_states = E_vals.shape[0]
    mag_n = np.zeros(n_states, dtype=float)
    for i in range(n_states):
        # eigenvector i:
        psi_i = V[:, i]
        # compute <psi_i|Sz_total|psi_i>
        # We can do a dot(psi_i.conj(), Sz_total.dot(psi_i)) or use H.matrix_ele
        M_i = np.vdot(psi_i, Sz_total.dot(psi_i)).real
        mag_n[i] = M_i

    M_exact = np.sum(mag_n * Boltz_factors) / Z

    return E_exact, M_exact


def validate_small_system_ed():
    """
    Compare SSE results vs exact diagonalization for L=4 or L=6, 
    at a few temperature points. 
    Print a table of T, E_SSE, E_ED, M_SSE, M_ED, etc.
    """

    L_small = 4
    J_val = 1.0
    bc = "periodic"
    M_val = 2 * L_small
    sweeps = 5000
    thermal = 1000
    temperatures_test = [0.5, 1.0, 2.0]  # example test points

    print("\nValidation for L={} using exact diagonalization:\n".format(L_small))
    print(" T    E_SSE     E_ED      M_SSE     M_ED")

    for T in temperatures_test:
        beta = 1.0 / T
        sse_data = run_sse_simulation(
            num_sites=L_small,
            beta=beta,
            M=M_val,
            num_sweeps=sweeps,
            J=J_val,
            periodic=(bc == "periodic")
        )

        e_sum_sse = sse_data["energy_sum"][thermal:]
        m_sse = sse_data["magnetization"][thermal:]
        E_sse = np.mean(e_sum_sse)
        M_sse = np.mean(m_sse)

        # ED
        E_ed, M_ed = exact_diagonalization_energy(L_small, T, J_val, bc=bc)

        print(" {T:4.1f}  {E_sse:8.4f}  {E_ed:8.4f}  {M_sse:8.4f}  {M_ed:8.4f}".format(
            T=T, E_sse=E_sse, E_ed=E_ed, M_sse=M_sse, M_ed=M_ed
        ))


#####################
# Main
#####################
def main():
    """
    1) Run SSE for a larger L=12 system across a temperature range.
    2) Save results to CSV.
    3) (Optional) Validate for small L=4 with ED.
    """
    # Part 1-2: SSE with bigger system
    results_list = []

    for T in temperatures:
        beta = 1.0 / T
        sim_data = run_sse_simulation(L, beta, M, total_sweeps, J=J, periodic=periodic)
        energies_sse = sim_data["energy_sse"][thermal_sweeps:]
        energies_sum = sim_data["energy_sum"][thermal_sweeps:]
        mags = sim_data["magnetization"][thermal_sweeps:]

        E_avg, spec_heat, M_avg, chi = compute_observables(energies_sse, energies_sum, mags, beta)
        results_list.append({"T": T, "E": E_avg, "C": spec_heat, "M": M_avg, "Chi": chi})
        print(f"L={L}, T={T:.3f}, E={E_avg:.4f}, C={spec_heat:.4f}, M={M_avg:.4f}, χ={chi:.4f}")

    # Save to CSV
    csv_path = os.path.join("data", f"results_temperature_sweep_L{L}.csv")
    with open(csv_path, "w", newline="") as csvfile:
        fieldnames = ["T", "E", "C", "M", "Chi"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results_list:
            writer.writerow(row)

    print(f"\nData has been saved to: {csv_path}")
    print("Now performing validation on a small system via ED...")

    # Part 3: Validate for L=4
    validate_small_system_ed()

if __name__ == "__main__":
    main()
