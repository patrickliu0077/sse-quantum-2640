"""
Plot Results from SSE Simulations
=================================

This script reads the CSV files produced by run_simulation.py and generates
Matplotlib plots of thermodynamic observables vs. temperature.

Usage:
  python src/plot_results.py

Configuration:
  - data_file: path to the CSV with columns [T, E, C, M, Chi].
  - Output: A PNG or PDF file with subplots for E(T), C(T), M(T), and Chi(T).

Feel free to tweak plotting styles, colors, marker types, etc.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

def plot_sse_results(data_file, output_file="results_plot.png"):
    """
    Reads CSV data with columns T, E, C, M, Chi and plots them as separate subplots.
    """
    T_vals = []
    E_vals = []
    C_vals = []
    M_vals = []
    Chi_vals = []

    with open(data_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            T_vals.append(float(row["T"]))
            E_vals.append(float(row["E"]))
            C_vals.append(float(row["C"]))
            M_vals.append(float(row["M"]))
            Chi_vals.append(float(row["Chi"]))

    T_vals = np.array(T_vals)
    E_vals = np.array(E_vals)
    C_vals = np.array(C_vals)
    M_vals = np.array(M_vals)
    Chi_vals = np.array(Chi_vals)

    fig, axes = plt.subplots(2, 2, figsize=(10,8))
    fig.suptitle("SSE Results for 1D Heisenberg Chain")

    # Subplot 1: Energy vs T
    ax = axes[0, 0]
    ax.plot(T_vals, E_vals, marker='o', color='blue', label='Energy')
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("Energy")
    ax.grid(True)
    ax.legend()

    # Subplot 2: Specific Heat vs T
    ax = axes[0, 1]
    ax.plot(T_vals, C_vals, marker='s', color='red', label='Specific Heat')
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("C(T)")
    ax.grid(True)
    ax.legend()

    # Subplot 3: Magnetization vs T
    ax = axes[1, 0]
    ax.plot(T_vals, M_vals, marker='^', color='green', label='Magnetization')
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("Magnetization")
    ax.grid(True)
    ax.legend()

    # Subplot 4: Susceptibility vs T
    ax = axes[1, 1]
    ax.plot(T_vals, Chi_vals, marker='d', color='purple', label='Susceptibility')
    ax.set_xlabel("Temperature (T)")
    ax.set_ylabel("χ(T)")
    ax.grid(True)
    ax.legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # leave space for suptitle
    plt.savefig(output_file, dpi=150)
    print(f"Plots saved to {output_file}")

if __name__ == "__main__":
    # Example usage
    data_file = os.path.join("data", "results_temperature_sweep_L12.csv")
    output_file = "results_plot_L12.png"
    plot_sse_results(data_file, output_file)
