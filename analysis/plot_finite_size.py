"""
Plot Finite-Size Analysis Results with Publication-Quality Styling
================================================================

Reads data from data/finite_size_results_extended.csv and produces publication-quality
plots of thermodynamic observables vs temperature for different system sizes.

Features:
- Professional color scheme
- LaTeX-style labels and fonts
- Consistent styling across all plots
- High-resolution output
- Proper axis scaling and grid styling
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
from datetime import datetime

# Create timestamped output directory
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
output_dir = os.path.join('results', f'analysis_{timestamp}')
os.makedirs(output_dir, exist_ok=True)

# Set up publication-quality plot styling
plt.style.use('default')
mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'text.usetex': True,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 16,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'lines.linewidth': 2.0,
    'lines.markersize': 8,
    'axes.linewidth': 1.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.prop_cycle': cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
})

def plot_finite_size(path="data/finite_size_results_publication.csv"):
    """Generate publication-quality plots from finite-size analysis data."""
    
    print(f"Reading data from: {path}")
    
    # Read and organize data
    data = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                "L": float(row["L"]),
                "T": float(row["T"]),
                "E": float(row["E"]),
                "C": float(row["C"]),
                "M": float(row["M"]),
                "Chi": float(row["Chi"]),
            })

    # Group by L
    data_by_L = {}
    for d in data:
        L = d["L"]
        if L not in data_by_L:
            data_by_L[L] = []
        data_by_L[L].append(d)
    
    # Sort each group by T
    for L in data_by_L:
        data_by_L[L].sort(key=lambda x: x["T"])

    # Plot Energy vs T
    fig, ax = plt.subplots(figsize=(8, 6))
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Evals = [r["E"] for r in rows]
        ax.plot(Tvals, Evals, 'o-', label=f'$L={int(L)}$', markersize=8)

    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Energy per site $(E/NJ)$')
    ax.set_title(r'Energy vs Temperature: Finite-Size Analysis')
    ax.legend(frameon=True, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "energy.pdf"))
    plt.close()

    # Plot Specific Heat
    fig, ax = plt.subplots(figsize=(8, 6))
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Cvals = [r["C"] for r in rows]
        ax.plot(Tvals, Cvals, 's-', label=f'$L={int(L)}$', markersize=8)

    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Specific Heat $(C/k_B)$')
    ax.set_title(r'Specific Heat vs Temperature')
    ax.legend(frameon=True, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "specific_heat.pdf"))
    plt.close()

    # Plot Magnetization
    fig, ax = plt.subplots(figsize=(8, 6))
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Mvals = [r["M"] for r in rows]
        ax.plot(Tvals, Mvals, '^-', label=f'$L={int(L)}$', markersize=8)

    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Magnetization $(M/M_{\mathrm{sat}})$')
    ax.set_title(r'Magnetization vs Temperature')
    ax.legend(frameon=True, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "magnetization.pdf"))
    plt.close()

    # Plot Susceptibility
    fig, ax = plt.subplots(figsize=(8, 6))
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Chivals = [r["Chi"] for r in rows]
        ax.plot(Tvals, Chivals, 'd-', label=f'$L={int(L)}$', markersize=8)

    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Susceptibility $(\chi J)$')
    ax.set_title(r'Magnetic Susceptibility vs Temperature')
    ax.legend(frameon=True, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "susceptibility.pdf"))
    plt.close()

    # Create a combined plot (2x2 subplots)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    # Energy (top left)
    ax = axes[0,0]
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Evals = [r["E"] for r in rows]
        ax.plot(Tvals, Evals, 'o-', label=f'$L={int(L)}$', markersize=6)
    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Energy per site $(E/NJ)$')
    ax.legend(frameon=True, fancybox=True)

    # Specific Heat (top right)
    ax = axes[0,1]
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Cvals = [r["C"] for r in rows]
        ax.plot(Tvals, Cvals, 's-', label=f'$L={int(L)}$', markersize=6)
    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Specific Heat $(C/k_B)$')
    ax.legend(frameon=True, fancybox=True)

    # Magnetization (bottom left)
    ax = axes[1,0]
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Mvals = [r["M"] for r in rows]
        ax.plot(Tvals, Mvals, '^-', label=f'$L={int(L)}$', markersize=6)
    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Magnetization $(M/M_{\mathrm{sat}})$')
    ax.legend(frameon=True, fancybox=True)

    # Susceptibility (bottom right)
    ax = axes[1,1]
    for L, rows in sorted(data_by_L.items()):
        Tvals = [r["T"] for r in rows]
        Chivals = [r["Chi"] for r in rows]
        ax.plot(Tvals, Chivals, 'd-', label=f'$L={int(L)}$', markersize=6)
    ax.set_xlabel(r'Temperature $(k_B T/J)$')
    ax.set_ylabel(r'Susceptibility $(\chi J)$')
    ax.legend(frameon=True, fancybox=True)

    fig.suptitle(r'Thermodynamic Properties of 1D Heisenberg Model', fontsize=16)
    plt.savefig(os.path.join(output_dir, "combined.pdf"))
    plt.close()

    print(f"\nPlots saved in directory: {output_dir}")
    print("Generated PDFs:")
    print("  - Individual plots: energy.pdf, specific_heat.pdf, magnetization.pdf, susceptibility.pdf")
    print("  - Combined plot: combined.pdf")

if __name__ == "__main__":
    plot_finite_size()
