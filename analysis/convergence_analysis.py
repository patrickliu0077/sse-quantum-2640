"""
Analyze Convergence and Low-Temperature Behavior
==============================================

This script analyzes the convergence of observables with system size
and examines the low-temperature behavior of the 1D Heisenberg model.

Features:
1. Plots observables vs 1/L to examine finite-size scaling
2. Focuses on low-T behavior where quantum effects dominate
3. Estimates infinite-size limits of observables
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime

# Create timestamped output directory
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
output_dir = os.path.join('results', f'convergence_{timestamp}')
os.makedirs(output_dir, exist_ok=True)

# Publication-quality plot settings
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
})

def load_data(filepath):
    """Load and organize data from CSV."""
    print(f"Reading data from: {filepath}")
    data = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                'L': float(row['L']),
                'T': float(row['T']),
                'E': float(row['E']),
                'C': float(row['C']),
                'M': float(row['M']),
                'Chi': float(row['Chi'])
            })
    return data

def plot_finite_size_scaling(data):
    """
    Plot 1/L scaling for various observables at fixed temperatures.
    This helps visualize how quickly we approach the thermodynamic limit.
    """
    # Group data by temperature
    temps = sorted(list(set(d['T'] for d in data)))
    sizes = sorted(list(set(d['L'] for d in data)))
    inv_L = [1/L for L in sizes]

    # Select a few representative temperatures
    selected_temps = [temps[0], temps[len(temps)//4], temps[len(temps)//2]]

    # Plot E vs 1/L
    fig, ax = plt.subplots(figsize=(8, 6))
    for T in selected_temps:
        E_vals = [d['E'] for d in data if abs(d['T'] - T) < 1e-6]  # Handle floating point comparison
        ax.plot(inv_L, E_vals, 'o-', label=f'$T/J={T:.2f}$')
    
    ax.set_xlabel(r'$1/L$')
    ax.set_ylabel(r'Energy per site $(E/NJ)$')
    ax.set_title(r'Finite-Size Scaling of Energy')
    ax.legend(frameon=True, fancybox=True)
    plt.savefig(os.path.join(output_dir, 'energy_scaling.pdf'))
    plt.close()

    # Plot C vs 1/L
    fig, ax = plt.subplots(figsize=(8, 6))
    for T in selected_temps:
        C_vals = [d['C'] for d in data if abs(d['T'] - T) < 1e-6]
        ax.plot(inv_L, C_vals, 's-', label=f'$T/J={T:.2f}$')
    
    ax.set_xlabel(r'$1/L$')
    ax.set_ylabel(r'Specific Heat $(C/k_B)$')
    ax.set_title(r'Finite-Size Scaling of Specific Heat')
    ax.legend(frameon=True, fancybox=True)
    plt.savefig(os.path.join(output_dir, 'specificheat_scaling.pdf'))
    plt.close()

def analyze_low_T_behavior(data):
    """
    Focus on the low-temperature regime where quantum effects dominate.
    Plot observables vs T for the largest system sizes.
    """
    # Get largest two system sizes
    sizes = sorted(list(set(d['L'] for d in data)))
    largest_sizes = sizes[-2:]

    # Filter for low temperatures
    temps = sorted(list(set(d['T'] for d in data)))
    low_T_cutoff = temps[len(temps)//3]  # Use bottom third of temperature range

    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    for L in largest_sizes:
        # Get data for this size
        size_data = [d for d in data if d['L'] == L and d['T'] <= low_T_cutoff]
        size_data.sort(key=lambda x: x['T'])  # Sort by temperature
        T_vals = [d['T'] for d in size_data]
        
        # Energy
        E_vals = [d['E'] for d in size_data]
        ax1.plot(T_vals, E_vals, 'o-', label=f'$L={int(L)}$')
        ax1.set_xlabel(r'$T/J$')
        ax1.set_ylabel(r'$E/NJ$')
        ax1.set_title('Low-T Energy')
        
        # Specific Heat
        C_vals = [d['C'] for d in size_data]
        ax2.plot(T_vals, C_vals, 's-', label=f'$L={int(L)}$')
        ax2.set_xlabel(r'$T/J$')
        ax2.set_ylabel(r'$C/k_B$')
        ax2.set_title('Low-T Specific Heat')
        
        # Magnetization
        M_vals = [d['M'] for d in size_data]
        ax3.plot(T_vals, M_vals, '^-', label=f'$L={int(L)}$')
        ax3.set_xlabel(r'$T/J$')
        ax3.set_ylabel(r'$M/M_{\mathrm{sat}}$')
        ax3.set_title('Low-T Magnetization')
        
        # Susceptibility
        Chi_vals = [d['Chi'] for d in size_data]
        ax4.plot(T_vals, Chi_vals, 'd-', label=f'$L={int(L)}$')
        ax4.set_xlabel(r'$T/J$')
        ax4.set_ylabel(r'$\chi J$')
        ax4.set_title('Low-T Susceptibility')
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.legend(frameon=True, fancybox=True)
        ax.grid(True, alpha=0.3)
    
    fig.suptitle('Low-Temperature Behavior of 1D Heisenberg Model', fontsize=16)
    plt.savefig(os.path.join(output_dir, 'lowT_analysis.pdf'))
    plt.close()

def main():
    # Load the extended finite-size data
    data_path = os.path.join('data', 'finite_size_results_publication.csv')
    data = load_data(data_path)
    
    # Create plots
    plot_finite_size_scaling(data)
    analyze_low_T_behavior(data)
    
    print(f"\nAnalysis complete. Results saved in: {output_dir}")
    print("Generated plots:")
    print("1. energy_scaling.pdf - Energy vs 1/L scaling")
    print("2. specificheat_scaling.pdf - Specific Heat vs 1/L scaling")
    print("3. lowT_analysis.pdf - Detailed low-temperature behavior")

if __name__ == '__main__':
    main()
