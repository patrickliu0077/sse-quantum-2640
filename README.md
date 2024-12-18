# Finite-Temperature Quantum Monte Carlo Simulations of 1D Heisenberg Spin Systems

This project implements the Stochastic Series Expansion (SSE) quantum Monte Carlo method to study finite-temperature properties of the 1D Heisenberg spin chain. The implementation focuses on thermodynamic properties like energy, specific heat, magnetization, and susceptibility across varying temperatures and system sizes.

## Features

### Core Implementation
- Stochastic Series Expansion (SSE) quantum Monte Carlo algorithm
- Periodic boundary conditions for the 1D Heisenberg model
- Efficient operator string updates with diagonal and off-diagonal moves
- Support for various system sizes (L=16 to 64)
- Temperature range from T=0.5J to 3.0J

### Analysis Capabilities
- Finite-size scaling analysis
- Low-temperature quantum effects
- Thermodynamic observables:
  * Energy per site (E/NJ)
  * Specific Heat (C/kB)
  * Magnetization (M/Msat)
  * Magnetic Susceptibility (χJ)
- Publication-quality plotting

### Technical Features
- Python implementation with NumPy/SciPy
- QuSpin integration for exact diagonalization comparisons
- Efficient data storage and analysis
- Modular code structure for easy modification

## Installation

1. Ensure Python 3.x is installed
2. Install required packages:
```bash
pip install numpy scipy matplotlib quspin
```

## Project Structure

```
project_root/
├── src/                    # Source code
│   ├── sse.py             # Core SSE implementation
│   ├── hamiltonian_quspin.py  # QuSpin interface
│   ├── run_simulation.py  # Main simulation runner
│   └── plot_results.py    # Basic plotting utilities
├── analysis/              # Analysis tools
│   ├── finite_size_analysis.py  # Finite-size scaling
│   ├── convergence_analysis.py  # Convergence studies
│   ├── statistics.py      # Statistical analysis
│   └── plot_finite_size.py      # Publication plots
├── data/                  # Simulation results
├── results/               # Generated plots and analysis
└── docs/                  # Documentation
```

## Usage Guide

### Running Simulations

1. Basic simulation with default parameters:
```bash
python src/run_simulation.py
```

2. Finite-size analysis with multiple system sizes:
```bash
python analysis/finite_size_analysis.py
```
This runs simulations for L = [16, 24, 32, 48, 64] with:
- 30,000 total Monte Carlo sweeps
- 6,000 thermalization sweeps
- 20 temperature points from T=0.5 to T=3.0

### Analyzing Results

1. Generate publication-quality plots:
```bash
python analysis/plot_finite_size.py
```
Creates plots in results/analysis_[timestamp]/:
- Energy vs Temperature
- Specific Heat vs Temperature
- Magnetization vs Temperature
- Susceptibility vs Temperature

2. Analyze convergence and scaling:
```bash
python analysis/convergence_analysis.py
```
Generates:
- Finite-size scaling plots
- Low-temperature behavior analysis
- System size convergence studies

### Data Format

Simulation results are stored in CSV format with columns:
- L: System size
- T: Temperature
- E: Energy per site
- C: Specific heat
- M: Magnetization
- Chi: Susceptibility

## Physical Results

The implementation successfully captures key physical features of the 1D Heisenberg model:

1. Energy:
- Systematic size dependence at low T
- Proper convergence at high T
- Smooth curves indicating good statistics

2. Specific Heat:
- Peak height scales with system size
- Clear low-temperature quantum effects
- Peak position shows systematic behavior

3. Magnetization:
- Fluctuates around zero (correct for Heisenberg model)
- Reduced fluctuations for larger systems
- Demonstrates quantum disorder

4. Susceptibility:
- Systematic decrease with system size
- Clear finite-size scaling
- Proper high-temperature behavior

## Performance Considerations

- Memory usage scales with system size L and expansion order M
- Computation time increases with:
  * System size L
  * Inverse temperature β
  * Number of Monte Carlo sweeps
- Parallel execution possible for different temperatures/sizes

## References

1. Sylju˚asen & Sandvik, "Quantum Monte Carlo with directed loops", Physical Review E, 66(4), 046701, 2002.
2. Sandvik, "Computational Studies of Quantum Spin Systems", AIP Conference Proceedings 1297, 135, 2010.

## Contributing

Feel free to submit issues, fork the repository, and create pull requests for any improvements.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
