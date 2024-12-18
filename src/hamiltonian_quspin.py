"""
Build Heisenberg Hamiltonian using QuSpin for small systems.

Note:
 QuSpin's spin_basis_1d expects the spatial length of the lattice
 as the first positional argument, e.g. spin_basis_1d(L, pauli=True).

Requires:
    pip install quspin

Usage:
    python hamiltonian_quspin.py
"""

import numpy as np
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian

def build_heisenberg_hamiltonian(L=6, J=1.0, bc='periodic'):
    """
    Constructs Heisenberg Hamiltonian for a 1D chain of length L (spin-1/2).
    H = J * sum_{i} [Sx_i Sx_{i+1} + Sy_i Sy_{i+1} + Sz_i Sz_{i+1}]
    with optional open or periodic boundary conditions.

    Parameters:
      L : int
         chain length
      J : float
         coupling constant
      bc : str
         'periodic' or 'open' for boundary conditions

    Returns:
      H : quspin operator
         The Heisenberg Hamiltonian
      basis : spin_basis_1d
         The spin basis
    """
    # Create a spin-1/2 basis for 1D chain.
    # Pass L as the first positional argument to spin_basis_1d.
    basis = spin_basis_1d(L, pauli=True)

    # Build bond list
    if bc == 'periodic':
        bonds = [(i, (i+1) % L) for i in range(L)]
    else:
        bonds = [(i, i+1) for i in range(L - 1)]

    # Heisenberg interactions
    # Sx_i Sx_j + Sy_i Sy_j + Sz_i Sz_j
    static_list = []
    static_list.append(["xx", [[J, i, j] for (i, j) in bonds]])
    static_list.append(["yy", [[J, i, j] for (i, j) in bonds]])
    static_list.append(["zz", [[J, i, j] for (i, j) in bonds]])

    # Build Hamiltonian
    H = hamiltonian(static_list, [], basis=basis, dtype=np.float64)
    return H, basis

if __name__ == "__main__":
    # Test usage
    L = 4
    J = 1.0
    bc = 'periodic'
    H, basis = build_heisenberg_hamiltonian(L=L, J=J, bc=bc)
    print(f"Heisenberg Hamiltonian built for L={L}, J={J}, bc={bc}")
    print("Dimension of Hilbert space:", basis.Ns)

    # Optionally do ED
    E_vals, E_vecs = H.eigh()
    print("Eigenvalues:\n", E_vals)
