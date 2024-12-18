"""
SSE Implementation with a (Simplified) Directed Loop Update
===========================================================

Overview
--------
This file extends the previous diagonal-update-based SSE approach by adding a rudimentary
Directed Loop update to handle off-diagonal spin-flip processes. In a full Heisenberg model,
the off-diagonal part arises from S^x S^x + S^y S^y terms. The SSE representation includes
operators that can flip neighboring spins. The "Directed Loop" or "Worm" algorithm ensures
ergodic sampling among all spin configurations.

Key Concepts:
-------------
1. We define two types of operators for each bond:
   - Diagonal:  SzSz
   - Off-Diagonal: S+S- + S-S+ (collectively responsible for flipping spins).

2. The Directed Loop algorithm:
   - Picks a random entry point (p in operator string), a random site. 
   - Attempts to propagate a "worm" through the operator string, flipping spins and transitioning
     diagonal operators to off-diagonal or vice versa. This is done while preserving detailed
     balance and ensuring the correct SSE weights.

3. Simplifications:
   - We only demonstrate the structure of a directed loop approach. A full correct
     implementation of directed loops can be significantly more involved, particularly for
     large systems or advanced features (like directed loop in the extended subspace). We'll
     outline a minimal version that shows how spins might be flipped as we traverse the
     operator string.
   - The exact acceptance probabilities for transitions can be found in references such as
     Syljuåsen & Sandvik (2002). Here, we provide a simplified approach to illustrate how
     one might insert off-diagonal operators that flip spins.

Note: To fully validate the loop method, you’d need to test on small systems and compare with
exact diagonalization. This code should serve as a starting point for a more robust SSE
off-diagonal update.

References:
-----------
- A. W. Sandvik, R. P. Singh, "Quantum Monte Carlo Methods," Lecture Notes, (various).
- O. F. Syljuåsen & A. W. Sandvik, "Quantum Monte Carlo with Directed Loops," Phys. Rev. E 66, 046701 (2002).
"""

import numpy as np

################################################################################
# Utility functions
################################################################################

def generate_bonds(num_sites, periodic=True):
    """
    Generate a list of nearest-neighbor bonds for a 1D chain.
    If periodic is True, connects site (num_sites-1) to site 0.
    """
    bonds = [(i, i+1) for i in range(num_sites - 1)]
    if periodic:
        bonds.append((num_sites - 1, 0))
    return bonds

def bond_energy_diagonal(spins, bond, J=1.0):
    """
    Diagonal term SzSz for the Heisenberg model. Each spin is ±1 for ±1/2 * 2.
    E_loc = J * spins[i]*spins[j], ignoring the 1/4 factor for simplicity.
    """
    i, j = bond
    return J * spins[i] * spins[j]

def can_flip(spins, bond):
    """
    For an off-diagonal operator to act, the spins on a bond must be opposite
    (e.g., +1, -1) for S+S- or S-S+ to flip them. Return True if they are flippable.
    """
    i, j = bond
    return spins[i] != spins[j]  # e.g., +1,-1 or -1,+1

def flip_spins(spins, bond):
    """
    Flip the spins on the given bond. This models an S+S- or S-S+ action.
    """
    i, j = bond
    spins[i] *= -1
    spins[j] *= -1


################################################################################
# SSE Operator String Representation
################################################################################
# We'll define possible operator types:
#  - "I": Identity
#  - "SzSz": Diagonal operator (Sz_i Sz_j)
#  - "Flip": Off-diagonal operator that flips spins on a bond (S+S- + S-S+)
################################################################################

def initialize_spins(num_sites):
    """Randomly initialize spin configuration in ±1 for each site."""
    spins = np.random.choice([+1, -1], size=num_sites)
    return spins

def initialize_operator_string(M):
    """
    Initialize operator string with length M, all set to identity:
        Each element is (op_type, bond_index).
    """
    return [("I", None) for _ in range(M)]


################################################################################
# Diagonal Update
################################################################################

def diagonal_update(operator_string, spins, beta, bonds, J=1.0):
    """
    SSE diagonal update (more accurate ratio than naive approach).
    For each operator position, we try to insert or remove a SzSz operator.
    """
    M = len(operator_string)
    n_non_identity = sum(1 for (op, _) in operator_string if op != "I")
    Nbonds = len(bonds)

    for p in range(M):
        op_type, b_idx = operator_string[p]
        # Recount n_non_identity for clarity
        n_non_identity = sum(1 for (op, _) in operator_string if op != "I")

        if op_type == "I":
            # Insert SzSz
            if n_non_identity < M:
                proposed_bond_idx = np.random.randint(Nbonds)
                e_loc = bond_energy_diagonal(spins, bonds[proposed_bond_idx], J=J)
                localWeight = np.exp(-beta * e_loc)

                numerator = beta * Nbonds * localWeight
                denominator = (M - n_non_identity)
                p_accept = min(1.0, numerator / denominator)

                if np.random.rand() < p_accept:
                    operator_string[p] = ("SzSz", proposed_bond_idx)

        elif op_type == "SzSz":
            # Remove SzSz
            e_loc = bond_energy_diagonal(spins, bonds[b_idx], J=J)
            localWeight = np.exp(-beta * e_loc)

            numerator = (M - (n_non_identity - 1))
            denominator = beta * Nbonds * localWeight
            p_accept = min(1.0, numerator / denominator)

            if np.random.rand() < p_accept:
                operator_string[p] = ("I", None)


################################################################################
# Directed Loop Update
################################################################################

def directed_loop_update(operator_string, spins, beta, bonds, J=1.0):
    """
    Simplified directed loop approach for off-diagonal operators:
      1. Choose a random position p in [0, M-1] and a random site "start_site".
      2. Move forward through the operator string, deciding whether to flip spins
         and transform diagonal ops <-> Flip ops if feasible. 
      3. Possibly wrap around to the beginning of the string if needed (cyclic).
      4. End when we return to (p, start_site).

    This is a minimal demonstration of the idea. A fully correct directed loop
    approach has more detailed acceptance steps matching the local operators,
    spin states, and SSE combinatorics from references like Syljuåsen & Sandvik.

    We'll do:
      - If we see "I", we might attempt to insert a "Flip" operator if the bond is flippable.
      - If we see "Flip", we can remove it or keep traveling. 
      - We track which bond is relevant by the local site or the bond in question.

    WARNING: This is not a complete, validated directed loop. It's an illustrative
    structure only to show how one might traverse the operator string flipping spins.
    """

    M = len(operator_string)
    if M == 0:
        return

    # 1. Choose random entry point and random site
    p = np.random.randint(M)
    start_site = np.random.randint(len(spins))

    # We store "current_site" in the loop
    current_site = start_site
    visited_count = 0

    while True:
        visited_count += 1
        if visited_count > 2 * M:
            # Bail out if we've looped too many times – safety check
            break

        op_type, bond_idx = operator_string[p]

        if op_type == "I":
            # Attempt to insert a Flip operator if the chosen bond is flippable
            # We'll pick a random bond that includes current_site
            possible_bonds = []
            for b_i, (s1, s2) in enumerate(bonds):
                if s1 == current_site or s2 == current_site:
                    possible_bonds.append(b_i)
            if possible_bonds:
                chosen_bond_idx = np.random.choice(possible_bonds)
                if can_flip(spins, bonds[chosen_bond_idx]):
                    # Acceptance ratio in a real directed loop is more involved
                    # We'll do a simplistic approach:
                    operator_string[p] = ("Flip", chosen_bond_idx)
                    # Immediately flip spins
                    flip_spins(spins, bonds[chosen_bond_idx])
                    # Update current_site to the other spin in that bond
                    s1, s2 = bonds[chosen_bond_idx]
                    current_site = s2 if current_site == s1 else s1

        elif op_type == "Flip":
            # Attempt to remove or keep traveling
            # Typically we also consider a scattering event to choose new direction
            # We'll remove it with some probability for illustration
            remove_prob = 0.5  # placeholder
            if np.random.rand() < remove_prob:
                # Remove the operator
                bond_to_unflip = bond_idx
                # Flip the spins back
                flip_spins(spins, bonds[bond_to_unflip])
                operator_string[p] = ("I", None)
                # current_site remains the same. 
                # In a real approach we'd handle loop direction changes

            else:
                # "Scatter" or continue traveling in the operator string
                # We'll flip spins again to model a "bounce" 
                flip_spins(spins, bonds[bond_idx])
                # Move current site
                s1, s2 = bonds[bond_idx]
                current_site = s2 if current_site == s1 else s1

        else:
            # If diagonal or anything else, we just continue
            # In a real approach, we consider if the diagonal operator affects the loop
            pass

        # Move to the next position in the string
        p = (p + 1) % M

        # Check if we've returned to the starting condition
        if p == 0 and current_site == start_site:
            break


################################################################################
# Measurement
################################################################################

def measure_sse_energy(operator_string, beta, bonds, spins, J=1.0):
    """
    SSE-based energy estimate:
      - n_diagonal ~ average number of diagonal operators (SzSz)
      - Off-diagonal operators also contribute to the expanded Hamiltonian
      In a full approach, we parse both diagonal and off-diagonal operators 
      to compute <n> / beta, etc. But we'll keep it simple.

    We'll do:
      E_sse ~ (n_diag / beta) ignoring constant shifts.
      e_sum = sum of diagonal bond energies from "SzSz" plus ignoring "Flip" 
      for now in the sum. 
    """
    M = len(operator_string)
    n_diagonal = sum(1 for (op, _) in operator_string if op == "SzSz")

    # Basic SSE formula
    e_sse = n_diagonal / beta

    # Sum local diagonal energies
    e_sum = 0.0
    for (op, b_idx) in operator_string:
        if op == "SzSz" and b_idx is not None:
            e_sum += bond_energy_diagonal(spins, bonds[b_idx], J=J)

    return e_sse, e_sum

def measure_magnetization(spins):
    """Average magnetization Mz."""
    return np.mean(spins)


################################################################################
# Main SSE Simulation
################################################################################

def run_sse_simulation(num_sites, beta, M, num_sweeps, J=1.0, periodic=True):
    """
    Main SSE simulation orchestrator:
        - Initialize spins and SSE operator string
        - For each sweep:
            1) Diagonal update
            2) Directed loop update (off-diagonal)
            3) Measure observables
    """
    spins = initialize_spins(num_sites)
    operator_str = initialize_operator_string(M)
    bonds = generate_bonds(num_sites, periodic=periodic)

    energies_sse = []
    energies_sum = []
    mags = []

    for _ in range(num_sweeps):
        diagonal_update(operator_str, spins, beta, bonds, J=J)
        directed_loop_update(operator_str, spins, beta, bonds, J=J)

        e_sse, e_sum = measure_sse_energy(operator_str, beta, bonds, spins, J=J)
        m = measure_magnetization(spins)

        energies_sse.append(e_sse)
        energies_sum.append(e_sum)
        mags.append(m)

    results = {
        "energy_sse": energies_sse,  # SSE-based estimator from n_diagonal
        "energy_sum": energies_sum,  # local sum of diagonal energies
        "magnetization": mags,
    }

    return results


if __name__ == "__main__":
    # Example usage:
    L = 6           # chain length
    beta = 1.0      # inverse temperature
    M = 2 * L       # expansion order
    sweeps = 200    # number of MC sweeps
    J = 1.0
    periodic = True

    sim_results = run_sse_simulation(L, beta, M, sweeps, J=J, periodic=periodic)
    e_sse_avg = np.mean(sim_results["energy_sse"])
    e_sum_avg = np.mean(sim_results["energy_sum"])
    m_avg = np.mean(sim_results["magnetization"])

    print(f"Simplified SSE with Directed Loop example:")
    print(f"  SSE-based energy estimator = {e_sse_avg:.4f}")
    print(f"  Sum of diagonal energies   = {e_sum_avg:.4f}")
    print(f"  Avg magnetization         = {m_avg:.4f}")
