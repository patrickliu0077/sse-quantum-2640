"""
Robust Statistical Analysis for SSE Data
========================================

This module demonstrates a simple Jackknife approach to estimate statistical 
uncertainties for observables (e.g., energy, magnetization) in the SSE data. 
You can extend this or use alternatives like Bootstrap.

Usage:
  python statistics.py

Description:
  1. Read an SSE CSV file (with columns T, E, C, M, Chi) for each temperature.
  2. Perform a Jackknife resampling for each observable (E, M, etc.) to estimate
     standard errors.
  3. Output the results to a new CSV with columns T, E, dE, C, dC, M, dM, Chi, dChi.

Command-line arguments or a config can be added for controlling file paths.
"""

import os
import csv
import numpy as np

def jackknife(data):
    """
    Given a 1D NumPy array of samples 'data', return the jackknife 
    average and standard error.

    The jackknife procedure:
      - Let N be the number of samples.
      - For i in 0..N-1, compute the mean of all samples except data[i].
      - The jackknife average is the mean of these partial means.
      - The standard error is sqrt((N-1)/N * sum( (x_i - x_jk_avg)^2 )), 
        where x_i are the partial means and x_jk_avg is their average.

    Returns (mean, stderror).
    """
    N = len(data)
    if N < 2:
        return np.mean(data), 0.0

    means = np.zeros(N, dtype=float)
    total_mean = np.mean(data)
    total_sum = np.sum(data)
    # Instead of recomputing the full mean for each iteration, we can do:
    # mean_except_i = (total_sum - data[i]) / (N-1)
    for i in range(N):
        means[i] = (total_sum - data[i]) / (N - 1)

    jk_avg = np.mean(means)
    var = np.sum((means - jk_avg)**2)
    # standard error
    jk_std = np.sqrt((N - 1) / N * var)
    return jk_avg, jk_std

def analyze_sse_data_with_jackknife(input_csv, output_csv):
    """
    Reads SSE results from input_csv, performs a Jackknife 
    for E, C, M, Chi, and writes output to output_csv 
    with columns T, E, dE, C, dC, M, dM, Chi, dChi.
    """

    # We'll store data for each T as a list of samples for E, C, M, Chi
    # However, the SSE CSV likely has only one row per T. That means we
    # DO NOT have multiple samples for the same T in the final aggregated CSV.
    #
    # To use jackknife properly, we need the time-series data for each T 
    # (like a list of energy values from each sweep). That is typically 
    # stored differently, e.g. "run_sse_simulation" logs each measurement 
    # in a separate file. If we only have aggregated means per T, 
    # there's no straightforward way to do a jackknife here.
    #
    # The correct approach: 
    # For each T, you might store the entire distribution of E, M, 
    # then load that distribution for the jackknife. We'll demonstrate 
    # how to handle a single aggregator CSV anyway, but it won't 
    # produce meaningful stats if each T only has one measurement.

    # We'll assume we have multiple lines per T, each line is one sample 
    # from the final phase of the simulation. That's not what run_simulation.py 
    # does currently. So consider this code a template.

    # Let's read everything first
    raw_rows = []
    with open(input_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # We'll assume each row is a sample 
            # (In practice, run_simulation might need to store each sweep's data.)
            raw_rows.append(row)

    # Group rows by T
    data_by_T = {}
    for row in raw_rows:
        T = float(row["T"])
        # E, C, M, Chi
        E = float(row["E"])
        C = float(row["C"])
        M = float(row["M"])
        Chi = float(row["Chi"])

        if T not in data_by_T:
            data_by_T[T] = {
                "E_list": [],
                "C_list": [],
                "M_list": [],
                "Chi_list": []
            }
        data_by_T[T]["E_list"].append(E)
        data_by_T[T]["C_list"].append(C)
        data_by_T[T]["M_list"].append(M)
        data_by_T[T]["Chi_list"].append(Chi)

    # Now compute jackknife for each T
    results = []
    for T, dat in sorted(data_by_T.items(), key=lambda x: x[0]):
        E_arr = np.array(dat["E_list"])
        C_arr = np.array(dat["C_list"])
        M_arr = np.array(dat["M_list"])
        Chi_arr = np.array(dat["Chi_list"])

        E_jk, dE_jk = jackknife(E_arr)
        C_jk, dC_jk = jackknife(C_arr)
        M_jk, dM_jk = jackknife(M_arr)
        Chi_jk, dChi_jk = jackknife(Chi_arr)

        results.append({
            "T": T,
            "E": E_jk, "dE": dE_jk,
            "C": C_jk, "dC": dC_jk,
            "M": M_jk, "dM": dM_jk,
            "Chi": Chi_jk, "dChi": dChi_jk
        })

    # Write to output
    with open(output_csv, "w", newline="") as f:
        fieldnames = ["T", "E", "dE", "C", "dC", "M", "dM", "Chi", "dChi"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"Jackknife analysis written to {output_csv}")

if __name__ == "__main__":
    # Example usage with hypothetical data
    # For it to work properly, the input CSV must contain multiple rows with the 
    # same T, each row representing one sample from the simulation.
    in_csv = os.path.join("data", "results_temperature_sweep_L12.csv")
    out_csv = os.path.join("data", "results_temperature_sweep_L12_jackknife.csv")
    analyze_sse_data_with_jackknife(in_csv, out_csv)
