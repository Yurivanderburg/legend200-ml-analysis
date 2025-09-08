"""
Summarize LEGEND sensitivity results from JSON files.

Usage
-----
python -m scripts.summarize_results --path path/to/results --output summary.json
"""

import os
import glob
import json
import argparse
import numpy as np
from tqdm import tqdm

Q_bb = 2039.0  # keV


def extract_from_file(res: dict, evts: dict) -> dict:
    """
    Extract key quantities from a single JSON result.

    Parameters
    ----------
    res: dict
        Single results dictionary
    evts: dict
        Single fake events dictionary

    Returns
    -------
    dict
        Dictionary with scalar values for this file.
    """
    # --- Global mode S and eff ---
    s_val = res["global_modes"]["S"] * 1e-27
    eff_val = res["global_modes"]["ε"]

    # --- CI on S ---
    ci_left = [ele["left"] for ele in res["ci_68"]["S"]]
    ci_right = [ele["right"] for ele in res["ci_68"]["S"]]

    # --- Energy resolution, ROI counts ---
    omega = res["global_modes"]["ω"][0]
    fwhm = 2.355 * omega
    roi_low, roi_high = Q_bb - fwhm, Q_bb + fwhm
    energies = [e["energy"] for e in evts["events"]]
    roi_count = sum(roi_low <= e <= roi_high for e in energies)

    # --- Efficiency CI ---
    eff_ci_eps = res["ci_68"]["ε"][0][0]

    # --- Quantiles ---
    S90 = res["quantile90"]["S"] * 1e-27
    b_val = res["refined_global_modes"]["B_l200a_all"]

    return {
        "global_mode_S": s_val,
        "global_mode_eff": eff_val,
        "ci_68_S_left": ci_left,
        "ci_68_S_right": ci_right,
        "ci_68_eff_left": eff_ci_eps["left"],
        "ci_68_eff_right": eff_ci_eps["right"],
        "roi_count": roi_count,
        "S90": S90,
        "global_mode_B": b_val,
    }


def summarize_results(path: str) -> dict:
    """
    Store key quantities in a new dictionary and store as JSON

    Parameters
    ----------
    path : str
        Base path to results JSON files

    Returns
    -------
    results_arr : dict
        Results dictionary (is also stored) 
    """
    mcmc_files = sorted(glob.glob(os.path.join(path, "sensitivity/mcmc_files/*.json")))
    event_files = sorted(glob.glob(os.path.join(path, "sensitivity/fake_data/*.json"))) 
    if not mcmc_files:
        raise FileNotFoundError(f"MCMC files not found.")
    elif not event_files:
        raise FileNotFoundError("Fake data files not found.")

    # Load and extract
    all_results = [
        extract_from_file(json.load(open(f)), json.load(open(f2)))
        for f, f2 in tqdm(zip(mcmc_files, event_files))
    ]

    # Convert to plain lists instead of numpy arrays
    results_arr = {
        k: [res[k] for res in all_results]
        for k in all_results[0].keys()
    }
    # Add summary values
    summary = {
        "global_mode_S": {
            "median": float(np.median(results_arr["global_mode_S"])),
            "mean": float(np.mean(results_arr["global_mode_S"])),
        },
        "global_mode_eff": {
            "median": float(np.median(results_arr["global_mode_eff"])),
            "mean": float(np.mean(results_arr["global_mode_eff"])),
        },        
        "S_68_CI": {
            "median_left": float(np.median(np.concatenate(results_arr["ci_68_S_left"]))),
            "median_right": float(np.median(np.concatenate(results_arr["ci_68_S_right"]))),
        },
        "ROI_counts": {
            "median": int(np.median(results_arr["roi_count"])),
            "mean": float(np.mean(results_arr["roi_count"])),
        },
        "eff_68_CI": {
            "median_left": float(np.median(results_arr["ci_68_eff_left"])),
            "median_right": float(np.median(results_arr["ci_68_eff_left"])),
        },
        "S90": {
            "low": float(np.percentile(results_arr["S90"], 16)),
            "med": float(np.percentile(results_arr["S90"], 50)),
            "high": float(np.percentile(results_arr["S90"], 84)),
        },
        "T90": {
            "low": float(1 / np.percentile(results_arr["S90"], 84)),
            "med": float(1 / np.percentile(results_arr["S90"], 50)),
            "high": float(1 / np.percentile(results_arr["S90"], 16)),
        },
        "global_mode_B": {
            "median": float(np.median(results_arr["global_mode_B"]))
        }
    }

    results_arr['summary'] = summary

    return results_arr

