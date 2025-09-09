"""
Calculate pulse shape discrimination efficiency double escape peaks.
"""


# --- Standard library ---
import os
import traceback
import glob

# --- Third party library ---
import numpy as np
import pickle
import tqdm
import awkward as ak
from scipy.optimize import minimize

# --- Project modules ---
from utils.data_io import load_pickle_file
from scipy.optimize import curve_fit
from utils.math import gauss_bkg



# ------------------------------------------------------------
# Core function
# ------------------------------------------------------------
def psd_eff_dep(config: dict, peak: str) -> dict:
    """Calculate the efficiency and its uncertainty for every detector in the Tl-208 DEP. """

    # Determine file paths
    if peak == 'Th228':
        file_in = sorted(glob.glob(config['DEP_Th228']['path_in']))
        file_out = os.path.join(config['global_path_out'], "Results_DEP_Th228.pkl")
    elif peak == 'Co56':
        file_in = config['DEP_Co56']['path_in']
        file_out = os.path.join(config['global_path_out'], "Results_DEP_Co56.pkl")

    data = load_pickle_file(file_in)

    # Create a new dictionary with detector names as keys (must switch from detid to name)
    if peak == 'Th228':
        det_ids = set(data.geds.rawid[:,0])
    else:
        det_ids = set(data.channel_id)
    
    detector_metadata = load_pickle_file(config['path_to_metadata'])

    name_rawid = {}
    for k,v in detector_metadata.items():
        name_rawid[int(v['daq']['rawid'])] = k


    # Loop over all detector ids -> one fit per detector
    final_results = {}
    for det_id in tqdm.tqdm(det_ids):

        try:
            final_results[name_rawid[int(det_id)]] = {}

            if peak == 'Th228':
                data_single_det = data[data.geds.rawid[:,0] == det_id]
            else:
                data_single_det = data[data.channel_id == det_id]

            # Calculate efficiencies for every detector
            try:
                if peak == 'Th228':
                    results = fit_Th228_peak(data_single_det)
                elif peak == 'Co56':
                    results = fit_Co56_peak(data_single_det)


                for model in config['models']:
                    results[model] = {}
                    eff, eff_err = amplitude_efficiency(results[model + "_pass"]['popt'][0], 
                                                np.sqrt(results[model + "_pass"]['pcov'][0,0]), 
                                                results[model + "_fail"]['popt'][0], 
                                                np.sqrt(results[model + "_fail"]['pcov'][0,0]))
                    results[model]['eff'] = eff
                    results[model]['eff_err'] = eff_err
                    results['n_data'] = len(data_single_det)

                final_results[name_rawid[int(det_id)]] = results
            except Exception:
                print("Fit failed for detector {}".format(det_id))
                print(traceback.format_exc())

        except:
            print("Failed for detector {}".format(det_id))

        
    # Create a summary entry -> perform a Chi2 fit with floating variance
    final_results['summary'] = {}
    # For every model, calculate total efficiency
    for model in config['models']:

        final_results['summary'][model] = {}
        effs = []
        eff_errs = []
        for k,v in sorted(final_results.items()):
            if k != 'summary':
                effs.append(v[model]['eff'])
                eff_errs.append(v[model]['eff_err'])

        final_results['summary'][model]['effs'] = np.array(effs)
        final_results['summary'][model]['eff_errs'] = np.array(eff_errs)

        # Minimize the χ² function
        initial_guess = [np.mean(effs), 0.02]
        result = minimize(chi_squared, initial_guess, args=(np.array(effs), np.array(eff_errs)), method='L-BFGS-B', bounds=[(0, 1), (0, None)])
        final_results['summary'][model]['fit_result'] = result

    with open(file_out, "wb") as file:
        pickle.dump(final_results, file)


    return final_results



def psd_eff_timevar(config: dict) -> dict:
    """Calculate the efficiency and its uncertainty for every run in the Tl-208 DEP. """

    data = load_pickle_file(sorted(glob.glob(config['DEP_Th228']['path_in'])))

    # Sort detectors by runs & Loop over all runs
    final_results = {}
    for run in tqdm.tqdm(set(data.run)):

        final_results[run] = {}
        data_single_run = data[data.run == run]

        # Calculate efficiencies for every detector
        try:
            results = fit_Th228_peak(data_single_run)

            for model in config['models']:
                results[model] = {}
                eff, eff_err = amplitude_efficiency(results[model + "_pass"]['popt'][0], 
                                            np.sqrt(results[model + "_pass"]['pcov'][0,0]), 
                                            results[model + "_fail"]['popt'][0], 
                                            np.sqrt(results[model + "_fail"]['pcov'][0,0]))
                results[model]['eff'] = eff
                results[model]['eff_err'] = eff_err
                results['n_data'] = len(data_single_run)

            final_results[run] = results
        except Exception:
            print("Fit failed for run {}".format(run))
            print(traceback.format_exc())

    # Create a summary entry -> perform a Chi2 fit with floating variance
    final_results['summary'] = {}
    # For every model, calculate total efficiency
    for model in config['models']:

        final_results['summary'][model] = {}
        effs = []
        eff_errs = []
        for k,v in sorted(final_results.items()):
            if k != 'summary':
                effs.append(v[model]['eff'])
                eff_errs.append(v[model]['eff_err'])

        final_results['summary'][model]['effs'] = np.array(effs)
        final_results['summary'][model]['eff_errs'] = np.array(eff_errs)

        # Minimize the χ² function
        initial_guess = [np.mean(effs), 0.02]
        result = minimize(chi_squared, initial_guess, args=(np.array(effs), np.array(eff_errs)), method='L-BFGS-B', bounds=[(0, 1), (0, None)])
        final_results['summary'][model]['fit_result'] = result

    file_out = os.path.join(config['global_path_out'], "Results_timevar.pkl")
    with open(file_out, "wb") as file:
        pickle.dump(final_results, file)

    return final_results


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def fit_Th228_peak(data: ak.Array, bin_keV: int = 1):
    """
    Fit the 208Tl DEP (around 1593 keV) for total and PSD selections.

    Parameters
    ----------
    data : ak.Array
        Input dataset with fields like `geds.energy`, `geds.psd.is_bb_like`,
        and `Transformer_v{1,2,3}.Prediction`.
    bin_keV : int, optional
        Histogram bin width in keV, by default 1.

    Returns
    -------
    dict
        Mapping of selection name -> dict with 'energies', 'counts', 'popt', 'pcov'.
    """
    results = {}

    # Define binning and total counts
    binning = np.arange(1570, 1615 + bin_keV, bin_keV)
    counts_total, bin_edges = np.histogram(ak.to_numpy(data.geds.energy[:, 0]), bins=binning)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    # total fit (use tight bounds around 1593 keV)
    bkg = float(np.mean(counts_total[:10]))
    amp = float(np.max(counts_total) - bkg)
    stp = float(np.mean(counts_total[:10]) - np.mean(counts_total[-10:]))
    p0 = [amp, 1593.0, 1.0, 0.0, bkg, stp]
    bounds = ([0.0, 1590.0, 0.0, -1e-4, 0.0, -1000.0], [np.inf, 1596.0, 5.0, 1e-4, 1e4, 1e5])

    popt_total, pcov_total = curve_fit(gauss_bkg, bin_centers, counts_total, p0=p0, bounds=bounds)

    results["Total"] = {
        "energies": ak.to_numpy(data.geds.energy[:, 0]),
        "counts": counts_total,
        "popt": popt_total,
        "pcov": pcov_total,
    }

    # Automate selections
    selections = [
        ("AoE_pass", data.geds.psd.is_bb_like[:, 0] == True),
        ("AoE_fail", data.geds.psd.is_bb_like[:, 0] == False),
        ("T1_pass", data.Transformer_v1.Prediction == "prob_SSE"),
        ("T1_fail", data.Transformer_v1.Prediction != "prob_SSE"),
        ("T2_pass", data.Transformer_v2.Prediction == "prob_SSE"),
        ("T2_fail", data.Transformer_v2.Prediction != "prob_SSE"),
        ("T3_pass", data.Transformer_v3.Prediction == "prob_SSE"),
        ("T3_fail", data.Transformer_v3.Prediction != "prob_SSE"),
    ]

    for name, condition in selections:
        energies = data[condition].geds.energy[:, 0]
        counts, popt, pcov = peak_fit(
            energies=energies,
            binning=binning,
            bin_centers=bin_centers,
            peak_energy=1593.0,
        )
        results[name] = {
            "energies": ak.to_numpy(energies),
            "counts": counts,
            "popt": popt,
            "pcov": pcov,
        }

    return results


def fit_Co56_peak(data: ak.Array, bin_keV: int = 1):
    """
    Fit the Co-56 DEP (around 2231.5 keV) for total and PSD selections.

    Parameters
    ----------
    data : ak.Array
        Input dataset with fields like `trapEmax_ctc_cal`, `AoE_Double_Sided_Cut`,
        and `Transformer_v{1,2,3}.Prediction`.
    bin_keV : int, optional
        Histogram bin width in keV, by default 1.

    Returns
    -------
    dict
        Mapping of selection name -> dict with 'energies', 'counts', 'popt', 'pcov'.
    """
    results = {}

    binning = np.arange(2215, 2247 + bin_keV, bin_keV)
    counts_total, bin_edges = np.histogram(ak.to_numpy(data.trapEmax_ctc_cal), bins=binning)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    bkg = float(np.mean(counts_total[:10]))
    amp = float(np.max(counts_total) - bkg)
    stp = float(np.mean(counts_total[:10]) - np.mean(counts_total[-10:]))
    p0 = [amp, 2231.5, 1.0, 0.0, bkg, stp]
    bounds = ([0.0, 2229.0, 0.0, -1e-4, 0.0, -1000.0], [np.inf, 2235.0, 5.0, 1e-4, 1e4, 1e5])

    popt_total, pcov_total = curve_fit(gauss_bkg, bin_centers, counts_total, p0=p0, bounds=bounds)

    results["Total"] = {
        "energies": ak.to_numpy(data.trapEmax_ctc_cal),
        "counts": counts_total,
        "popt": popt_total,
        "pcov": pcov_total,
    }

    selections = [
        ("AoE_pass", data.AoE_Double_Sided_Cut == True),
        ("AoE_fail", data.AoE_Double_Sided_Cut == False),
        ("T1_pass", data.Transformer_v1.Prediction == "prob_SSE"),
        ("T1_fail", data.Transformer_v1.Prediction != "prob_SSE"),
        ("T2_pass", data.Transformer_v2.Prediction == "prob_SSE"),
        ("T2_fail", data.Transformer_v2.Prediction != "prob_SSE"),
        ("T3_pass", data.Transformer_v3.Prediction == "prob_SSE"),
        ("T3_fail", data.Transformer_v3.Prediction != "prob_SSE"),
    ]

    for name, condition in selections:
        energies = data[condition].trapEmax_ctc_cal
        counts, popt, pcov = peak_fit(
            energies=energies,
            binning=binning,
            bin_centers=bin_centers,
            peak_energy=2231.5,
        )
        results[name] = {
            "energies": ak.to_numpy(energies),
            "timestamp": ak.to_numpy(data[condition].timestamp),
            "counts": counts,
            "popt": popt,
            "pcov": pcov,
        }

    return results


def amplitude_efficiency(A_surv: float, A_surv_err: float, A_fail: float, A_fail_err: float):
    """
    Compute efficiency from surviving and failing amplitudes (background-reduced).

    Parameters
    ----------
    A_surv, A_fail : float
        Amplitudes for events passing/failing the PSD selection.
    A_surv_err, A_fail_err : float
        Uncertainties of the amplitudes.

    Returns
    -------
    eff : float
        Efficiency A_surv / (A_surv + A_fail).
    eff_err : float
        Propagated uncertainty.
    """
    eff = A_surv / (A_surv + A_fail)
    eff_err = eff * (1.0 - eff) * np.sqrt((A_surv_err / A_surv) ** 2 + (A_fail_err / A_fail) ** 2)
    return eff, eff_err


def chi_squared(params: list | np.ndarray, eps_i: np.ndarray, 
                sigma_eps_i: np.ndarray) -> float:
    """
    Compute the chi-squared function with floating mean and systematic width.

    Parameters
    ----------
    params : list or array
        [mean_efficiency, delta_stab_p]
    eps_i : np.ndarray
        Array of efficiency measurements eps_i
    sigma_eps_i : np.ndarray
        Array of statistical uncertainties sig_eps,i
    
    Returns
    -------
    float
        The χ² value
    """
    mean_eff, delta_stab_p = params

    if delta_stab_p < 0:
        return np.inf  # avoid non-physical negative widths

    total_var = sigma_eps_i**2 + delta_stab_p**2
    chi2_terms = ((eps_i - mean_eff)**2) / total_var + np.log(total_var)
    return np.sum(chi2_terms)



def peak_fit(
    energies: np.ndarray,
    binning: np.ndarray,
    bin_centers: np.ndarray,
    peak_energy: float,
):
    """
    Fit Gaussian + linear background + step function to a peak histogram.

    Parameters
    ----------
    energies : np.ndarray
        Event energy array to be histogrammed.
    binning : np.ndarray
        Bin edges for the histogram.
    bin_centers : np.ndarray
        Precomputed centers of the bins (same length as counts).
    peak_energy : float
        Expected peak position (keV) used to initialize the fit.

    Returns
    -------
    counts : np.ndarray
        Histogram counts corresponding to `binning`.
    popt : np.ndarray
        Best-fit parameters [A, mu, sig, a, b, d].
    pcov : np.ndarray
        Covariance matrix of the fit parameters (6x6).
    """
    counts, _ = np.histogram(energies, bins=binning)

    # crude initial guesses
    bkg = np.mean(counts[:10])
    amp = np.max(counts) - bkg
    stp = np.mean(counts[:10]) - np.mean(counts[-10:])

    p0 = [amp, peak_energy, 1.0, 0.0, bkg, stp]
    bounds = (
        [0.0, peak_energy - 3.0, 0.0, -1e-6, 0.0, -np.inf],
        [np.inf, peak_energy + 3.0, 5.0, 1e-6, np.inf, np.inf],
    )

    popt, pcov = curve_fit(gauss_bkg, bin_centers, counts, p0=p0, bounds=bounds)
    return counts, popt, pcov