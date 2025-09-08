"""
Utils for determining the Pulse Shape Discrimination efficiency at 
$Q_{\beta \beta}$ for LEGEND-200 ML Analysis.
"""

# --- Standard library ---
import os
import pickle

# --- Third-party ---
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

# --- Project modules ---
from utils.math import linear, linear_with_err



# ------------------------------------------------------------
# Core function
# ------------------------------------------------------------
def psd_eff_qbb(
    config: dict,
    efficiencies_Th228: dict,
    efficiencies_2vbb: dict,
    efficiencies_Co56: dict,
    efficiencies_Th228_timevar: dict,
):
    """
    Combine Th-228 DEP, 2vbb window, and Co-56 DEP into final PSD efficiency.
    """
    data_x = np.array([1593.5, 1150.0, 2231.5])
    efficiency = {}

    for model in config['models']:
        efficiency[model] = {}

        mean_th228, sig_th228 = mean_and_error_with_delta(
            eps=efficiencies_Th228["summary"][model]["effs"],
            sig=efficiencies_Th228["summary"][model]["eff_errs"],
            delta=efficiencies_Th228["summary"][model]["fit_result"].x[1],
        )
        mean_co56, sig_co56 = mean_and_error_with_delta(
            eps=efficiencies_Co56["summary"][model]["effs"],
            sig=efficiencies_Co56["summary"][model]["eff_errs"],
            delta=efficiencies_Co56["summary"][model]["fit_result"].x[1],
        )

        data_y = np.array(
            [
                mean_th228,
                float(efficiencies_2vbb["window1"][model]["eff"]),
                mean_co56,
            ]
        )
        sigma_y = np.array(
            [
                sig_th228,
                float(efficiencies_2vbb["window1"][model]["eff_err"]),
                sig_co56,
            ]
        )

        popt, pcov = curve_fit(
            linear,
            data_x,
            data_y,
            sigma=sigma_y,
            absolute_sigma=True,
            p0=[-1e-4, 1.0],
            bounds=([-np.inf, 0.0], [0.0, 1.0]),
        )
        efficiency[model]["fit"] = popt

        eff_qbb, eff_qbb_err = linear_with_err(x=2039.0, fit_param=popt, covariance_matrix=pcov)

        # Save new dictionary with results
        efficiency[model]["values"] = {"qbb": eff_qbb, "2vbb_diff": -0.007}
        efficiency[model]["errors"] = {
            "qbb": eff_qbb_err,
            "timevar": efficiencies_Th228_timevar["summary"][model]["fit_result"].x[1],
            "noise": 0.0065,
            "2vbb_diff": 0.02,
        }

        efficiency[model]["efficiency"] = efficiency[model]["values"]["qbb"] + efficiency[model]["values"]["2vbb_diff"]
        efficiency[model]["uncertainty"] = np.sqrt(
            efficiency[model]["errors"]["qbb"] ** 2
            + efficiency[model]["errors"]["timevar"] ** 2
            + efficiency[model]["errors"]["noise"] ** 2
            + efficiency[model]["errors"]["2vbb_diff"] ** 2
        )
        
        # Also store Th228 and Co56 and 2vbb individual results
        efficiency[model]['Th228'] = {}
        efficiency[model]['Th228']['values'] = efficiencies_Th228['summary'][model]
        efficiency[model]['Th228']['efficiency'] = mean_th228
        efficiency[model]['Th228']['uncertainty'] = sig_th228


        efficiency[model]['Co56'] = {}
        efficiency[model]['Co56']['values'] = efficiencies_Th228['summary'][model]
        efficiency[model]['Co56']['efficiency'] = mean_co56
        efficiency[model]['Co56']['uncertainty'] = sig_co56

        efficiency[model]['2vbb'] = efficiencies_2vbb["window1"][model]

    # Store results
    file_out = os.path.join(config['global_path_out'], "Results_combined.pkl")
    with open(file_out, "wb") as file:
        pickle.dump(efficiency, file)

    return efficiency, data_y, sigma_y


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def get_uncertainty_band(x: np.ndarray, popt: np.ndarray, pcov: np.ndarray) -> np.ndarray:
    """
    Propagate fit parameter covariance to model uncertainty for gauss_bkg.

    Parameters
    ----------
    x : np.ndarray
        Points where the band should be evaluated.
    popt : np.ndarray
        Best-fit parameters [A, mu, sig, a, b, d].
    pcov : np.ndarray
        6x6 covariance matrix of the parameters.

    Returns
    -------
    y_err : np.ndarray
        One-sigma uncertainty at each x (same shape as x).

    Notes
    -----
    This implements standard linear error propagation:
        sig_y^2 = J · Cov · J^T
    with J = ∂f/∂θ evaluated at the best-fit parameters.
    """
    A, mu, sig, a, b, d = popt
    x = np.asarray(x, dtype=float)

    # Precompute terms
    z = (x - mu) / (np.sqrt(2.0) * sig)
    exp_term = np.exp(-((x - mu) ** 2) / (2.0 * sig**2))
    erfc_term = erfc(z)
    erf_gauss = (2.0 / np.sqrt(np.pi)) * np.exp(-z**2)  # derivative of erfc

    # Partials of f(x) = A*exp(...) + a*x + b + 0.5*d*erfc(z)
    dA = exp_term
    dmu = A * exp_term * (x - mu) / (sig**2) + 0.5 * d * erf_gauss * (1.0 / (np.sqrt(2.0) * sig))
    dsig = A * exp_term * ((x - mu) ** 2) / (sig**3) + 0.5 * d * erf_gauss * (z / sig)
    da = x
    db = np.ones_like(x)
    dd = 0.5 * erfc_term

    # Stack gradients -> shape (N, 6)
    J = np.stack([dA, dmu, dsig, da, db, dd], axis=-1)

    # σ_y^2 = J * Cov * J^T  -> compute efficiently with einsum
    y_var = np.einsum("ni,ij,nj->n", J, pcov, J, optimize=True)
    y_var = np.maximum(y_var, 0.0)  # numeric safety
    return np.sqrt(y_var)




def mean_and_error_with_delta(eps, sig, delta):
    """
    Weighted mean with per-point sigma and an additional common systematic delta.

    Parameters
    ----------
    eps : array_like
        Measurements.
    sig : array_like
        Statistical uncertainties for each measurement.
    delta : float
        Additional (common) systematic added in quadrature.

    Returns
    -------
    mu : float
        Weighted mean.
    sigma_mu : float
        Uncertainty of the weighted mean.
    """
    eps = np.asarray(eps, float)
    sig = np.asarray(sig, float)
    v = sig**2 + float(delta) ** 2
    w = 1.0 / v
    mu = np.sum(w * eps) / np.sum(w)
    sigma_mu = 1.0 / np.sqrt(np.sum(w))
    return mu, sigma_mu




