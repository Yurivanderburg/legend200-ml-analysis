"""
Calculate pulse shape discrimination efficiency for 2vbb events. 
Two windows are used:
1) [1000, 1300] keV
2) [1530, 1635] keV \ (42K gamma lines)
"""



# --- Standard library ---
import numpy as np
import pickle
import os

# --- Project modules ---
from utils.data_io import load_pickle_file


# ------------------------------------------------------------
# Core function
# ------------------------------------------------------------
def psd_eff_2vbb(config: dict) -> dict:
    """ Determine PSD eff for 2vbb events"""

    data = load_pickle_file(config['2vbb']['path_in'])
    final_results = {}

    # Perform again a sanity check - should not do anything but just to be sure !
    data = data[
        data.geds.psd.is_good &
        ~data.coincident.muon & 
        ~data.coincident.muon_offline & 
        ~data.coincident.spms 
    ]

    # Define two windows: Window 1 between 1000 and 1300 keV
    # Define two windows: Window 2 between 1530 and 1635 keV, but excluding gamma lines 42K
    energy_windows = {}
    energy_windows['window1'] = ((data.geds.energy > 1000) & (data.geds.energy < 1300))
    energy_windows['window2'] = ((data.geds.energy > 1530) & 
                                (data.geds.energy < 1725) &
                                ((data.geds.energy < 1583) | (data.geds.energy > 1598)) &
                                ((data.geds.energy < 1616) | (data.geds.energy > 1635)))

    for window, energy_mask in energy_windows.items():
        final_results[window] = {}

        data_windowed = data[energy_mask]

        for model in config['models']:

            final_results[window][model] = {}

            if model == 'AoE':
                eff, eff_err = get_2vbb_eff(
                                T_p = len(data_windowed[data_windowed.geds.psd.is_bb_like == True]), 
                                T_f = len(data_windowed[data_windowed.geds.psd.is_bb_like == False]),
                                l_B = config['2vbb']['l_B'],
                                l_Bp = config['2vbb']['l_Bp'], 
                                T_p_err = np.sqrt(len(data_windowed[data_windowed.geds.psd.is_bb_like == True])),
                                T_f_err = np.sqrt(len(data_windowed[data_windowed.geds.psd.is_bb_like == False])),
                                l_B_err = config['2vbb']['l_B_err'], 
                                l_Bp_err = config['2vbb']['l_Bp_err']
                                )
            else:
                model_no = model[-1]
                eff, eff_err = get_2vbb_eff(
                                T_p = len(data_windowed[data_windowed['Transformer_v' + str(model_no)].Prediction == 'prob_SSE']), 
                                T_f = len(data_windowed[data_windowed['Transformer_v' + str(model_no)].Prediction != 'prob_SSE']),
                                l_B = config['2vbb']['l_B'],
                                l_Bp = config['2vbb']['l_Bp'],
                                T_p_err = np.sqrt(len(data_windowed[data_windowed['Transformer_v' + str(model_no)].Prediction == 'prob_SSE'])),
                                T_f_err = np.sqrt(len(data_windowed[data_windowed['Transformer_v' + str(model_no)].Prediction != 'prob_SSE'])),
                                l_B_err = config['2vbb']['l_B_err'], 
                                l_Bp_err = config['2vbb']['l_Bp_err']
                                )

            final_results[window][model]['eff'] = eff
            final_results[window][model]['eff_err'] = eff_err

    # Store dictionary as pickle file
    file_out = os.path.join(config['global_path_out'], "Results_2vbb.pkl")
    with open(file_out, "wb") as file:
        pickle.dump(final_results, file)

    return final_results


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def get_2vbb_eff(
    T_p: float,
    T_f: float,
    l_B: float,
    l_Bp: float,
    T_p_err: float,
    T_f_err: float,
    l_B_err: float,
    l_Bp_err: float,
):
    """
    SSE survival fraction of the 2vbb signal with error propagation.

    Parameters
    ----------
    T_p, T_f : float
        Counts passing/failing the cut.
    l_B : float
        Fraction of background in the window.
    l_Bp : float
        Fraction of background passing the cut.
    T_p_err, T_f_err, l_B_err, l_Bp_err : float
        Uncertainties of the above.

    Returns
    -------
    eff : float
        Estimated SSE survival fraction.
    eff_err : float
        Propagated uncertainty.
    """
    eff = (T_p / (T_p + T_f) - l_B * l_Bp) / (1.0 - l_B)

    denom = (T_p + T_f) ** 2 * (1.0 - l_B)
    error_term_T_p = T_f / denom
    error_term_T_f = -T_p / denom
    error_term_l_B = (T_p / (T_p + T_f) - (l_Bp * (1.0 - l_B) + l_B * l_Bp)) / ((1.0 - l_B) ** 2)
    error_term_l_Bp = -l_B / (1.0 - l_B)

    eff_err = np.sqrt(
        (error_term_T_p * T_p_err) ** 2
        + (error_term_T_f * T_f_err) ** 2
        + (error_term_l_B * l_B_err) ** 2
        + (error_term_l_Bp * l_Bp_err) ** 2
    )
    return eff, eff_err
