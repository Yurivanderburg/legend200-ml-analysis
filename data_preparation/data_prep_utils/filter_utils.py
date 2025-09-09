"""
Utility functions for the filtering datasets with filter_dataset.py

"""

# Standard library
import os
import sys
import argparse
import logging
import datetime as dt

# Third-party
import numpy as np
import awkward as ak
import pandas as pd
from legendmeta import LegendMetadata
from scipy.signal import find_peaks

# Project modules




def attach_labels(data: ak.Array, ortec_type: bool, config: dict)-> ak.Array:
    """ 
    Attaches labels to the data. This is achieved by looping calling create_label_mask()
    and looping over the resulting dictionary. Finally merges the datasets and returns 
    one awkward array with all labelled events.

    Parameters
    ----------
    data: ak.Array
        awkward array of where labels should be attached
    ortec_type: bool
        Whether the detectord was produced by Ortec or not

    Returns
    -------
    merged_data: ak.Array
        awkward array with the labels attached
    """ 
    masks = create_label_mask(data=data, ortec_type=ortec_type, config=config)

    data_list = []
    for label, mask in masks.items():

        data_subset = data[mask]
        data_subset['label'] = label
        data_list.append(data_subset)
    
    merged_data = ak.concatenate(data_list)          

    return merged_data



def create_label_mask(data: ak.Array, ortec_type: bool, config: dict) -> dict:
    """ 
    Creates one label mask for the four different event types considered (MSE, SSE, pcontact, 
    ncontact). Also distinguishes different detector types and manufacturers. If specified in 
    the config, it will also apply an energy cut on the n-contact events (because they are 
    unreliable for low energies - so they will only be considered above a certain energy value.)

    Parameters
    ----------
    data: ak.Array
        awkward array with the events
    ortec_type: bool
        Whether the detectord was produced by Ortec or not
    config: dict
        Config file as parsed

    Returns
    -------
    mask_dict: dict
        dictionary of the four masks, labels as keys
    """

    if ortec_type == True:

        mask_SSE = (
            data.AoE_Double_Sided_Cut & data.LQ_Cut
        )

        mask_MSE = (
            ~data.AoE_Low_Cut & ~data.is_surface #& data.AoE_High_Side_Cut 
        )

        mask_ncontact = ( # determine on low-energy peak
            data.AoE_High_Side_Cut & ~data.AoE_Low_Cut & ~data.LQ_Cut & data.is_surface
        )

        mask_pcontact = (
            data.AoE_Low_Cut & ~data.LQ_Cut
        )
    else:
        mask_SSE = (
            data.AoE_Double_Sided_Cut
        )

        mask_MSE = (
            ~data.AoE_Low_Cut & ~data.is_surface # & data.AoE_High_Side_Cut 
        )

        mask_ncontact = ( # determine on low-energy peak
            data.AoE_High_Side_Cut & ~data.AoE_Low_Cut & data.is_surface
        )

        mask_pcontact = (
            ~data.AoE_High_Side_Cut
        )

    if config['FILTER']['energy_cut_on_nconact'] == True:
        mask_ncontact_energycut = (data.peak == 'DEP_Tl')
        mask_ncontact = mask_ncontact & mask_ncontact_energycut


    mask_dict = {'SSE': mask_SSE, 'MSE': mask_MSE, 'ncontact': mask_ncontact, 'pcontact': mask_pcontact}

    return mask_dict



def peak_finding(data: ak.Array, logger: logging.Logger, config: dict) -> ak.Array:
    """
    Peak-finding algorithm. Performs a peak-finding on the smeared derivatives of every event.
    Calls weighted_moving_average() to smear the signal. Several scipy.find_peaks() are stored
    as new fields, most importantly the FWHM of the largest peak and the number of peaks. 

    Parameters
    ----------
    data; ak.Array
        Labelled array containing events 
    logger: logging.Logger
        Logger
    config: dict
        Config file as parsed

    Returns
    -------
    data_out: ak.Array
        Array with new fields
    """
    # Initialize new fields
    peak_val = []
    scp_npeaks = []
    scp_peaklist = []
    scp_width = [] 

    # Loop over all events
    for i, ele in enumerate(data):
        waveform = ele.waveform
        wf_smeared = weighted_moving_average(waveform, half_win_samps=2, gauss_sigma_samps=1.5)
        grad_smeared = np.gradient(wf_smeared)
        grad_smeared_norm = grad_smeared/np.max(grad_smeared)

        # Find-peaks function
        peaks, properties = find_peaks(
                    x = grad_smeared_norm,
                    prominence = config['FILTER']['find_peaks_SSE']['prominence'],
                    distance = config['FILTER']['find_peaks_SSE']['distance'], 
                    width = config['FILTER']['find_peaks_SSE']['width'], 
                    height = config['FILTER']['find_peaks_SSE']['height']
                    )
                        
        # Attach peak & width of highest peak
        if len(peaks) != 0:
            scp_peaklist.append(list(peaks))
            scp_width.append(properties['widths'][np.argmax(properties['prominences'])])
        else:
            scp_peaklist.append([])
            scp_width.append(0)
        
        # Attach n-peaks and baseline list
        scp_npeaks.append(int(len(peaks)))
        
        # Attach peak-validation label
        if ele.label == 'MSE':
            if len(peaks) > 1:
                peak_val.append(bool(True))
            else:
                peak_val.append(bool(False))
        else:
            if len(peaks) == 1:
                peak_val.append(bool(True))
            else: 
                peak_val.append(bool(False))

    # Attach new fields
    data['scp_peaklist'] = ak.enforce_type(scp_peaklist, 'var * int64')
    data['peak_val'] = ak.enforce_type(peak_val, 'bool')
    data['scp_npeaks'] = ak.enforce_type(scp_npeaks, 'int32')
    data['scp_widths'] = ak.enforce_type(scp_width, 'float32')

    return data



def remove_outliers(data: ak.Array, logger: logging.Logger, filename_meta: str, 
                    config: dict) -> ak.Array:
    """
    Remove all events outside a certain standard deviation (specified in the config) of 
    a field. This is achieved by first looping over the detectors and fields, calculating
    mean and standard deviation. In a second loop, every event is then evaluated for all 
    outlier_cuts specified in the config, and the masks are stored in a dictionary. For 
    reference, this dictionary is stored as a pickle file. Lastly, the new fields are 
    attached to the data. 
    NOTE: There for sure is a more elegant way to do this.

    Parameters
    ----------
    data: ak.Array
        Labelled akward array containing events
    logger: logging.Logger
        Logger to be able to log infos
    filename_meta: str
        Name of the pickle output file.  
    config: dict
        Config file as parsed

    Returns
    -------
    data: ak.Array
        Array containing new fields
    """

    # Parameters
    dict_calc = {}
    mask_dict = {}

    # Use smaller variable names for readabilty
    outlier_cuts = config['FILTER']['outlier_cuts']
    n_sig = config['FILTER']['cut_value_outliers']
    n_sig_width = config['FILTER']['cut_value_width']

    for cut in outlier_cuts:
        mask_dict[cut] = []

    # Calculate new risetime fields 
    data['rt1'] = data.tp_90 - data.tp_10
    data['rt2'] = data.tp_50 - data.tp_10  

    # Create dictionary with mean and std of every attribute
    for ch_id in set(data.channel_id):
        dict_calc[str(ch_id)] = {}
        for label in set(data.label):
            dict_calc[str(ch_id)][label] = {}
            dict_calc[str(ch_id)][label]['n'] = len(data[(data.channel_id == ch_id) & (data.label == label)])

            for cut in outlier_cuts:
                dict_calc[str(ch_id)][label][f"{cut}_mean"] = np.mean(data[(data.channel_id == ch_id) & 
                                                                       (data.label == label)][cut])
                dict_calc[str(ch_id)][label][f"{cut}_std"] = np.std(data[(data.channel_id == ch_id) & 
                                                                     (data.label == label)][cut])
                                                                    

    # Loop over all the events to determine if it is an outlier
    for event in data:
        ch_id = str(event.channel_id)

        for cut in outlier_cuts:

            if cut == 'scp_widths':
                # We only want to exclude MSE events that are too narrow
                if event.label == 'MSE':
                    mask = (event[cut] > dict_calc[ch_id][event.label][f"{cut}_mean"]
                        - n_sig_width * dict_calc[ch_id][event.label][f"{cut}_std"])
                # For all other labels, we only want to exclude events that are too wide (merged pulses)
                else:
                    mask = (event[cut] < dict_calc[ch_id][event.label][f"{cut}_mean"]
                        + n_sig_width * dict_calc[ch_id][event.label][f"{cut}_std"])
            
            else:
    
                mask = ((event[cut] < dict_calc[ch_id][event.label][f"{cut}_mean"]
                            + n_sig * dict_calc[ch_id][event.label][f"{cut}_std"]) & 
                        (event[cut] > dict_calc[ch_id][event.label][f"{cut}_mean"]
                            - n_sig * dict_calc[ch_id][event.label][f"{cut}_std"]))
                

            mask_dict[cut] += [mask]

    # Add new column
    for cut in outlier_cuts:
        data[f"{cut}_mask"] = mask_dict[cut]    

    
    # Remove outliers if not debugging
    if config['debug_mode'] == False:
        n_start = len(data)
        for cut in outlier_cuts:
            data = data[data[f"{cut}_mask"] == True]
        n_end = len(data)
        percent_removed = np.round(100*(n_start-n_end)/(n_start), 2)
        logger.info(f"Outlier cuts removed {n_start-n_end} ({percent_removed}%) events.")
    else:
        logger.info(f"Outlier cuts are not removed.")

    return data



def filter_metadata_short(legend_metadata: dict) -> dict:
    """
    Load Legend Metadata and select only valid_detector_types. Returns a dictionary
    which has the channel id's as keys. Very convinent to quickly grab the most important
    metadata based on a channel id. 

    Parameters
    ----------
    legend_metadata: dict
        Instance LegendMetadata which contains LEGEND metadata

    Returns
    -------
    detectors: dict
        Dict containing name and id of the expected detectors
    """
    detectors = {}

    # Loop over legend metadata dictionary
    for key, item in legend_metadata.channelmap().items():

        # Skip non-geds detector systems
        if not item['system'] == 'geds':
            continue
        detector_dict = {}

        # Add important fields
        detector_dict['name'] = key
        detector_dict['type'] = item['type']
        detector_dict['manufacturer'] = item['production']['manufacturer']
        detector_dict['processable'] = item['analysis']['processable']
        detector_dict['lq'] = item['analysis']['psd']['status']['lq']
        detectors[f"ch{item['daq']['rawid']}"] = detector_dict 

    return detectors



def weighted_moving_average(arr, half_win_samps:int, gauss_sigma_samps:float, center:bool=True):
    '''
    This function produce a smoothed array of the input array (the "arr" input parameter)
    with a rolling window and a gaussian kernel. This function uses the functionalities 
    provided by Pandas.

    Parameters
    ----------
    arr
        The original array of amplitudes that should be smoothed. 
    half_win_samps: int
        The half_width of the window in samples units around the central sample
    Gauss_sigma_samps: float
        Width of the gaussian used as smoothing kernel (convolution kernel)
    center: bool
        If True the rolling weighted mean is computed at middle sample of the sliding window

    Returns
    -------
    output: np.ndarray
        The smoothed array of amplitudes with the same shape as the one at the input. 
                           It returns None if the input array (or sequence) is not convertible to a 1D array.
    '''
    arr = np.array(arr)
    if(arr.ndim!=1):
        print('\nERROR --> weighted_moving_average: The input 1D array parameter has a wrong shape:  \
              {}!'.format(arr.shape), file=sys.stderr)
        return None
    
    ser = pd.Series(arr)
    
    output = ser.rolling(window=2*half_win_samps+1, center=center, closed='both', min_periods=half_win_samps+1, win_type='gaussian').mean(std=gauss_sigma_samps )
    return output
