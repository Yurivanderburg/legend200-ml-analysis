"""
Utility functions for the data preparation script prepare_dataset.py

"""

# Standard library
import os
import sys
import argparse
import logging
import datetime as dt
import time
import glob

# Third-party
import numpy as np
import awkward as ak
import pandas as pd
from legendmeta import LegendMetadata
from scipy.signal import find_peaks
import lgdo.lh5 as lh5
from pygama.pargen.data_cleaning import get_tcm_pulser_ids



def load_filenames(period: str, run: str, logger: logging.Logger, config: dict) -> dict:
    """
    Collect input file paths for a given period and run.

    Searches RAW and all configured TIER subdirectories for `.lh5` files,
    and optionally a PAR `.json` file if energy_mode == 'peaks'. 
    Logs errors if expected files are missing.

    Parameters
    ----------
    period : str
        Period identifier (e.g. 'p03').
    run : str
        Run identifier (e.g. 'r001').
    logger : logging.Logger
        Logger
    config: dict
        Config file as parsed

    Returns
    -------
    filenames_dict : dict
        Dictionary mapping tier names ('RAW', 'PHT', 'PSP', etc.) to lists of file paths.
        Includes 'PAR' if energy_mode == 'peaks'.
    validated_filenames : bool
        True if all required files are found, False otherwise.
    """
    filenames_dict = {}
    validated_filenames = True

    # Collect RAW tier files
    filenames_dict['RAW'] = sorted(
        glob.glob(os.path.join(
            config['DIRS']['LEGEND_RAW_DIR'],
            config['DIRS']['TYPE'], period, run, '*.lh5'
        ))
    )

    # Collect files for each additional tier
    for tier, tier_subdir in config['DIRS']['TIER_SUBDIRS'].items():
        path = os.path.join(
            config['DIRS']['LEGEND_DATA_DIR'],
            tier_subdir,
            config['DIRS']['TYPE'], period, run, '*.lh5'
        )
        filenames_dict[tier] = sorted(glob.glob(path))
    
    # Validate: all required tiers must have files
    for tier, tier_files in filenames_dict.items():
        if not tier_files:
            logger.error(f"Missing {tier} files for {period} {run}. Skipping run.")
            validated_filenames = False

    # For 'peaks' mode: ensure exactly one PAR file exists
    if config['energy_mode'] == 'peaks':
        filenames_dict['PAR'] = glob.glob(os.path.join(
            config['DIRS']['LEGEND_DATA_DIR'],
            config['DIRS']['PAR_SUBDIR'],
            config['DIRS']['TYPE'], period, run, '*.json'
        ))
        if len(filenames_dict['PAR']) != 1:
            logger.error(f"Missing or multiple PAR files for {period} {run}. Skipping run.")
            validated_filenames = False

    return filenames_dict, validated_filenames



def load_data(pht_files: list, raw_files: list, tcm_files: list, psp_files: list, 
              par_data: dict, channel_id: str, logger: logging.Logger, config: dict) -> ak.Array:
    """ 
    Loads desired waveforms and higher-tier information stored in RAW/PHT files
    into an awkward array. 

    First loads the HIT event (pht_files), then pulser events (tcm_files). Selects 
    the correct energy-calibration parameter (par_data), before calling 
    load_peak_mask() to apply energy cuts (as defined in the config).
    With the correct mask ids, only load the necessary raw events. 
    Additionaly load psp data for DCR parameters. Also calls get_peak_region().
    Finally merge all attributes into one awkward array, and return that array. 

    Parameters
    ----------
    pht_files: list
        Filenames that contain higher-level hit events
    raw_files: list
        Filename that contain raw waveforms
    tcm_files:list
        Filenames to identify pulser events
    psp_files: dict
        Filenames that contain various DSP parameters
    par_data: dict
        Contains run information for energy resolutions
    channel_id: str
        Channel id of the detector used. Will be added as label
    logger: logging.Logger
        Logger to be able to log infos
    config: dict
        Config file as parsed

    Returns
    -------
    merged_data: ak.Array
        Array containg all pht parameters and waveforms   
    """
    # Create LH5Store instance
    store = lh5.LH5Store()

    # Load PHT data (HIT events)
    pht_data, n_pht = store.read(f"ch{channel_id}/hit/", pht_files)
    pht_data = pht_data.view_as("ak")
    logger.info(f"{np.round(n_pht/1e3)}k events found before cuts.")

    # Identify pulser events (TCM). If no TCM files, default to all False.
    if tcm_files:
        _, mask_pulser = get_tcm_pulser_ids(tcm_files, f"{channel_id}", 5)
    else:
        mask_pulser = ak.Array([False] * len(pht_data))

    # Build event mask depending on energy selection mode
    if config["energy_mode"] == "peaks":
        # Use PAR file peaks for selection
        relevant_par_data = par_data[f"ch{channel_id}"]
        total_mask = load_peak_mask(
            pht_data, relevant_par_data, mask_pulser, logger=logger, config=config
        )
        ids = np.where(total_mask)[0].to_numpy()

    elif config["energy_mode"] == "manual_peaks":
        # Manually select events around specified peaks
        total_peak_mask = ak.Array([False] * len(pht_data))
        qc_mask = quality_cut_mask(pht_data)

        for i, peak in enumerate(config["manual_peaks"].values()):
            peak_mask = (
                (pht_data.trapEmax_ctc_cal > float(peak) - config["manual_peak_region"])
                & (pht_data.trapEmax_ctc_cal < float(peak) + config["manual_peak_region"])
            )
            total_peak_mask = peak_mask if i == 0 else (total_peak_mask | peak_mask)

        total_mask = total_peak_mask & ~mask_pulser & qc_mask
        ids = np.where(total_mask)[0].to_numpy()

    else:
        # Select by continuous energy region
        energy_mask = (
            (pht_data.trapEmax_ctc_cal > float(config["continuous_region"][0]))
            & (pht_data.trapEmax_ctc_cal < float(config["continuous_region"][1]))
        )
        qc_mask = quality_cut_mask(pht_data)
        total_mask = energy_mask & ~mask_pulser & qc_mask
        ids = np.where(total_mask)[0].to_numpy()

    # Load RAW data (waveforms and baseline) for selected events
    time_ = time.time()
    raw_data, _ = store.read(
        f"ch{channel_id}/raw/",
        raw_files,
        idx=ids,
        field_mask=["waveform_windowed", "baseline", "timestamp"],
    )
    raw_data = raw_data.view_as("ak")
    logger.info(f"Loading raw file took {np.round(time.time() - time_, 0)} seconds.")

    # Optionally load PSP data (not for PHY data)
    if config["DIRS"]["TYPE"] != "phy":
        psp_data, _ = store.read(f"ch{channel_id}/dsp/", psp_files, idx=ids)
        psp_data = psp_data.view_as("ak")
        psp_fields = psp_data.fields

    # Apply event mask to PHT data
    pht_data = pht_data[total_mask]

    # Add RAW information
    pht_data["baseline"] = raw_data["baseline"]
    pht_data["waveform"] = raw_data["waveform_windowed"]["values"]
    pht_data["channel_id"] = channel_id

    # Add PSP information (if present)
    if config["DIRS"]["TYPE"] != "phy":
        for field in psp_fields:
            if field != "timestamp":
                pht_data[field] = psp_data[field]

    # Attach peak labels (for peak modes)
    if config["energy_mode"] == "peaks":
        peak_list = []
        for peak_name, peak in config["peaks"].items():
            peak_mask = get_peak_region(
                pht_data, relevant_par_data, peak, config=config
            )
            subset = pht_data[peak_mask]
            subset["peak"] = peak_name
            peak_list.append(subset)
        merged_data = ak.concatenate(peak_list)

    elif config["energy_mode"] == "manual_peaks":
        peak_list = []
        for peak_name, peak in config["manual_peaks"].items():
            peak_mask = (
                (pht_data.trapEmax_ctc_cal > float(peak) - config["manual_peak_region"])
                & (pht_data.trapEmax_ctc_cal < float(peak) + config["manual_peak_region"])
            )
            subset = pht_data[peak_mask]
            subset["peak"] = peak_name
            peak_list.append(subset)
        merged_data = ak.concatenate(peak_list)

    else:
        merged_data = pht_data

    return merged_data



def quality_cut_mask(data: ak.Array) -> ak.Array:
    """ Python implementation of the LEGEND-200 quality cuts.
    
    Parameters
    ----------
    data: ak.Array
        Array on which the QC should be applied

    Returns
    -------
    mask: ak.Array
        A boolean mask of the input array
    """

    mask_qc = (
            ## Initial cuts
            ~data.is_discharge & # Get rid of waveforms with discharge
            ~data.is_saturated &  # Get rid of saturated events

            # Baseline cuts
            (data.is_valid_bl_slope == 1) & # get rid of rising & falling baseline slopes
            (data.is_valid_bl_slope_rms == 1) & # removes events with discontinuities in baseline

            # Tail cuts
            data.is_valid_tail_rms & # removes events with discontinuities in tail
            ~data.is_noise_burst & # Removes noise bursts

            # Empty energy cuts -> Not sure the trapTmax/Tmin are useful??
            ~data.is_valid_cuspEmax &
            data.is_valid_cuspEmin &
            data.is_valid_trap_tpmax &
            ~data.is_valid_trap_tpmin &
            # TODO: Do we need is_valid_trapEmax ?

            # Pulse start cuts
            data.is_valid_t0 & # eliminate early and late triggers
            data.is_valid_rt & # cut on rise times
            data.is_valid_dteff &# cut on effective drift time

            # Newly implemented cuts -> accurding to DSP slide
            ~data.is_downgoing_baseline &
            ~data.is_neg_energy & 
            ~data.is_negative & 
            ~data.is_upgoing_baseline &
            ~data.is_negative_crosstalk_old &
            data.is_valid_baseline_old & 
            data.is_valid_ediff & 
            data.is_valid_efrac & 
            data.is_valid_tmax &
            data.is_valid_neg_current &
    
            # Honestly not sure about those two cuts
            data.bl_pileup_cut &
            data.tail_pileup_cut
        )
    
    return mask_qc



def load_peak_mask(pht_data: ak.Array, relevant_par_data: dict, mask_pulser: ak.Array, 
                   logger: logging.Logger, config: dict) -> ak.Array:
    """
    Define the quality cuts and apply both quality cuts and pulser cuts. Afterwards, all peaks 
    defined in the config are filtered. Calls get_peak_region for the region around the peak energy. 
    Finally stores a total mask which is used to filter the pht data. This makes it possible to 
    only load relevant raw data, which greatly increases runtime.  

    Parameters
    ----------
    pht_data: ak.Array
        Array of physical hit events
    relevant_par_data: dict
        Additional parameters 
    mask_pulser: ak.Array
        Boolean mask of pulser events
    logger: logger.Logger
        Logger.
    config: dict
        Config file as parsed

    Returns
    -------
    total_peak_mask: ak.Array
        Boolean mask - whether or not the index is relevant 
    """

    mask_qc = quality_cut_mask(pht_data)
  
    # Loop over peaks
    # NOTE: WE CAN NOT APPLY ANY MASKS YET - we need the correct indicies to load raw data!!
    total_peak_mask = ak.Array(['False']*len(pht_data))
    for i, peak in enumerate(config['peaks'].values()):
        # Get peak mask
        peak_mask = get_peak_region(pht_data, relevant_par_data, peak, config=config)

        # Merge masks: OR operation because we start with FALSE values per default
        # and peaks are non-overlapping.
        if i==0:
            total_peak_mask = peak_mask
        else:
            total_peak_mask = (total_peak_mask | peak_mask)

        final_peak_mask = (total_peak_mask & ~mask_pulser & mask_qc)
        
    return final_peak_mask



def get_peak_region(data: ak.Array, relevant_par_data: dict, peak_config: float,
                    config: dict) -> ak.Array:
    """
    Select an energy region of the data that will be selected. In keV.
    
    Parameters
    ----------
    data: ak.Array
        Data on which the mask is based on
    relevant_par_data: dict
        Par data that contains the fwhm values
    peak_config: float
        Energy peak as defined in the config file
    config: dict
        Config file as parsed

    Returns
    -------
    peak_mask: ak.Array
        Mask of the selected peak region
    """
    # Create a dictionary which has the real-data peaks as values and removed-decimals as keys:
    peak_dict_data = {}
    for peak_in_data in relevant_par_data['results']['ecal']['trapEmax_ctc_runcal']['pk_fits'].keys():
        peak_dict_data[str(peak_in_data).split('.')[0]] = peak_in_data
    
    # Grab the key, i.e. the matching peak, e.g. '2614' -> '2614.533'
    matching_peak = peak_dict_data[str(peak_config)]
    fwhm = (relevant_par_data['results']['ecal']['trapEmax_ctc_runcal']['pk_fits'][str(matching_peak)]['fwhm_in_kev']) 
    cut_value = (config['cut_value_sigmas']*fwhm)/(2*np.sqrt(2*np.log(2))) # SCALE to std deviations

    # Calculate mask
    peak_mask =  ((data.trapEmax_ctc_cal > (float(matching_peak) - cut_value)) &  
                (data.trapEmax_ctc_cal < (float(matching_peak) + cut_value)))
    
    return peak_mask



def filter_metadata(legend_metadata: dict) -> dict:
    """
    Load Legend Metadata and select only valid_detector_types. 
    Also creates a list with all expected rawids of all detector types.

    Parameters
    ----------
    legend_metadata: dict
        Instance LegendMetadata which contains LEGEND metadata

    Returns
    -------
    detectors: dict
        Metadata containing name and id of the expected detectors
    """

    detectors = {}
    icpc_detectors = []
    ppc_detectors = []
    bege_detectors = []

    # Loop over legend metadata dictionary
    for key, item in legend_metadata.channelmap().items():

        # Skip non-geds detector systems
        if not item['system'] == 'geds':
            continue
        detector_dict = {}

        # Add important fields
        detector_dict['rawid'] = item['daq']['rawid']
        detector_dict['type'] = item['type']
        detector_dict['manufacturer'] = item['production']['manufacturer']
        detector_dict['processable'] = item['analysis']['processable']
        detector_dict['lq'] = item['analysis']['psd']['status']['lq']
        detectors[key] = detector_dict

        # Create lists with rawids
        if item['type'] == 'icpc':
            icpc_detectors.append(item['daq']['rawid'])
        elif item['type'] == 'ppc':
            ppc_detectors.append(item['daq']['rawid'])
        elif item['type'] == 'bege':
            bege_detectors.append(item['daq']['rawid'])
        elif item['type'] == 'coax':
            pass # skip for now
        else:
            pass # dont need this

        # Fill lists into dictionary
        detectors['rawid_lists'] = {
            'ICPC': icpc_detectors,
            'PPC': ppc_detectors,
            'BEGe': bege_detectors,
            }

    return detectors



def filter_valid_detectors(legend_metadata: dict, timestamp: dt.datetime, detector_type: str) -> dict:
    """ 
    Filter LEGEND metadata for a specific timepoint to get the active
    and usable detectors at that point in time.

    Parameters
    ----------
    legend_metadata: dict
        Instance LegendMetadata which contains LEGEND metadata
    timestamp: dt.datetime
        Datetime point at which the metadata should be evaluated
    detector_type: str
        Detector type as parsed or defined in the config

    Returns
    -------
    detectors: dict
         Dictionary of Dictionary that contains relevant at t = timestamp
    """

    chmap = legend_metadata.channelmap(timestamp)
    detectors = {}

    # Loop over chmap dictionary
    for name, item in chmap.items():
        inner_dict = {}
        # Filter geds and specific detector type
        if item.system == 'geds':
            if item.type == detector_type.lower():
                
                # Apply 'useable' label
                if ((item.analysis.processable == True) & (item.analysis.usability == 'on')):
                    inner_dict['useable'] = True
                else:
                    inner_dict['useable'] = False

                # Add some more metadata that will be useful
                inner_dict['rawid'] = item.daq.rawid
                inner_dict['manufacturer'] = item.production.manufacturer
                
                # Fill inner dictionary into outer dictionary
                detectors[name] = inner_dict

    return detectors
