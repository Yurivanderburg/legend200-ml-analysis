# Example config for data preparation of high-level physics data (CAL tier)

GLOBAL_PARAM = {
    # ------------------------------------------------------------
    # General
    # ------------------------------------------------------------
    'detector_type': 'ICPC',      # options: 'ICPC', 'BEGe', 'Coax'
    'debug_mode': False,          # skip cuts for inspection
    'remove_pulser': False,       # remove pulser events if True
    'SSH_AUTH': '',

    # ------------------------------------------------------------
    # Directories
    # ------------------------------------------------------------
    'DIRS': {
        'LEGEND_RAW_DIR': 'path/to/raw/waveforms',
        'LEGEND_SKM_DIR': 'path/to/calibration/data',
        'PERIODS': ['p09'],   # periods to include
        'output_dir': 'output/path',
        'filtered_subdir': 'FILTERED',
        'lmeta': 'path/to/Legendmetadata/'


    },

    # ------------------------------------------------------------
    # Energy selection
    # ------------------------------------------------------------
    'energy_mode': 'continuum',   # options: 'continuum', 'peaks', 'manual_peaks'

    'peaks': {                    # used if energy_mode = 'peaks'
        'DEP_Tl': 1592,
        'FEP_Bi': 1620,
        'SEP_Tl': 2103,
        'FEP_Tl': 2614,
    },

    'manual_peaks': {             # used if energy_mode = 'manual_peaks'
        # Source: https://journals.aps.org/pr/pdf/10.1103/PhysRev.184.1111
        '8 -> 2': 1771.42,
        '10 -> 2': 2025.35,
        '11 -> 2': 2034.92,
        '7 -> 1': 2598.57,
        '9 -> 1': 3202.19,
        '10 -> 1': 3253.64,
        '11 -> 1': 3273.19,
    },
    'manual_peak_region': 2.5,     # keV around manual peaks
    'continuous_region': [0, 5000],  # keV, used if energy_mode = 'continuum'

    # ------------------------------------------------------------
    # Cuts and labels
    # ------------------------------------------------------------
    'add_labels': False,
    'cut_value_sigmas': 3,        # std devs for selection window
    'cut_value_sigma_DEP': 5,     # exception for DEP
    'valid_detector_types': ['PPC', 'ICPC', 'BEGe'],
    'exclude_period_run': ['p05_r000'],  # bad periods to exclude

    # ------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------
    'FILTER': {
        'detector_types_to_filter': ['ICPC'],  # list or 'all'
        'periods_to_filter': ['p09'],          # list or 'all'
        'merge_surface_labels': True,
        'replace_files': True,
        'energy_cut_on_ncontact': False,

        # Properties for outlier cuts
        'outlier_cuts': [
            'rt1', 'rt2', 'tp_80', 'tp_50', 'tp_10',
            'AoE_Corrected', 'LQ_Classifier', 'scp_widths'
        ],
        'cut_value_outliers': 3,
        'cut_value_width': 2,

        # scipy.find_peaks parameters
        'find_peaks_MSE_p_n': {
            'prominence': 1e-3, 'distance': 2, 'width': 1, 'height': 0.15,
        },
        'find_peaks_SSE': {
            'prominence': 0.01, 'distance': 2, 'width': 1, 'height': 0.15,
        },
    },

    # ------------------------------------------------------------
    # Saving
    # ------------------------------------------------------------
    'SAVE': {
        'remove_fields': False,
        'peakfind_filter': False,
        'as_pickle': True,
    }
}
