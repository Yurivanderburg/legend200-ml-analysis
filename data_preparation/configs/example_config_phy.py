# Example config for data preparation of physics data

GLOBAL_PARAM = {
    # ------------------------------------------------------------
    # General
    # ------------------------------------------------------------
    'detector_type': 'ICPC',
    'debug_mode': False,
    'remove_pulser': False,
    'SSH_AUTH': '',


    # ------------------------------------------------------------
    # Directories
    # ------------------------------------------------------------
        'DIRS': {
            'LEGEND_RAW_DIR': 'path/to/raw',
            'LEGEND_SKM_DIR': 'path/to/skm',
            'PERIODS': ['p03', 'p04', 'p06'],
            'output_dir': 'path/to/output',
            'filtered_subdir': 'FILTERED',
            'lmeta': 'path/to/lmeta'

    },

    # ------------------------------------------------------------
    # Energy selection
    # ------------------------------------------------------------
    'energy_mode': 'continuum',  # options: 'continuum', 'peaks', 'manual_peaks'

    'peaks': {                   # used if energy_mode = 'peaks'
        'DEP_Tl': 1592,
        'FEP_Bi': 1620,
        'SEP_Tl': 2103,
        'FEP_Tl': 2614,
    },

    'manual_peaks': {
        '8 -> 2': 1771.42,
        '10 -> 2': 2025.35,
        '11 -> 2': 2034.92,
        '7 -> 1': 2598.57,
        '9 -> 1': 3202.19,
        '10 -> 1': 3253.64,
        '11 -> 1': 3273.19,
    },

    'manual_peak_region': 2.5,       # keV around manual peaks
    'continuous_region': [2114, 2118],  # keV

    # ------------------------------------------------------------
    # Cuts and labels
    # ------------------------------------------------------------
    'add_labels': False,
    'cut_value_sigmas': 3,
    'cut_value_sigma_DEP': 5,
    'valid_detector_types': ['PPC', 'ICPC', 'BEGe'],
    'exclude_period_run': [],

    # ------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------
    'FILTER': {
        'detector_types_to_filter': ['ICPC'],   # or 'all'
        'periods_to_filter': ['p03'],
        'merge_surface_labels': True,
        'replace_files': True,
        'energy_cut_on_ncontact': False,

        'outlier_cuts': [
            'rt1', 'rt2', 'tp_80', 'tp_50', 'tp_10',
            'AoE_Corrected', 'LQ_Classifier', 'scp_widths'
        ],
        'cut_value_outliers': 3,
        'cut_value_width': 2,

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
