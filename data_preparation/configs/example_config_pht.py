# ============================================================
# Config file for data preparation
# ============================================================
# This file defines global parameters for preparing LEGEND-200
# waveform datasets. Copy this template and edit paths/parameters
# for your dataset (see configs/README.md for examples).
# ============================================================

GLOBAL_PARAM = {

    # ------------------------------------------------------------
    # General settings
    # ------------------------------------------------------------
    'detector_type': 'ICPC',       # Valid: 'ICPC', 'BEGe', 'Coax'
    'debug_mode': True,            # Skip cuts for inspection (True/False)
    'SSH_AUTH': '',

    # ------------------------------------------------------------
    # Directories
    # ------------------------------------------------------------
    'DIRS': {
        'LEGEND_DATA_DIR': 'path/to/datadir',
        'TYPE': 'cal',             # Dataset type: 'cal' or 'phy'

        # Subdirectories by data tier
        'TIER_SUBDIRS': {
            'PHT': 'path/to/pht', # higher-level attributes
            'TCM': 'path/to/tcm', # pulser events
            'PSP': 'path/to/psp', # DCR param -> Not yet used
        },
        'LEGEND_RAW_DIR': 'path/to/raw/waveforms',
        'PAR_SUBDIR': 'path/to/par',     # Calibration parameters
        'PERIODS': ['p03', 'p04'],       # Which data-taking periods to include
        'output_dir': 'path/to/output',
        'filtered_subdir': 'filtered',   # Subdir for filtered outputs
        'lmeta': 'path/to/Legendmetadata/'
    },

    # ------------------------------------------------------------
    # Energy selection
    # ------------------------------------------------------------
    'energy_mode': 'continuum',    # Options: 'continuum', 'peaks'

    'peaks': {                     # Only used if energy_mode='peaks'
        'DEP_Tl': 1592,
        'FEP_Bi': 1620,
        'SEP_Tl': 2103,
        'FEP_Tl': 2614
    },

    'continuous_region': [0, 5000],  # Used if energy_mode='continuum' (keV)

    # ------------------------------------------------------------
    # Labelling options
    # ------------------------------------------------------------
    'add_labels': False,
    'cut_value_sigmas': 3,         # Std. devs. for selection window
    'cut_value_sigma_DEP': 5,      # Exception for DEP

    'valid_detector_types': ['PPC', 'ICPC', 'BEGe'],
    'exclude_period_run': ['p05_r000'],  # Skip bad periods

    # ------------------------------------------------------------
    # Filtering parameters
    # ------------------------------------------------------------
    'FILTER': {
        'detector_types_to_filter': ['ICPC'],  # List or 'all'
        'periods_to_filter': ['p03'],          # List or 'all'
        'merge_surface_labels': False,
        'replace_files': True,
        'energy_cut_on_ncontact': False,

        'outlier_cuts': ['rt1', 'rt2', 'tp_80', 'tp_50',
                         'tp_10', 'AoE_Corrected', 'LQ_Classifier', 'scp_widths'],
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
    # Saving options
    # ------------------------------------------------------------
    'SAVE': {
        'remove_fields': False,
        'peakfind_filter': False,
        'as_pickle': True,
    }
}
