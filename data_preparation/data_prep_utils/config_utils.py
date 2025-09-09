"""
Utility functions for handling configuration files.
"""
import importlib.util
import logging
import os



def load_config(path: str) -> dict:
    """
    Load a Python config file that defines GLOBAL_PARAM.

    Parameters
    ----------
    path : str
        Path to the config file (e.g. configs/config_co56_example.py).

    Returns
    -------
    dict
        The GLOBAL_PARAM dictionary from the config file.

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    AttributeError
        If the config file does not define GLOBAL_PARAM.
    """
    spec = importlib.util.spec_from_file_location("config", path)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)

    if not hasattr(config, "GLOBAL_PARAM"):
        raise AttributeError("Config file must define GLOBAL_PARAM")

    return config.GLOBAL_PARAM


 
def check_config(config_file: dict, logger: logging.Logger) -> bool:
    """
    Check the config file for the most important errors.

    Parameters
    ----------
    config_file: dict
        The configuration file to check
    logger: logging.Logger
        Logging info, warnings and errors.

    Returns
    -------
    config_validated
        Boolean value whether or not the config file is valid.
    """
    # Define variables
    config_validated = True


    # Info: Print general info
    if config_file['debug_mode'] == True:
        logger.info('You chose the debug mode. Not all cuts will be applied.')

    # Validation 1: Will events be stored? -> critical 
    if config_file['SAVE']['as_pickle'] == False:
        logger.error('No events will be saved. Please modify the config file.')
        config_validated = False

    # Validation 2: Will files be replaced?
    if config_file['FILTER']['replace_files'] == False:
        logger.warning('No files will be replaced. Events might not be stored!')

    # Validation 3: Is the input path valid? -> critical
    if not os.path.exists(config_file['DIRS']['LEGEND_RAW_DIR']):
        logger.error('The path LEGEND_RAW_DIR is not valid. Files can not be loaded.')
        config_validated = False

    # Validation 4: Is the given detector_type in the filtering part valid? -> critical
    if isinstance(config_file['FILTER']['detector_types_to_filter'], list):
        for ele in config_file['FILTER']['detector_types_to_filter']:
            if ele not in config_file['valid_detector_types']:
                config_validated = False

    # Validation 5: Make sure output directory exists (or create it)
    os.makedirs(config_file['DIRS']['output_dir'], exist_ok=True)
    
    return config_validated

