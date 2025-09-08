import numpy as np

config = {

    'save_figures': True,
    'calculate_efficiencies': False, # re-calculate efficiencies
    '2vbb':
    {
        'path_in': 'path/in/',
        'l_B': 0.087, 
        'l_B_err': 0.017, 
        'l_Bp': 0.5,
        'l_Bp_err': 0.01/np.sqrt(12),
     }, 

    'DEP_Th228':
    {
        'path_in': 'path/in/',
    },

    'DEP_Co56':
    {
        'path_in': 'path/in/',
    },

    'models': ['AoE', 'T1', 'T2'],
    'global_path_out': '/path/to/output_dir/', 
    'plot_dir': '/path/to/plots/',
    'path_to_metadata': '/path/to/metadata/',
}