# Visualization

This directory contains the code used to generate figures presented in my master's thesis.  
**Note**: The scripts here require prepared data to run. You can re-create the steps to obtain 
that data with `data_preparation` and a Transformer model. Each script corresponds to one or 
more figures, and can be run independently to reproduce the plots.  

---

## Scripts

* Calibration spectrum of Th-228
* Energy spectrum - A/E vs. Transformer on 2vbb events
* Model performance - Transformer classification scores and plots
* UMAP clustering - Unsupervised learnig for labelling


## Usage

### Generate a Single Figure
Run from root (inside legend200-ml-analysis/) any script directly, e.g.:
```bash
python -m visualization.scripts.calibration_spectrum
```

By default, outputs are saved in the same working directory as `thesis_figures/`.

### Generate all Figures
Run (again from root):
```bash
python -m visualization.make_all_figures
```


## Notes

* Always run all scripts from project root.

* The data used for the plots here is not public, so you need to 
    either ask for access or re-create your own.
