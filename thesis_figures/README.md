# Thesis Figures

This directory contains the code used to generate one of the figures presented in my master's thesis.  
Each script corresponds to one or more figures, and can be run independently to reproduce the plots.  

---

## Scripts

* Activation functions
* Charge trapping in HPGe detectors
* X-Ray attenuation in Germanium
* Charge collection in HPGe detectors
* 0vbb and 2vbb decay spectra
* Loss functions
* Positional encoding used in the Transformer
* Waveform example -> A/E



## Usage

### Generate a Single Figure
Run from root (inside legend200-ml-analysis/) any script directly, e.g.:
```bash
python -m thesis_figures.scripts.03_attenuation
```

By default, outputs are saved in the working directory or as defined in each script.

### Generate all Figures
Run (again from root):
```bash
python -m thesis_figures.plot_all
```

## Run Animations

Manim scripts are stored in `manim/`. Example:

```bash
manim -pql manim/MLP_diagram.py MLPDiagram
```

For documentation, check out [Manim](https://github.com/manimCommunity/manim)


## Notes

* Always run all scripts from project root.

* All data used for plotting (NIST attenuation tables, 
    Ge-76 spectra, waveform examples) is stored in utils/.

* Shared helper functions are in utils/figures_utils.py.

* Animations are separated from static figures for clarity.