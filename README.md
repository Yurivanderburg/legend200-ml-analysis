# LEGEND-200 ML Analysis

This repository contains the code and utilities used in the 
**LEGEND-200 machine learning–based pulse-shape discrimination (PSD) efficiency and sensitivity studies**
, as well as data preparation and plotting scripts for the associated thesis work.

---

## Repository Structure
```
├── data_preparation/              # Pre-process data for ML training
│   ├── configs/                   # Configuration files <- **Edit here**
│   ├── data_prep_utils/           # Utilities for processing data
│   ├── scripts/                   # Scripts <- **Run here**
│   └── slurm_jobs/                # SLURM submission .sh files for the python code
├── logs/                          # Logging for all python scripts 
├── physics_analysis/              # Core physics analysis
│   ├── psd_efficiency/            # Pulse Shape Discrimination at $Q_{\beta \beta}$
│   ├── pulse_simulation/          # Pulse Shape Simulation codes 
│   └── sensitivity_study/         # Bayesian sensitivity study for $0 \nu \beta \beta$ decay half-life
├── plotting/                      # *WIP*
├── thesis/                        # LaTeX files for the Master's thesis
├── defense/                       # LaTeX files for the Master's thesis defense talk
└── utils/                         # Various utilities
```

---

## Getting Started

### Requirements
- Python 3.11+ (tested)
- Recommended: create a virtual environment

Install dependencies:
```bash
pip install -r requirements.txt
```


### Data Preparation

Contains the pipeline for preparing LEGEND-200 waveform data for
training an encoder-only Transformer for event classification. It handles 
pre-loading, filtering, dataset construction, data cleaning, format conversion, 
and logging. The goal is to produce standardized, ML-ready datasets that can 
be used by the transformer model and physics analysis modules. 


### Physics Analysis

Contains the code physics analyses performed in the scope of this work. 
They include work on Pulse Shape Simulation (with LegendGeSim.jl), some
MaGe (GEANT4 based) simulations and custom made scripts to combine both. 
Also contains the code for estimating the PSD efficiency at $Q_{\beta \beta}$
and analysis for the Bayesian sensitivity study to estimate the 
$0 \nu \beta \beta$ decay half-life of the experiment, based on the PSD 
efficiencies. 
