# Physics Analysis 

This part contains the main physics analyses for the LEGEND-200 ML 
study, including **PSD efficiency evaluation**, **pulse simulation**, 
and **sensitivity studies**. Each submodule is designed to be run 
standalone but shares utilities where possible.

---

## PSD Efficiency

- **Purpose:** Compute PSD survival fractions and uncertainties for physics 
    data and simulated peaks, then interpolate efficiencies to the double-beta decay Q-value.
- **Key components:**
  - `run_psd_efficiency.py` — main driver script.
  - `config/psd_config.py` — example configuration of datasets and analysis parameters.
  - `psd_utils/efficiency/`:
    - `double_escape_peaks.py` — DEP efficiency fits.
    - `two_neutrino_bb.py` — $2 \nu \beta \beta$ efficiency calculation.
    - `interpolation_to_qbb.py` — extrapolation of efficiencies to Qββ.
  - `psd_utils/visualization.py` — plotting utilities.

**Run example:**
```bash
python -m physics_analysis.psd_efficiency.run_psd_efficiency
```

---

## Pulse simulation

- **Purpose:** Generate and simulate detector geometry, physics events and waveforms.
- **Key components:**
  - `config/complete_sis.json` - Source Insertion System geometry for GEANT4 simulation.
  - `config/legend_settings.mac` — Macro file for GEANT4 simulation
  - `scripts/create_pet_file.py` - Position, Energy, Time file for GEANT4 -> LegendGeSim.jl
  - `scripts/simulate_waveforms.jl` - Julia Code for simulating waveforms (based on LegendGeSim.jl)
  - `slurm_jobs/submit_remage_job.sh` - Slurm script to submit the GEANT4 event simulation to a Cluster

### Notes:
* Requires [LegendGeSim.jl](https://github.com/legend-exp/LegendGeSim.jl) and [MaGe](https://arxiv.org/abs/0802.0860)


---

## Sensitivity study

- **Purpose:** Quantify the sensitivity improvement of ML-based PSD methods compared to baseline cuts (AoE). 
    Also creates various plots for visualization.
- **Key components:**
  - `config/` - Config files used for ZeroNuFit.jl fits (see documentation there)
  - `sensitivity_utils/collect_json.py` — Calculate statistics and summary of Monte Carlo fits.
  - `sensitivity_utils/visualize_sensitivity.py` — Plotting utilities
  - `slurm_jobs/` - Submit ZeroNuFit jobs to cluster
  - `run_analysis.py` - Main entry to calculate summary and create plots.


### Notes:
* Requires [ZeroNuFit.jl](https://github.com/legend-exp/ZeroNuFit.jl) 
* First, perform the fits with ZeroNuFit, and then calculate the summary with:

```bash
python -m physics_analysis.sensitivity_study.run_analysis --path path/to/results --run_calc yes
```