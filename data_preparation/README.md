# Data Preparation

This module contains the pipeline for preparing LEGEND-200 waveform data for
training an encoder-only Transformer for event classification. It handles pre-loading, 
filtering, dataset construction, data cleaning, format conversion, and logging.

The goal is to produce standardized, ML-ready datasets that can be used by the
transformer model and physics analysis modules. 

---



## Typical workflow

1. **Prepare dataset**

   ```
   python -m data_preparation.scripts.prepare_pht_dataset --config data_preparation/configs/example_config_pht.py
   ```

2. **Filter dataset**

   ```bash
   python -m data_preparation.scripts.filter_pht_dataset --config data_preparation/configs/config_pht.py 
   ```

3. **Convert to HDF5**

   ```
   python -m data_preparation.scripts.pickle_to_hdf5.py input.pkl output.hdf5
   ```

Logs are written to `logs/` automatically.

---

# Workflows

|         Script           | Config file    | Purpose                                  |
|--------------------------|----------------|------------------------------------------|
| `prepare_phy_dataset`    | config_phy     | Load physics data (skm + waveform)       | 
| `prepare_evt_dataset`    | config_evt     | Load cal. data (evt tier + waveform)     | 
| `prepare_pht_dataset`    | config_pht     | Load cal. data (hit tier + waveform)     | 
| `filter_pht_dataset`     | config_phy     | Filter & label hit tier data             |


** You need to always run scripts directly from ```legend200-ml-analysis/```, e.g.,: 
```bash
python -m data_preparation.scripts.prepare_phy_dataset --config data_preparation/configs/config_phy.py 
```
---

## Notes

* LEGEND-200 data is **not included** in this repository.
* Only HIT-tier level data is used for training and therefore needs the 
labeling from `filter_dataset`; Both physics data and high-level EVT-tier already
have calibrated and approved PSD labels (is_bb_like). 

