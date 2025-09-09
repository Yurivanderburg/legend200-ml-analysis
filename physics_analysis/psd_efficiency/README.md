# Pulse Shape Discrimination efficiency at DEP, $2 \nu \beta \beta$ events and $Q_{\beta \beta}$


## Typical workflow

1. Edit configuration files:
    - psd_config.json: Specify settings and input/output paths

2. Run main analysis:
    - Calculates PSD efficiency at three energies ([1000, 1300], 1592.5 and 2231.0 keV)

    ```bash
    python -m physics_analysis.sensitivity_study.run_psd_efficiency
    ```