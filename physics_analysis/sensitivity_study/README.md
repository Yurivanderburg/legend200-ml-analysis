# Bayesian Sensitivity study with ZeroNuFit.jl for $0 \nu \beta \beta$


## Typical workflow

1. Clone ZeroNuFit.jl:

    ```bash
    git clone https://github.com/legend-exp/ZeroNuFit.jl.git
    ```

2. Edit configuration files:
    - fake_config.json: Fit configuration
    - fake_partitions.json: Define analysis window & nuisance parameters

3. Run ```sensitivity.jl```, e.g. via slurm scripts:

    ```bash
    sbatch submit_job_AoE.sh
    ```

4. Calculate summary and create plots run_analysis.py:

    ```bash
    python -m physics_analysis.sensitivity_study.run_analysis --path path/to/results --run_calc yes
    ```