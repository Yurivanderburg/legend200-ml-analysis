# Pulse Shape Simulation


## Typical workflow

1. **Detector simulation**

* Simulate the detector geometry for a given a source-configuration file:

    ```legend-pygeom-l200 --visualize --fiber-modules segmented --config config/complete_sis.json path/to/geometry_file.gdml```


2. **Energy deposition simulation**

* Simulate physics events with ReMaGe for given a geometry file and a given macro file:

    ```remage --threads $N_THREADS -g path/to/geometry_file.gdml -o $OUTPUT_FILE macro_file.mac```


3. **Create PET file**

* Convert GEANT4 output into PET file:

    ```python -m physics_analysis.pulse_simulation.create_pet_file```


4. **Waveform simulation**

* For a given PET file (specified in the julia file), simulate raw waveforms:

    ```julia simulate_waveforms.jl```
