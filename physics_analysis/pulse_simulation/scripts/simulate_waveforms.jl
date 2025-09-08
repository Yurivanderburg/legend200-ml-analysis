using LegendGeSim
using HDF5
using LegendHDF5IO


# Parameters and paths
detector_metadata_filename = "path/to/legend-detectors/"
pet_filename = "path/to/pet_file"
output_filename = "path/to/output"



# Settings
environment_settings = Dict(
    "crystal_temperature_in_K" => 77,
    "medium" => "vacuum",
)

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "", # a non-empty string will cache the simulation results
)

daq_settings = Dict(

    # Settings for V02160A detector
    "preamp" => Dict(
        "type" => "generic",
        "t_decay_in_us" => 250, # from V04545A HADES data # changed from 43
        "t_rise_in_ns" => 100, # by eye
        "gain_ADC_eV" => 0.0154068, #estimated from previous run
        "offset_in_ADC" => 15000, # from V04545A HADES data mean() of baseline
        "noise_sigma_in_keV" => 1 # by eye
    ),
    "fadc" => Dict(
        "type" => "generic",
        "sampling_interval" => 16 # ns, like L200 data
    ),
    "trigger" => Dict(
        "type" => "trapezoidal",
        "window_lengths" => [250,250,250],
        "threshold" => 9 # keV
    ),
    "daq" => Dict(
        "type" => "generic",
        "nsamples" => 1400, # windowed waveforms in L200
        "baseline_length" => 375 # by eye from data
    )
)

noise_model = Dict(
    "type" => "sim"
)


## Simulation
raw_table = LegendGeSim.simulate_raw(detector_metadata_filename, pet_filename, environment_settings, simulation_settings, daq_settings, noise_model)  #; n_waveforms=1000)


## Store file - only half the length because both p+ and n+ signals are stored
half_length = convert(Int32, length(raw_table)/2)


h5open(output_filename, "w") do f
    LegendHDF5IO.writedata(f, "raw", raw_table[1:half_length])
end