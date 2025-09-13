"""
Plot all figures in scripts/ for figure plotting
"""

from thesis_figures.scripts import (
    activation_functions,
    charge_trapping,
    attenuation,
    charge_collection,
    double_beta_spectra,
    loss_functions,
    positional_encoding,
    waveform_example,
) 





def main():
    print("Generating all figures...")

    activation_functions.main()
    charge_trapping.main()
    attenuation.main()
    charge_collection.main()
    double_beta_spectra.main()
    loss_functions.main()
    positional_encoding.main()
    waveform_example.main()

    print("All figures generated successfully.")


if __name__ == "__main__":
    main()
