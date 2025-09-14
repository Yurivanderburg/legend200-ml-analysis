"""
Plot all figures in scripts/ for figure plotting
"""

from visualization.scripts import (
    calibration_spectrum,
    energy_spectrum_2vbb,
    model_performance,
    umap_clustering
) 



def main():
    print("Generating all figures...")

    calibration_spectrum.main()
    energy_spectrum_2vbb.main()
    model_performance.main()
    umap_clustering.main()

    print("All figures generated successfully.")


if __name__ == "__main__":
    main()
