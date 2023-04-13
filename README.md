# EBSD Data Processing with MTEX
This code demonstrates how to process EBSD (Electron Backscatter Diffraction) data using the MTEX toolbox in MATLAB. It includes various calculations and visualizations such as grain reconstruction, texture index, and orientation distribution function (ODF) plots.

# Prerequisites
- MATLAB (version R2017a or later)
- MTEX toolbox (version 5.2.2 or later)
# Getting Started
1. Clone or download this repository to your local machine.
2. Open MATLAB and navigate to the cloned/downloaded repository.
3. Update the adress and sample variables to match your file directory path and sample ID.
4. Run the EBSD_data_processing.m script.
# Contents
The EBSD_data_processing.m script includes the following steps:
1. Load EBSD data and calculate the Burgers orientation between Beta and Alpha phases.
2. Calculate grains from the indexed EBSD data and reconstruct the parent grains.
3. Calculate grain boundary votes and reconstruct parent grains.
4. Merge similar grains and inclusions.
5. Plot the Alpha phase with its orientations and overlay the parent grains' boundaries.
6. Calculate and plot the Beta phase ODF.
7. Calculate and plot the Beta phase texture index.
8. Calculate and plot the Alpha phase ODF and save as a TIFF image file.
9. Calculate and plot the texture index of the Alpha phase.
# Acknowledgments
This code is based on the MTEX tutorial examples and documentation. Special thanks to Dr. Franz Roters and the MTEX development team for creating this powerful toolbox.
