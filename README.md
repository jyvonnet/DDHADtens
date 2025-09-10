# 3D version
This document provides instructions to run the Matlab and Python implementation of DDHADtens 3D, which computes internal damage variables from collections of 6x6 elastic tensors.
# Matlab version
Run run DDHADtens3D.m. All files should be in the same directory, including RVEA_3D.mat... RVED_3D.mat
# Python version
Requirements
- Python 3.9 or later
- Required libraries:
  * numpy
  * matplotlib
  * scipy

To install the dependencies, run:
    pip install numpy matplotlib scipy
Input Data
The program requires MATLAB .mat files containing the elastic tensor data. The expected format is a variable named 'CC_save2' of size (6,6,N), where N is the number of time/load steps.

Available example files (to be placed in the same directory as the script):
- RVEA_3D.mat
- RVEB_3D.mat
- RVEC_3D.mat
- RVED_3D.mat
Usage
1. Open a terminal in the directory containing DDHADV3_v1.py and the .mat files.
2. Run the script with:
       python DDHADV3_v1.py
3. Choose the database to load (1–4) when prompted.
4. The program will:
   - Compute eta^0 and kappa^0 functions from the initial tensor
   - Construct the operator C_tilde
   - Compute internal variables alpha for each tensor in the database
   - Plot:
       * Evolution of alpha_i
       * Comparison of original vs reconstructed tensor components
       * 3D damage surface d(theta,phi)
       * Reduced internal variables beta_i (via PCA)

Output
The program displays several figures during execution. No files are written to disk by default. If desired, you can save results or figures by adding save commands in the script.
