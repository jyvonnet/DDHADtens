# 2D version
This document provides instructions to run the Matlab and Python implementation of DDHADtens 2D, which computes internal damage variables from collections of 3x3 matrices.
# Matlab version
Run DDHADtens2D.m. All files should be in the same directory, including RVEA_2D.mat... RVEE_2D.mat
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
The program requires MATLAB .mat files containing the elastic tensor data. The expected format is a variable named 'CC_save2' of size (3,3,N), where N is the number of time/load steps.

Available example files (to be placed in the same directory as the script):
- RVEA_2D.mat
- RVEB_2D.mat
- RVEC_2D.mat
- RVED_2D.mat
- RVEE_2D.mat
Usage
1. Open a terminal in the directory containing DDHAD2D_v6.py and the .mat files.
2. Run the script with:
       python DDHAD2D_v6.py
3. Choose the database to load (1–5) when prompted.
4. The program will:
   - Compute eta^0 and kappa^0 functions from the initial tensor
   - Construct the operator C_tilde
   - Compute internal variables alpha for each tensor in the database
   - Plot:
       * Evolution of alpha_i
       * Comparison of original vs reconstructed tensor components
       * 2D damage surface d(theta)
       * Reduced internal variables beta_i (via PCA)
Output
The program displays several figures during execution. No files are written to disk by default. If desired, you can save results or figures by adding save commands in the script.

