AYBCO Calorimetric Analysis – 04/06

Files:
- results_YBCO.ipynb: Main Jupyter notebook containing the analysis steps, visualizations, and summarized results.
- Python scripts: Python scripts used to perform the calculations in a modular way.
- H20250604(P250415 HPHT 475C 6GPa 3.5h YBCO7).dat: Experimental data file used for the analyses.

This folder gathers the analyses carried out this week on the heat capacity data of an optimally doped YBCO sample.

Objectives:
- Visualize the experimental data to identify thermal anomalies.
- Detect a Schottky-type contribution.
- Estimate the Debye temperature through linear and non-linear fits.
- Extract the associated physical parameters: γ, β, T_Debye, etc.

Required dependencies :
- Python > 3.8
- Libraries : numpy, matplotlib.pyplot, scipy.optimize, fitutils