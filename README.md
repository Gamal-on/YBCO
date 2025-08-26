# YBCO Heat Capacity Analysis

Processing and analysis of low-temperature specific-heat measurements on YBa₂Cu₃O₇₊ₓ samples (optimally doped and heavily overdoped). The repository contains reproducible fitting routines, notebooks for data visualization, and the extracted fit parameters. 

![Specific heat of YBCO : ](images/data_plot_paper.png)

## Installation and Execution

Clone the repository and install dependencies using:

```bash
pip install -r requirements.txt
```

## Objective:

- Visualize experimental specific-heat data for three samples (YBCO_HPHT, YBCO_P240917, YBCO_ISIS).

- Compare theoretical models (Debye — with/without integral approximation, Debye+T⁵, Einstein, Schottky term) and fitting methods (least-squares / gradient descent, Monte-Carlo).

- Quantify the effect of heavy overdoping on the electronic (γ) and phononic (β, θD) contributions.

## Files:

- Python scripts : main files to perform the fit

- plot_data.ipynb: Notebook to visualize the experimental data of the three measured samples: optimally doped YBCO synthesized at high temperature and high pressure (YBCO_HPHT), and overdoped YBCO (YBCO_P240917 and YBCO_ISIS)

- results_YBCO_HPHT.ipynb: Notebook containing the linear fits for undoped YBCO (reference)

- results_YBCO_P24.ipynb: Notebook containing the linear fits for overdoped YBCO (x ≈ 0.3)

- results_YBCO_ISIS.ipynb: Notebook containing the linear fits for overdoped YBCO, cold‐pressed (x ≈ 0.3)

- results_values.md : Markdown file containing all the parameters determined by fitting the data

- image.ipynb : To plot the data on the same figure as previous reported ones

## Notes on method

- Analysis focuses on low-T range (< 20 K) to separate electronic (γT) and phononic (βT³ …) contributions; Schottky anomalies can dominate below ≈4 K and bias γ.

- Different fitting strategies (linear fits in selected T windows vs full nonlinear fits including Schottky and full Debye integral) can yield systematically different γ and β — consult the notebooks and results_values.md for comparisons.