import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import fitutils as ft
from schottky_analysis import schottky, T_max, E_exp

# Data

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature

squared_temperature = temperature**2  # K**2
C_div_T = sample_HC/temperature  # mJ/K**2.mol
err_C_divT = err_sample_HC/temperature

# Constantes et tableaux

k = 1.380649e-23
r = 8.31446261815324  # J/mol.K

# Fit function


def fit_func(x, beta, gamma, n):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x
    y = E_exp/(k*np.sqrt(x))
    schottky = (y**2) * np.exp(y)/((np.exp(y) + 1)**2)
    return phonon + gamma + n*r*1e3*schottky/np.sqrt(x)

# Non linear fit


def nonlinear_fit(N):
    fit = opt.curve_fit(fit_func, squared_temperature[0:N], C_div_T[0:N],
                        sigma=err_C_divT[0:N], absolute_sigma=True)

    return fit[0]


def main():
    print(nonlinear_fit(65, 1e-2))


if __name__ == "__main__":
    main()
