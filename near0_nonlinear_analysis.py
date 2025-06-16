import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import fitutils as ft
import tools
from schottky_analysis import schottky, T_max, E_exp, n_exp

# Data

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature


squared_temperature = temperature**2  # K**2
C_div_T = sample_HC/temperature  # mJ/K**2.mol
err_C_divT = err_sample_HC/temperature
err_squared_temperature = 2*temperature*err_temperature

# Constantes et tableaux

k = 1.380649e-23
r = 8.31446261815324  # J/mol.K

# Fit function


def fit_func(x, beta, gamma, n, E):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x
    y = E/(k*np.sqrt(x))
    schottky = (y**2) * np.exp(y)/((np.exp(y) + 1)**2)
    return phonon + gamma + n*r*1e3*schottky/np.sqrt(x)

# Intervals


def interval(a, b):
    """Give the wanted values in squared temperature and C/T, a, b : lower and higher bounds of temperature (K)"""
    temperature_bounded, C_div_T_bounded = tools.tab_interval(
        temperature, C_div_T, a, b)
    temperature_bounded, squared_temperature_bounded = tools.tab_interval(
        temperature, squared_temperature, a, b)
    temperature_bounded, err_squared_temperature_bounded = tools.tab_interval(
        temperature, err_squared_temperature, a, b)
    temperature_bounded, err_C_div_T_bounded = tools.tab_interval(
        temperature, err_C_divT, a, b)
    return temperature_bounded, squared_temperature_bounded, C_div_T_bounded, err_squared_temperature_bounded, err_C_div_T_bounded

# Non linear fit


def nonlinear_fit(a, b):
    temperature_bounded, squared_temperature_bounded, C_div_T_bounded, err_squared_temperature_bounded, err_C_div_T_bounded = interval(
        a, b)
    fit = opt.curve_fit(fit_func, squared_temperature_bounded, C_div_T_bounded, bounds=([0.1, 0, 1e-3, 1e-23], [1, 40, 1e-2, 1e-22]),
                        sigma=err_C_div_T_bounded, absolute_sigma=True)
    return fit[0]


def plot_fit(a, b):
    beta, gamma, n, E = nonlinear_fit(a, b)
    temperature_bounded, squared_temperature_bounded, C_div_T_bounded, err_squared_temperature_bounded, err_C_div_T_bounded = interval(
        a, b)
    plt.figure()
    plt.plot(squared_temperature_bounded, fit_func(
        squared_temperature_bounded, beta, gamma, n, E), "-y", label="fit")
    plt.plot(squared_temperature_bounded, C_div_T_bounded,
             ".g", label="C/T xperimental")
    plt.grid(True)
    plt.xlabel('T² (K²)')
    plt.ylabel('C/T (mJ/K².mol)')
    plt.legend()
    plt.show()

# Debye temperature


def debye_temperature(a, b, N=78e23):
    """
    Calculate the Debye temperature from the fit parameters
    Returns the Debye temperature in K and gamma in J/K².mol"""
    beta = nonlinear_fit(a, b)[0]*1e-3  # conversion en J/K⁴.mol
    pi4 = np.pi**4
    theta_D = (N*k*pi4*12)/(5*beta)  # en K³
    return np.cbrt(theta_D)


def main():
    plot_fit(0, 10)
    print(nonlinear_fit(0, 10))
    print(debye_temperature(0, 10))


if __name__ == "__main__":
    main()
