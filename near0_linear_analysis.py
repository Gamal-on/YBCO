import numpy as np
import matplotlib.pyplot as plt
import tools
import fitutils as ft
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


# Substracted values

def interval(a, b,  E=E_exp, n=n_exp):
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


def plot_substracted(a, b,  E=E_exp, n=n_exp):
    """    Plot the substracted values of C/T - C_schottky vs T² in a given interval """
    temperature_bounded, squared_temperature_bounded, C_div_T_bounded = interval(
        a, b, E, n)[0:3]
    C_div_T_substracted = C_div_T_bounded - \
        schottky(temperature_bounded, E, n)/temperature_bounded
    return squared_temperature_bounded, C_div_T_substracted

# Linear fit : (C/T - C_schottky) (T²)


def linear_fit(a, b, E=E_exp, n=n_exp):
    """
    Perform a linear fit on the substracted values of C/T - C_schottky vs T²
    Returns the fit parameters (beta, gamma, n) in (mJ/K⁴.mol, mJ/K².mol)"""
    err_squared_temperature_bounded, err_C_div_T_bounded = interval(
        a, b, E, n)[3:5]
    fit = ft.linfitxy(
        plot_substracted(a, b, E, n)[0], plot_substracted(a, b, E, n)[
            1], err_squared_temperature_bounded, err_C_div_T_bounded, plot=True,
        markercolor="r", linecolor="orange")
    plt.show()
    return fit


# Debye temperature and gamma

def debye_temperature(a, b, E=E_exp, n=n_exp, N=78e23):
    """
    Calculate the Debye temperature and gamma from the linear fit parameters
    Returns the Debye temperature in K, gamma in J/K².mol and their respectiv errors"""
    beta, gamma, u_beta, u_gamma = linear_fit(
        a, b, E, n)*1e-3  # conversion en J
    pi4 = np.pi**4
    theta_D = (N*k*pi4*12)/(5*beta)  # en K³
    u_theta_D = np.cbrt(theta_D) * u_beta/(3*beta)
    temp_debye = np.cbrt(theta_D)
    return temp_debye, gamma, u_theta_D, u_gamma


def main():
    print(debye_temperature(0, 10, E=E_exp, n=n_exp))


if __name__ == "__main__":
    main()
