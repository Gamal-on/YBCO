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


def plot_fit(N, beta, gamma, n):
    plt.figure()
    plt.plot(squared_temperature[0:N], fit_func(
        squared_temperature[0:N], beta, gamma, n), "g--", label="Fit")
    plt.plot(squared_temperature[0:N], C_div_T[0:N],
             ".b", label="C/T xperimental")
    plt.grid(True)
    plt.xlabel(r'T² (K²)')
    plt.ylabel(r'C/T (mJ/K².mol)')
    plt.legend()
    plt.show()

# Debye temperature


def debye_temperature(N):
    """
    Calculate the Debye temperature from the fit parameters
    Returns the Debye temperature in K and gamma in J/K².mol"""
    beta = nonlinear_fit(N)[0]*1e-3  # conversion en J/K⁴.mol
    pi4 = np.pi**4
    theta_D = (r*pi4*12)/(5*beta)  # en K³
    return np.cbrt(theta_D)


# Main function to run the analysis and plot the fit

def final():
    N = int(input("Enter the number of data points to fit (N): "))
    beta, gamma, n = nonlinear_fit(N)  # en mJ
    plot_fit(N, beta, gamma, n)
    print("Debye temperature:", debye_temperature(N),
          "K", "Gamma : ", gamma, "mJ/K².mol", "n : ", n),


def main():
    final()


if __name__ == "__main__":
    main()
