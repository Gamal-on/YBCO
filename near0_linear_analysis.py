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

# Constantes et tableaux

k = 1.380649e-23
r = 8.31446261815324  # J/mol.K


# Substracted values


def plot_substracted(N, E=E_exp, n=n_exp):
    """
    Plot the substracted values of C/T - C_schottky vs T²"""
    C_divT_substracted = C_div_T - schottky(temperature,
                                            E_exp, n)/temperature  # mJ/K**2.mol
    plt.figure()
    plt.plot(squared_temperature[0:N], C_divT_substracted[0:N], 'g.',
             label="C/T - C_schottky")
    plt.xlabel("T² (K²)")
    plt.ylabel("C/T - C_schottky (mJ/K².mol)")
    plt.title("C/T - C_schottky vs T²")
    plt.show()


# Linear fit : (C/T - C_schottky) (T²)

def linear_fit(N, E=E_exp, n=n_exp):
    """
    Perform a linear fit on the substracted values of C/T - C_schottky vs T²
    Returns the fit parameters (beta, gamma, n) in mJ/K⁴.mol"""
    C_divT_substracted = C_div_T - schottky(temperature,
                                            E_exp, n)/temperature  # mJ/K**2.mol
    fit = ft.linfitxy(
        squared_temperature[0:N], C_divT_substracted[0:N], err_temperature[0:N], err_C_divT[0:N])
    return fit


def plot_linear_fit(N, E=E_exp, n=n_exp):
    """
    Plot the linear fit of the substracted values of C/T - C_schottky vs T²"""
    C_divT_substracted = C_div_T - schottky(temperature,
                                            E, n)/temperature  # mJ/K**2.mol
    fit = ft.linfitxy(squared_temperature[0:N], C_divT_substracted[0:N], err_temperature[0:N], err_C_divT[0:N],
                      plot=True, linecolor="red", markercolor="orange")
    plt.grid(True)
    plt.show()
    return fit


# Debye temperature and gamma

def debye_temperature(N, E=E_exp, n=n_exp):
    """
    Calculate the Debye temperature and gamma from the linear fit parameters
    Returns the Debye temperature in K and gamma in J/K².mol"""
    beta, gamma = linear_fit(N, n)[0:2]*1e-3  # conversion en J/K⁴.mol
    pi4 = np.pi**4
    theta_D = (r*pi4*12)/(5*beta)  # en K³
    return np.cbrt(theta_D), gamma


# Main function to run the analysis and plot results

def final():
    """
    Main function to run the analysis and plot results."""
    N = int(input("Enter the number of data points to consider (N): "))
    plot_linear_fit(N, E=E_exp, n=n_exp)
    debye_temp = debye_temperature(N, E=E_exp, n=n_exp)
    print("Debye temperature:", debye_temp[0],
          "K", "Gamma:", debye_temp[1], "J/K².mol")
    return debye_temp[0], debye_temp[1]


def main():
    final()


if __name__ == "__main__":
    main()
