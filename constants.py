import numpy as np
import matplotlib.pyplot as plt
import tools
import fitutils as ft

# Data

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature

squared_temperature = temperature**2  # K**2
C_div_T = sample_HC/temperature  # mJ/K**2.mol
err_C_divT = err_sample_HC/temperature
err_squared_temperature = 2*temperature*err_temperature

# Constants

k = 1.380649e-23
r = 8.31446261815324
N = 78.286e23

# Fitted parameters

# Determined by curve fit directly
n_curve_fit = 9.997e-3
E_curve_fit = 1.0966e-22  # J

# Deduced form Schottky analysis
T_exp = 2.9461005
E_exp_8 = 1.0681622692127778e-22  # with NR and 8.67
E_exp_10 = 1.1549944696942883e-22  # with mean and 10.148
n_exp = 0.0016374152345313721

# Determined by quadratic fit
beta_quadratic = 0.46318268614826286
gamma_quadratic = 0.06319440532206101
n_quadratic = 0.012979953809295845
E_quadratic = 1.1584047800714808e-22
alpha_quadratic = 0.00046894875156273756

# Debye temperature


def debye_temperature(beta):
    """
    Calculate the Debye temperature and gamma from the linear fit parameters
    Returns the Debye temperature in K, gamma in J/K².mol and their respectiv errors"""
    pi4 = np.pi**4
    factor = N*12*k*pi4/5
    temp_debye = np.cbrt(factor/(beta*1e-3))
    return temp_debye


def main():
    print(debye_temperature(0.46))


if __name__ == "__main__":
    main()
