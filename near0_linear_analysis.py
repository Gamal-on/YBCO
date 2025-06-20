import numpy as np
import matplotlib.pyplot as plt
import tools
import fitutils as ft
import schottky_analysis
import constants

# Substraction of Schottky contribution


def schottky_substraction(E, n, x, y):
    """Return the data C/T with Schottky contribution Cs(T) substracted :
    E : scalar, n : scalar, x = temperature, y = C/T"""
    substracted = y - schottky_analysis.schottky(x, E, n)/x
    return substracted

# Linear fit with Monte Carlo


def linear_fit(a, b, E, n, x_carre, y, err_x=constants.err_squared_temperature, err_y=constants.err_sample_HC):
    """
    Perform a linear fit on the substracted values of C/T - C_schottky vs T²
    Returns the fit parameters (beta, gamma, n) in (mJ/K⁴.mol, mJ/K².mol) 
    a, b : bounds, x = T² , y = C/T, err_x = u(T²), err_y = u(HC)"""
    y = schottky_substraction(E, n)
    x_interval, y_interval = tools.tab_interval(x, y, a, b)
    err_x_interval, err_y_interval = err_x[0:len(
        x_interval)], err_y[0:len(y_interval)]
    fit = ft.linfitxy(x_interval, y_interval, err_x_interval, err_y_interval)
    return fit


# Debye temperature and gamma

def debye_temperature(a, b, E, n, N=78e23):
    """
    Calculate the Debye temperature and gamma from the linear fit parameters
    Returns the Debye temperature in K, gamma in J/K².mol and their respectiv errors"""
    beta, gamma, u_beta, u_gamma = linear_fit(
        a, b, E, n)*1e-3  # conversion en J
    pi4 = np.pi**4
    theta_D = (N*constants.k*pi4*12)/(5*beta)  # en K³
    u_theta_D = np.cbrt(theta_D) * u_beta/(3*beta)
    temp_debye = np.cbrt(theta_D)
    return temp_debye, gamma, u_theta_D, u_gamma

# Plotting results


def plot_results(a, b, E, n, x=constants.squared_temperature,
                 err_x=constants.err_squared_temperature, err_y=constants.err_sample_HC):
    fit = linear_fit(a, b, E, n, x=constants.squared_temperature,
                     err_x=constants.err_squared_temperature, err_y=constants.err_sample_HC)
    a, b = fit[0:2]
    x_data, y_data = tools.tab_interval(
        constants.squared_temperature, constants.C_div_T, a, b)
    plt.figure()
    plt.plot()
    plt.plot()
    plt.grid(True)
    plt.show()


def main():
    print(linear_fit(0, 10, constants.E_curve_fit, constants.n_curve_fit))


if __name__ == "__main__":
    main()
