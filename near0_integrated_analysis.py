import numpy as np
import constants as cnt
import matplotlib.pyplot as plt
import tools
import schottky_analysis as sch
from functools import lru_cache
from scipy.integrate import quad
from scipy.optimize import curve_fit


def integrand(t):
    # t scalar
    return t**4 * np.exp(t) / ((np.exp(t)-1)**2)

# Fonction pour I(y) avec mise en cache


@lru_cache(maxsize=1000)
def I_of_y_scalar(y):
    val, err = quad(integrand, 0, y)
    return val


def I(y_array):
    result = []
    for yi in y_array:
        result.append(I_of_y_scalar(float(yi)))
    return result


def model_integrated_schottky_substracted(x, theta, gamma):
    """
    Modèle f(x) = gamma + x * I(theta / x)
    x: array-like positif
    theta, gamma: scalaires
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = theta / np.sqrt(x_arr)
    I_vals = I(y_arr)
    return gamma + x_arr * I_vals * 9e3*cnt.k*cnt.N/(theta**3)


def model_integrated_schottky(x, theta, gamma, E, n):
    """
    Modèle f(x) = gamma + x * I(theta / x) + Cs(sqrt(x))
    x: array-like positif
    theta, gamma: scalaires
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = theta / np.sqrt(x_arr)
    I_vals = I(y_arr)
    return gamma + x_arr * I_vals * 9e3*cnt.k*cnt.N/(theta**3) + sch.schottky(np.sqrt(x), E, n)/np.sqrt(x)


def fit_data(xdata, ydata, model, p0=None, bounds=None):
    """
    Ajuste les données (xdata, ydata) au modèle, retourne (popt, pcov).
    - xdata, ydata: array-like de même longueur.
    - p0: liste ou tuple [theta0, gamma0] d'initialisation. Par défaut [1.0, np.mean(ydata)].
    - bounds: bornes pour θ et γ : ([theta_min, gamma_min], [theta_max, gamma_max]). Par défaut aucun.

    Exemple d'utilisation :
        popt, pcov = fit_data(xdata, ydata, p0=[2.0, 0.5], bounds=([0, -np.inf], [np.inf, np.inf]))
    """
    x_arr = np.asarray(xdata, dtype=float)
    y_arr = np.asarray(ydata, dtype=float)
    if x_arr.shape != y_arr.shape:
        raise ValueError("xdata and ydata must have same shape.")
    # Initial guess
    if p0 is None:
        theta0 = 1.0
        gamma0 = np.mean(y_arr)
        p0 = [theta0, gamma0]
    popt, pcov = curve_fit(model, x_arr, y_arr, p0=p0,
                           bounds=bounds or (-np.inf, np.inf))
    return popt, pcov


def fit_integrand_schottky_subtracted(a, b, x_carre, x, y, E, n, p0, bounds):
    """Return the results of optimization of C/T(T²), linear form (Schottky contribution substracted)
    a, b : bounds (not squared)
    x_carre : T²
    x : T
    y : C/T
    E, n : Schottky parameters"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a**2, b**2)
    x_interval, x_interval = tools.tab_interval(x, x, a, b)
    y_data = y_interval - sch.schottky(x_interval, E, n)/x_interval
    fit = fit_data(x_carre_interval, y_data,
                   model_integrated_schottky_substracted, p0=p0, bounds=bounds)
    return fit


def fit_integrand_schottky(a, b, x_carre, y, p0, bounds):
    """Return the results of optimization of C/T(T²), Schottky contribution
    a, b : bounds (not squared)
    x_carre : T²
    x : T
    y : C/T"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a**2, b**2)
    fit = fit_data(x_carre_interval, y_interval,
                   model_integrated_schottky, p0=p0, bounds=bounds)
    return fit


def main():
    print(fit_integrand_schottky_subtracted(0, 20, cnt.squared_temperature, cnt.temperature,
          cnt.C_div_T, cnt.E_curve_fit, cnt.n_curve_fit, p0=[300, 0], bounds=([300, 0], [500, 10])))
    print(fit_integrand_schottky(0, 20, cnt.squared_temperature, cnt.C_div_T, p0=[
          300, 0, 9.7e-23, 1e-23], bounds=([300, 0, 9.7e-23, 1e-3], [500, 10, 1.2e-22, 5e-2])))


if __name__ == "__main__":
    main()
