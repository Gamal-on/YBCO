import numpy as np
import constants as cnt
import matplotlib.pyplot as plt
import tools
import schottky_analysis as sch
from functools import lru_cache
from scipy.integrate import quad
from scipy.optimize import curve_fit

# Définition de l'intégrande


def integrand(t):
    # t scalar
    return t**4 * np.exp(t) / (np.expm1(t)**2)

# Fonction pour I(y) avec mise en cache


@lru_cache(maxsize=1000)
def I_of_y_scalar(y):
    # y est un scalaire Python float >= 0
    val, err = quad(lambda t: integrand(t), 0, y, epsabs=1e-8, epsrel=1e-8)
    return val


def I_of_y(y_array):
    """
    Calcule I(y) pour un tableau ou scalaire y_array.
    Utilise la fonction scalaires mise en cache pour accélérer.
    """
    y_arr = np.atleast_1d(y_array)
    result = np.empty_like(y_arr, dtype=float)
    for i, yi in enumerate(y_arr):
        if yi < 0:
            raise ValueError(f"y must be >=0, got {yi}")
        result[i] = I_of_y_scalar(float(yi))
    if np.isscalar(y_array):
        return result.item()
    return result


def model(x, theta, gamma):
    """
    Modèle f(x) = gamma + x * I(theta / x)
    x: array-like positif
    theta, gamma: scalaires
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = theta / np.sqrt(x_arr)
    I_vals = I_of_y(y_arr)
    return gamma + x_arr * I_vals * 9*cnt.k*cnt.N/(theta**3)


def fit_data(xdata, ydata, p0=None, bounds=None):
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
    # Appel curve fit
    popt, pcov = curve_fit(model, x_arr, y_arr, p0=p0,
                           bounds=bounds or (-np.inf, np.inf))
    return popt, pcov


def main():
    squared_temperature_bounded, C_div_T_bounded = tools.tab_interval(
        cnt.squared_temperature, cnt.C_div_T, 0, 400)
    temperature_bounded, err_sample_bounded = tools.tab_interval(
        cnt.temperature, cnt.err_sample_HC, 0, 20)
    ydata = C_div_T_bounded - \
        sch.schottky(temperature_bounded, cnt.E_curve_fit,
                     cnt.n_curve_fit)/temperature_bounded
    print(fit_data(squared_temperature_bounded,
          ydata, p0=[300,-1] , bounds=([300, -1], [600, 10]))[0])


if __name__ == "__main__":
    main()
