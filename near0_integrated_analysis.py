import numpy as np
import constants as cnt
from scipy.integrate import quad
import scipy.optimize as sc
import matplotlib.pyplot as plt
import tools
import schottky_analysis as sch

# Fonction à intégrer :


def integrand(x):
    """Fonction à intégrer"""
    exp_x = np.exp(x)
    denom = (exp_x - 1)**2
    return (x**4 * exp_x) / denom

# Calcul de l'intégrale I(y)


def debye_integral(y):
    """Calcule l'intégrale de 0 à y de la fonction integrand(x)."""
    result, _ = quad(integrand, 0, y)
    return result


# Vecteurisation de la fonction I pour traiter les tableaux
I_vec = np.vectorize(debye_integral)

# Modèle complet


def model_integral_schottky(x, gamma, theta, E, n):
    """
    Modèle complet:
    Paramètres:
    x      : variable indépendante, T²
    Cs     : contribution de schottky
    gamma  : paramètre constant
    theta  : paramètre de l'intégrale
    """
    y = theta / np.sqrt(x)  # Calcul de y
    I_vals = I_vec(y)       # Calcul de l'intégrale
    term = 9*cnt.N*cnt.k * (x / (theta**3)) * I_vals
    return sch.schottky(np.sqrt(x), E, n) + gamma + term


def model_integral_schottky_substracted(x, gamma, theta):
    """
    Modèle linéaire :

    Paramètres:
    x      : variable indépendante, T²
    gamma  : paramètre constant
    theta  : paramètre de l'intégrale
    """
    y = theta / np.sqrt(x)  # Calcul de y
    I_vals = I_vec(y)       # Calcul de l'intégrale
    term = 9*cnt.N*cnt.k * (x / (theta**3)) * I_vals
    return gamma + term
    return

# Fonction wrapper pour curve fit


def model_integral_schottky_wrap(x, gamma, theta, E, n):
    return model_integral_schottky(x, gamma, theta, E, n)


def model_integral_schottky_substracted_wrap(x, gamma, theta):
    return model_integral_schottky_substracted(x, gamma, theta)

# Fit avec curve fit


def fit_integral_schottky(a, b, x, y):
    """Warning : squared temperature is used"""
    x_bounded, y_bounded = tools.tab_interval(x, y, a, b)
    fit = sc.curve_fit(model_integral_schottky_wrap, x_bounded, y_bounded)
    return fit[0]


def fit_integral_schottky_substracted(a, b, x, x_carre, y):
    """Warning : squared temperature is used"""
    mask = (x >= a) & (x <= b)
    x_carre_bounded, y_bounded = tools.tab_interval(x_carre, y, a**2, b**2)
    x_bounded = x[mask]
    y_schottky_subs_bounded = y_bounded - \
        sch.schottky(x_bounded, cnt.E_optic, n=cnt.n_optic)
    fit = sc.curve_fit(model_integral_schottky_substracted_wrap,
                       x_carre_bounded, y_schottky_subs_bounded, bounds=([0, 300], [10, 600]))
    return fit[0]


def main():
    print(fit_integral_schottky_substracted(
        0, 20, cnt.temperature, cnt.squared_temperature, cnt.C_div_T))


if __name__ == "__main__":
    main()
