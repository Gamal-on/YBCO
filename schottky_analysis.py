import numpy as np
import matplotlib.pyplot as plt
import tools
import near0_nonlinear_acoustic
import constants as cnt



beta, gamma, n, E = cnt.beta_optic, cnt.gamma_quadratic, cnt.n_exp, cnt.E_exp_8


def schottky(T, E, n):
    """Calculate the Schottky anomaly in mJ/K.mol"""
    """T: temperature in Kelvin, E: energy in Joules, n: number of particles, k: Boltzmann constant"""
    x = (E)/(cnt.k*T)
    cs = (x**2)*(np.exp(x)/(1+np.exp(x))**2)
    return n*cnt.r*1e3*cs


def dev_schottky(T, E, n):
    """Calculate the derivative of the Schottky anomaly"""
    """T: temperature in Kelvin, E: energy in Joules, n: number of particles"""
    a = E / cnt.k
    exp_at = np.exp(a / T)
    num = (2 * T + a) * exp_at + (2 * T - a) * (exp_at ** 2)
    denom = (T ** 4) * (1 + exp_at) ** 3
    return - n * cnt.k * (a ** 2) * num / denom


# Relation entre T_max (température à laquelle schottky est maximale) et E (delta d'énergie entre les deux niveaux du système)
# en recherchant les zeros de la dérivée : méthode dichotomie et Newton-Raphson
# On trouve un alpha = 2.5127

def f(x):
    """Function to find the zero in the derivative of the Schottky anomaly"""
    return 2+x+(2-x)*np.exp(x)


def fp(x):
    """Derivative of the function to find the zero in the derivative of the Schottky anomaly"""
    return 1 + 2*np.exp(x) - np.exp(x) - x*np.exp(x)


def alpha():
    """Calculate the alpha value as k*T_max*alpha = E"""
    a = tools.dicho_1(f, 1e-5, 0, 10)[1]
    b = tools.resol_nr(f, fp, 1e5, 2, 1e-5)[0]
    return b

# Recherche de T_max dans nos données expérimentales, comme étant la temérature à laquelle la dérivée de la Schottky est nulle.
# Intervalle de recherche entre 0 et 10 K²
# T_max trouvée : 2.9461005 K


def max_schottky(x, y, min, max):
    """Find the maximum of the Schottky anomaly in a given interval.
    x: array of temperatures in Kelvin, y: array of Schottky anomaly values, min: minimum temperature, max: maximum temperature
    Returns the temperature at which the maximum occurs and the maximum value."""
    x_interval, y_interval = tools.tab_interval(x, y, min, max)
    maxi, i = tools.maximum(y_interval)
    x_maxi = x_interval[i]
    return x_maxi, maxi

# n parameter determination

def n_det(x, y, E):
    w = E/(cnt.k*x)
    a = y*1e-3/cnt.r
    b = 1/(x*w**2)
    c = ((np.exp(w) + 1)**2)/np.exp(w)
    return a*b*c


def n_experimental(x, y, E):
    n_values = []
    for i in range(0, len(x)):
        n_values.append(n_det(x[i], y[i], E))
    return np.mean(n_values)

# Results


def results_shottky():
    T_max = max_schottky(cnt.temperature, cnt.C_div_T, 0, 3)[0]
    E_exp = cnt.k*alpha()*max_schottky(cnt.temperature, cnt.C_div_T, 0, 3)[0]
    n_exp = n_experimental(cnt.temperature[15:25], cnt.C_div_T[15:25], E_exp)
    return T_max, E_exp, n_exp


def main():
    pass


if __name__ == "__main__":
    main()
