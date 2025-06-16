import numpy as np
import matplotlib.pyplot as plt
import tools

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

# Functions


def schottky(T, E, n=1e-2):
    """Calculate the Schottky anomaly in mJ/K.mol"""
    """T: temperature in Kelvin, E: energy in Joules, n: number of particles, k: Boltzmann constant"""
    x = (E)/(k*T)
    cs = (x**2)*(np.exp(x)/(1+np.exp(x))**2)
    return n*r*1e3*cs


def dev_schottky(T, E, n=1):
    """Calculate the derivative of the Schottky anomaly"""
    """T: temperature in Kelvin, E: energy in Joules, n: number of particles"""
    a = E / k
    exp_at = np.exp(a / T)
    num = (2 * T + a) * exp_at + (2 * T - a) * (exp_at ** 2)
    denom = (T ** 4) * (1 + exp_at) ** 3
    return - n * k * (a ** 2) * num / denom


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
    # b = tools.resol_nr(f, fp, 100, 2, 1e-5)[0]
    return a

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

# Plot Schottky théorique aec les paramètres expérimentaux
# n pris aléatoirement entre 1e-3 et 1e-2


def plot_schottky(T):
    n = float(input("Enter the number of defaults : "))
    N = int(input("Enter the number of values to plot : "))
    E_exp = k*alpha()*max_schottky(temperature, C_div_T, 0, 3)[0]
    plt.figure()
    plt.plot(T[0:N], schottky(temperature, E_exp, n)[0:N], "g--")
    plt.ylabel("C/T (mJ/K².mol)")
    plt.xlabel("T (K) (or T²)")
    plt.title("Schottky anomaly, plot with experimental parameters")
    plt.grid(True)
    plt.show()

# n parameter determination


def n_det(x, y, E):
    w = E/(k*x)
    a = y*1e-3/r
    b = 1/(x*w**2)
    c = ((np.exp(w) + 1)**2)/np.exp(w)
    return a*b*c


def n_experimental(x, y, E):
    n_values = []
    for i in range(0, len(x)):
        n_values.append(n_det(x[i], y[i], E))
    return np.mean(n_values)

# Results


T_max = max_schottky(temperature, C_div_T, 0, 3)[0]
E_exp = k*alpha()*max_schottky(temperature, C_div_T, 0, 3)[0]
n_exp = n_experimental(temperature[15:25], C_div_T[15:25], E_exp)


def main():
    print(E_exp)


if __name__ == "__main__":
    main()
