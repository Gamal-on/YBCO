import numpy as np
import matplotlib.pyplot as plt
import fitutils as ft
from tools import tab_interval, maximum, schottky, plot

# Constants

k = 1.380649e-23
alpha = 2.399357318878174

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature

squared_temperature = temperature**2 #K**2
C_div_T = sample_HC/temperature #mJ/K**2.mol
err_C_divT = err_sample_HC/temperature

# Extremum Schottky anomaly

def max_schottky(x, y, min, max) :
    x_interval, y_interval = tab_interval(x, y, min, max)
    maxi, i = maximum(y_interval)
    x_maxi = x_interval[i]
    return x_maxi, maxi

def energie(x, y, min, max) :
    x_maxi, maxi = max_schottky(x, y, min, max)
    return alpha*np.sqrt(x_maxi)*k

# C/T without Schottky contribution

E = energie(squared_temperature, C_div_T, 0, 12)
schottky_anomaly = schottky(temperature, E, n=3e24)
schottky_CdivT = C_div_T - schottky_anomaly

# Trouver n = nombre d'atomes

# Linear fit

def ajustement(x, y, err_x, err_y) :
    fit = ft.linfitxy(x, y, err_x, err_y, plot=True)
    plt.show()
    return fit

def main() :
    plot

if __name__ == "__main__":
    main()