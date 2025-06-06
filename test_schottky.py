import numpy as np
import matplotlib.pyplot as plt

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature
from tools import resol_nr, dicho_1

squared_temperature = temperature**2 #K**2
C_div_T = sample_HC/temperature #mJ/K**2.mol
err_C_divT = err_sample_HC/temperature

# Constantes et tableaux 

k = 1.380649e-23
delta = 2.9461005*k
temp = np.arange(0, 3, 1e-3)

def schottky(T, E, n=1) :
    x = (E)/(k*T)
    cs = k*(x**2)*(np.exp(x)/(1+np.exp(x))**2)
    return n*cs

def dev_schottky(T, E, n=1) :
    a = E / k
    exp_at = np.exp(a / T)        
    num = (2 * T + a) * exp_at + (2 * T - a) * (exp_at ** 2)
    denom = (T ** 4) * (1 + exp_at) ** 3
    return - n * k * (a ** 2) * num / denom

def plot_schottky(T, E, n=1) :
    plt.figure()
    plt.plot(T, schottky(T,E, n=1), "r")
    plt.grid(True)
    plt.show()

def plot_dev_schottky(T, E, n=1) :
    plt.figure()
    plt.plot(T, dev_schottky(T, E, n=1), "r")
    plt.grid(True)
    plt.show()

# Relation entre T* et E en recherchant les zeros de la dérivée : E = k*T*2.399357318878174, méthode dichotomie

def f(x) :
    return 2+x+(2-x)*np.exp(x)

def fp(x) :
    return 1 + 2*np.exp(x) - np.exp(x) - x*np.exp(x)

def main() :
    print(dicho_1(f, 1e-10, 2.3,2.5))
    

if __name__ == "__main__":
    main()

