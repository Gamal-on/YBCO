import numpy as np
import matplotlib.pyplot as plt

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature

squared_temperature = temperature**2 #K**2
C_div_T = sample_HC/temperature #mJ/K**2.mol
err_C_divT = err_sample_HC/temperature

# Plot of HC/T (T**2)

def plot_HC_divT() :
    plt.figure()
    plt.plot(squared_temperature, C_div_T, ".g")
    plt.errorbar(squared_temperature, C_div_T, err_C_divT, err_temperature, "g+")
    plt.ylabel("C/T (mJ/K².mol)")
    plt.xlabel("T² (K²)")
    plt.show()

# Plot of HC/T (T**2) near 0

def tab_near0() :
    """Return HC/T and T² arrays with only the values T<20 K"""
    temp_near0 = []
    idx_near0 = 0
    for valeur in squared_temperature : 
        if valeur < 400 : 
            temp_near0.append(valeur)
            idx_near0 += 1
    C_div_T_near0 = C_div_T[0:idx_near0]
    squared_temp_near0 = np.array(temp_near0)
    err_C_divT_near0 = err_C_divT[0:idx_near0]
    err_temp_near0 = err_temperature[0:idx_near0]
    return squared_temp_near0, C_div_T_near0, err_C_divT_near0, err_temp_near0

def plot_near0() :
    """Plot of HC/T(T²) near 0"""
    x, y, err_y, err_x = tab_near0()
    plt.plot(x, y, "g.")
    plt.errorbar(x, y, err_y, err_x, fmt="g+")
    plt.show()

def tab_interval(a, b) :
    """Return HC/T and T² arrays with only the values a<T<b K"""
    min = a**2
    max = b**2
    idx_min = 0
    idx_max = 0
    tab = []
    for valeur in squared_temperature : 
        if valeur < min : 
            idx_min += 1
            idx_max += 1
        elif valeur == min :
            tab.append(valeur)
            idx_max += 1
        elif min < valeur < max : 
            tab.append(valeur)
            idx_max += 1
        elif valeur == max :
            tab.append(valeur)
            idx_max += 1
    C_div_T_interval = C_div_T[idx_min : idx_max]
    squared_temperature_interval = np.array(tab)
    err_C_divT_interval = err_C_divT[idx_min : idx_max]
    err_squared_temperature_interval = err_temperature[idx_min : idx_max]
    #print(len(C_div_T_interval), len(squared_temperature_interval))
    return squared_temperature_interval, C_div_T_interval, err_C_divT_interval, err_squared_temperature_interval

def plot() :
    a = float(input("Inferior bound : "))
    b = float(input("Superior bound : "))
    x, y, err_y, err_x = tab_interval(a,b)
    plt.plot(x, y, ".r")
    plt.errorbar(x, y, err_y, err_x, "g+")
    plt.grid(True)
    plt.xlabel("T² (K²)")
    plt.ylabel("C/T (mJ/K².mol)")
    plt.show()

def main() :
    plot()

if __name__ == "__main__":
    main()