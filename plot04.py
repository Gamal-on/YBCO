import numpy as np
import matplotlib.pyplot as plt
import tools

from data import temperature
from data import sample_HC
from data import err_sample_HC
from data import err_temperature

squared_temperature = temperature**2  # K**2
C_div_T = sample_HC/temperature  # mJ/K**2.mol
err_C_divT = sample_HC*((err_temperature/temperature) **
                        2 + (err_sample_HC/sample_HC)**2)
err_squared_temperature = (temperature)*2*err_temperature

# Plot of HC(T)


def plot_hc_vs_T():
    plt.figure()
    plt.plot(temperature, sample_HC, ".g", label="experimental data")
    plt.errorbar(temperature, sample_HC, err_sample_HC, err_temperature, "+g")
    plt.grid(True)
    plt.xlabel("T (K)")
    plt.ylabel("C (mJ/K.mol)")
    plt.title("Heat capacity function of temperature of YBCO optimally doped")
    plt.legend()
    plt.show()


# Plot of HC/T (T**2)

def plot_hc_vs_Tsquared():
    plt.figure()
    plt.plot(squared_temperature, C_div_T, ".g", label="experimental data")
    plt.errorbar(squared_temperature, C_div_T, err_C_divT,
                 err_squared_temperature, "+g")
    plt.grid(True)
    plt.xlabel("T² (K)")
    plt.ylabel("C/T (mJ/K².mol)")
    plt.legend()
    plt.show()

# Plot of HC/T (T²) near 0


def plot_hc_near0():
    a = int(input("Input the lower bound of T² : "))
    b = int(input("Input the higher bound of T² : "))
    x, y = tools.tab_interval(squared_temperature, C_div_T, a, b)
    x_err, y_err = err_squared_temperature[0:len(x)], err_C_divT[0:len(y)]
    plt.figure()
    plt.plot(x, y, ".g", label="experimental data")
    plt.errorbar(x, y, y_err, x_err, "+g")
    plt.grid(True)
    plt.xlabel("T² (K)")
    plt.ylabel("C/T (mJ/K.mol)")
    plt.title("C/T near 0 of YBCO optimally doped")
    plt.legend()
    plt.show()


def main():
    plot_hc_near0()


if __name__ == "__main__":
    main()
