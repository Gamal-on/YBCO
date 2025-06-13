from data import err_temperature
from data import err_sample_HC
from data import sample_HC
from data import temperature
import numpy as np
import matplotlib.pyplot as plt
import tools


# Constants

k = 1.380649e-23
alpha = 2.399357318878174


squared_temperature = temperature**2  # K**2
C_div_T = sample_HC/temperature  # mJ/K**2.mol
err_C_divT = err_sample_HC/temperature

# Entropy : integration


def entropie(a, b):
    return tools.integrate(temperature, C_div_T, a, b)


def main() :
    entropie(0, 100)

if __name__ == "__main__":
    main()