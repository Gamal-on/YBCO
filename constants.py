import numpy as np
import tools

# Constants

k = 1.380649e-23
r = 8.31446261815324
N = 78.286e23

# Data of the HPHT sample

from data_HPHT import err_temperature_HPHT
from data_HPHT import err_sample_HC_HPHT
from data_HPHT import sample_HC_HPHT
from data_HPHT import temperature_HPHT

squared_temperature_HPHT = temperature_HPHT**2  # K**2
hc_div_temp_HPHT = sample_HC_HPHT/temperature_HPHT  # mJ/K**2.mol
err_hc_div_temp_HPHT = err_sample_HC_HPHT/temperature_HPHT
err_squared_temperature_HPHT = 2*temperature_HPHT*err_temperature_HPHT

# Fitted parameters for the HPHT sample

# Determined by curve fit directly
n_curve_fit_HPHT = 9.997e-3
E_curve_fit_HPHT = 1.0966e-22  # J

# Deduced form Schottky analysis
T_exp_HPHT = 2.9461005
E_exp_8_HPHT = 1.0681622692127778e-22  # with NR and 8.67
E_exp_10_HPHT = 1.1549944696942883e-22  # with mean and 10.148
n_exp_HPHT = 0.0016374152345313721

# Determined by optic fit
beta_optic_HPHT = 0.44633863261624984
gamma_optic_HPHT = 1.9612364149478514
n_optic_HPHT = 0.010999999999999998
E_optic_HPHT = 1.1183498428073282e-22
alpha_optic_HPHT = 0.0005023827216230027


# Data of ISIS sample

from data_ISIS import err_temperature_ISIS
from data_ISIS import err_sample_HC_ISIS
from data_ISIS import sample_HC_ISIS
from data_ISIS import temperature_P24

squared_temperature_ISIS = temperature_P24**2  # K**2
hc_div_temp_ISIS = sample_HC_ISIS/temperature_P24  # mJ/K**2.mol
err_hc_div_temp_ISIS = err_sample_HC_ISIS/temperature_P24
err_squared_temperature_ISIS = 2*temperature_P24*err_temperature_ISIS

# Data of P240 sample

from data_P24 import temperature_P24
from data_P24 import err_temperature_P24
from data_P24 import sample_HC_P24
from data_P24 import err_sample_HC_P24

squared_temperature_P24 = temperature_P24**2  # K**2
hc_div_temp_P24 = sample_HC_ISIS/temperature_P24  # mJ/K**2.mol
err_hc_div_temp_P24 = err_sample_HC_ISIS/temperature_P24
err_squared_temperature_P24 = 2*temperature_P24*err_temperature_ISIS


def main():
    pass


if __name__ == "__main__":
    main()
