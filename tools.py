import numpy as np
import matplotlib.pyplot as plt

# Function to choose a certain interval


def tab_interval(x, y, min, max):
    """Choose bounds in tab1, return the corresponding values in tab 1 and tab 2"""
    idx_min = 0
    idx_max = 0
    tab = []
    for valeur in x:
        if valeur < min:
            idx_min += 1
            idx_max += 1
        elif valeur == min:
            tab.append(valeur)
            idx_max += 1
        elif min < valeur < max:
            tab.append(valeur)
            idx_max += 1
        elif valeur == max:
            tab.append(valeur)
            idx_max += 1
    y_interval = y[idx_min: idx_max]
    x_interval = np.array(tab)
    return x_interval, y_interval

# Find maximum


def maximum(tab):
    idx_max = 0
    max = tab[0]
    for i in range(0, len(tab)):
        if tab[i] > max:
            max = tab[i]
            idx_max = i
    return max, idx_max


test = np.array([2, 56, 2, 45, 234, 89, 6])

# Schottky anomaly


def schottky(T, E, n=1, k=1.380649e-23):
    x = (E)/(k*T)
    cs = k*(x**2)*(np.exp(x)/(1+np.exp(x))**2)
    return n*cs

# Plot


def plot(x, y, min, max):
    x_interval, y_interval = tab_interval(x, y, min, max)
    plt.figure()
    plt.plot(x_interval, y_interval, ".g")
    plt.grid(True)
    plt.show()

# Newton Raphson


def resol_nr(f, fp, N, xstart, eps):
    i = 0
    x2 = xstart
    condition = True
    while condition:
        i += 1
        x1 = x2
        x2 = x2-f(x2)/fp(x2)
        condition = np.abs(x2-x1) > eps and N < i
    return x2, f(x2)

# Dichotomie


def dicho_1(f, eps, a, b):
    borne_min = a
    borne_max = b
    milieu = (b-a)/2
    res = f(milieu)
    while np.abs(res) > eps:
        resmin = f(borne_min)
        if resmin * res > 0:
            borne_min = milieu
        else:
            borne_max = milieu
        milieu = (borne_min + borne_max)/2
        res = f(milieu)
    return res, milieu

# Integrate


def simpson(f, a, b, n):
    somme = 0
    dx = (b-a)/n
    for i in range(0, n):
        x = a + i*dx
        x1 = x + dx
        somme += (f(x) + 4*f((x + x1)/2) + f(x1))
    return (b-a)/(6*n) * somme


def main():
    print(maximum(test))


if __name__ == "__main__":
    main()
