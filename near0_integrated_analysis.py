import numpy as np
import constants as cnt
from scipy.integrate import quad
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import tools

# 1. Définition de l'intégrande de Debye


def integrande_debye(t):
    """Fonction intégrande pour l'intégrale de Debye."""
    # Gestion numérique de la singularité en t=0
    if abs(t) < 1e-12:
        return 0.0  # Limite t->0 : t^4/(e^t-1)^2 * e^t → 0
    return t**4 * np.exp(t) / (np.exp(t) - 1)**2

# 2. Modèle de capacité thermique de Debye


def modele_debye(T, C_scale, Td):
    """
    Modèle de Debye pour la capacité thermique.

    Args:
        T : Température (ou array de températures)
        C_scale : Facteur d'échelle (9*N*k)
        Td : Température de Debye (paramètre d'ajustement)

    Returns:
        C_v : Capacité thermique prédite
    """
    def integree(x):
        """Calcule l'intégrale pour une valeur x = Td/T"""
        result, _ = quad(integrande_debye, 0, x)
        return result

    # Vectorisation pour gérer les arrays
    if isinstance(T, (int, float)):
        T = np.array([T])

    C_v = np.zeros_like(T)
    for i, temp in enumerate(T):
        x = Td / temp
        integrale = integree(x)
        C_v[i] = C_scale * (temp/Td)**3 * integrale

    return C_v


# 3. Chargement des données expérimentales
T_data, C_data = tools.tab_interval(cnt.squared_temperature, cnt.C_div_T, 0, 400) # Températures en Kelvin

# 4. Ajustement du modèle aux données
# Estimation initiale des paramètres [C_scale, Td]
p0 = [1e-3, 300]

# Bornes pour éviter des valeurs non physiques
# [C_scale_min, Td_min], [C_scale_max, Td_max]
bornes = ([1e-6, 50], [1e-1, 2000])

params_opt, covariance = curve_fit(
    modele_debye,
    T_data,
    C_data,
    p0=p0,
    bounds=bornes,
    maxfev=5000  # Augmenter le nombre d'évaluations si nécessaire
)

# 5. Extraction des paramètres optimisés
C_scale_opt, Td_opt = params_opt
print(f"Paramètres optimisés:", params_opt)
print(f"Facteur d'échelle (C_scale) = {C_scale_opt:.4e}")
print(f"Température de Debye (Td) = {Td_opt:.1f} K")
print(f"Rapport signal/bruit: {np.diag(covariance)}")

# 6. Visualisation des résultats
T_fit = np.linspace(min(T_data), max(T_data), 100)
C_fit = modele_debye(T_fit, C_scale_opt, Td_opt)

plt.figure(figsize=(10, 6))
plt.scatter(T_data, C_data, label='Données expérimentales', s=50)
plt.plot(T_fit, C_fit, 'r-', label=f'Modèle Debye: Td={Td_opt:.1f}K')
plt.xlabel('Température (K)', fontsize=12)
plt.ylabel('Capacité thermique $C_v$', fontsize=12)
plt.title('Ajustement au modèle de Debye', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
