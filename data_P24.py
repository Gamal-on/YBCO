import numpy as np
import matplotlib.pyplot as plt


# Constants

mY = 88.905
mBa = 137.3270
mCu = 63.5460
mO = 15.9994
masse_molaire = 666.22  # g.mol
masse_sample = 6.81e-3  # g

# File reading

data_file = open("fork_P240917.txt", "r",  encoding="latin-1")
lignes = data_file.readlines()
data_file.close()

idx_data = 12  # line of "data"

# Header line after "Data"
header_line = lignes[idx_data + 1].strip()
noms_colonnes = [c.strip() for c in header_line.split(",")]
donnees_lignes = lignes[idx_data + 2:]

# Dictionnary empty with lists for colonns
colonnes = {nom: [] for nom in noms_colonnes}

# Fill the dictionnary
for ligne in donnees_lignes:
    if not ligne.strip():
        continue  # ignorer lignes vides
    champs = ligne.strip().split(",")
    if len(champs) != len(noms_colonnes):
        continue  # ignorer lignes incomplètes

    for nom, valeur in zip(noms_colonnes, champs):
        try:
            colonnes[nom].append(float(valeur))
        except ValueError:
            colonnes[nom].append(np.nan)  # ou None si tu préfères

# Arrays:

temperature_P24 = np.array(colonnes["Sample Temp (Kelvin)"])
# final heat capacity (without the addenda), mJ/K.mol
sample_HC_P24 = np.array(colonnes["Samp HC (mJ/mole-K)"])

# np.array(colonnes["Samp HC Err (mJ/mole-K)"])  # error on final HC
err_sample_HC_P24 = 0.02*sample_HC_P24
# no error in the data concerninf temperature
err_temperature_P24 = 1e-2*np.ones(len(err_sample_HC_P24))


def main():
    pass


if __name__ == "__main__":
    main()
