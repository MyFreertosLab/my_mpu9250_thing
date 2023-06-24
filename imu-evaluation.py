import csv
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def leggi_dati_da_csv(file_path):
    dati = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for riga in reader:
            dati.append({campo: float(riga[campo]) for campo in reader.fieldnames})
    return dati

def calcola_media_varianza_colonne(dati):
    media_colonne = {}
    varianza_colonne = {}
    median_colonne = {}
    distribution_colonne = {}
    valori_colonne = {}
    
    colonne_speciali = ['MX', 'MY', 'MZ']
    colonne_normali = [colonna for colonna in dati[0].keys() if colonna not in colonne_speciali]
    
    for colonna in colonne_normali:
        valori_colonna = [riga[colonna] for riga in dati]
        media_colonna = np.mean(valori_colonna)
        varianza_colonna = np.var(valori_colonna)
        median_colonna = np.median(valori_colonna)
        distribution_colonna = stats.describe(valori_colonna)
        
        media_colonne[colonna] = media_colonna
        varianza_colonne[colonna] = varianza_colonna
        median_colonne[colonna] = median_colonna
        distribution_colonne[colonna] = distribution_colonna
        valori_colonne[colonna] = valori_colonna
    
    for colonna in colonne_speciali:
        valori_colonna = [riga[colonna] for riga in dati if riga['MV'] == 1]
        media_colonna = np.mean(valori_colonna)
        varianza_colonna = np.var(valori_colonna)
        median_colonna = np.median(valori_colonna)
        distribution_colonna = stats.describe(valori_colonna)
        
        media_colonne[colonna] = media_colonna
        varianza_colonne[colonna] = varianza_colonna
        median_colonne[colonna] = median_colonna
        distribution_colonne[colonna] = distribution_colonna
        valori_colonne[colonna] = valori_colonna
    
    return media_colonne, varianza_colonne, median_colonne, distribution_colonne, valori_colonne

file_csv = "examples/01-imu-raw-data.csv"
dati = leggi_dati_da_csv(file_csv)
media_colonne, varianza_colonne, median_colonne, distribution_colonne, valori_colonne = calcola_media_varianza_colonne(dati)

print("valori[AZ]: ", valori_colonne["AZ"])

#print("Media e varianza delle colonne:")
for colonna, media in media_colonne.items():
    varianza = varianza_colonne[colonna]
    median = median_colonne[colonna]
    distribution = distribution_colonne[colonna]
    print(f"{colonna}: Media={media}, Varianza={varianza}, Mediana={median}, Distribution={distribution}")
    plt.hist(valori_colonne[colonna], bins='auto', alpha=0.7)

    # Aggiunta delle linee verticali per i parametri
    plt.axvline(media, color='r', linestyle='dashed', linewidth=2, label='Media')
    plt.axvline(median, color='g', linestyle='dashed', linewidth=2, label='Mediana')
    plt.legend()

    # Aggiunta del testo per i parametri
    testo = f"{colonna}\nMedia: {media:.2f}\nVarianza: {varianza:.2f}\nSkewness: {distribution.skewness:.2f}"
    plt.text(0.98, 0.80, testo, ha='right', va='top', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

    # Mostra il grafico
    plt.show()

