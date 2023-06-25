import csv
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

def generate_curve(data):
  # Calcola l'istogramma
  hist, bin_edges = np.histogram(data, bins='auto', density=True)

  # Calcola il centro di ogni intervallo di bin
  bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
  print("bin_centers", bin_centers)

  # Normalizza le frequenze
  hist_norm = hist / np.sum(hist)

  # Definisci la funzione di densità di probabilità della distribuzione normale
  def norm_pdf(x, mu, sigma):
      return stats.norm.pdf(x, loc=mu, scale=sigma)

  # Stimazione dei parametri utilizzando la regressione dei minimi quadrati
  params, _ = curve_fit(norm_pdf, bin_centers, hist_norm, method='lm')

  # Estrai i parametri stimati
  mu_estimato = params[0]
  sigma_estimato = params[1]

  # Genera la curva della distribuzione normale prevista
  x = np.linspace(min(bin_centers), max(bin_centers), len(data))
  y_pred = norm_pdf(x, mu_estimato, sigma_estimato)
  return x, y_pred,mu_estimato,sigma_estimato


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
#    x, y,mu_estimato,sigma_estimato = generate_curve(valori_colonne[colonna])
    plt.hist(valori_colonne[colonna], bins='auto', alpha=0.7, density=True)
#    plt.plot(x, y, 'r-', label='Distribuzione Normale')

    xmin, xmax = plt.xlim()
    print("xmin, xmax: ", xmin, xmax)
    xlen = 100 #len(valori_colonne[colonna])
    # Genera un array di valori x per la curva della distribuzione normale
    x = np.linspace(xmin, xmax, xlen)
    # Calcola i valori della funzione di densità di probabilità della distribuzione normale
    pdf = stats.norm.pdf(x, loc=media, scale=np.sqrt(varianza))
    # Traccia la curva della distribuzione normale
    plt.plot(x, pdf, 'r', label='Distribuzione Normale')


    # Aggiunta delle linee verticali per i parametri
    plt.axvline(media, color='r', linestyle='dashed', linewidth=2, label='Media')
    plt.axvline(median, color='g', linestyle='dashed', linewidth=2, label='Mediana')
    plt.legend()

    # Aggiunta del testo per i parametri
    testo = f"{colonna}\nMedia: {media:.2f}\nVarianza: {varianza:.2f}\nSkewness: {distribution.skewness:.2f}"
    plt.text(0.98, 0.80, testo, ha='right', va='top', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

    # Test di Kolmogorov-Smirnov
    D, p_value = stats.kstest(valori_colonne[colonna], 'norm', args=(media, np.sqrt(varianza)))
    print(f"Test di Kolmogorov-Smirnov: D={D}, p-value={p_value}")
    D, p_value = stats.kstest((valori_colonne[colonna]-media)/np.sqrt(varianza), 'norm')
    print(f"Test di Kolmogorov-Smirnov: D={D}, p-value={p_value}")


    # Mostra il grafico
    plt.show()

