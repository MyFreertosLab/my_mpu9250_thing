import csv
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

#print("Media e varianza delle colonne:")
for colonna, media in media_colonne.items():
    varianza = varianza_colonne[colonna]
    median = median_colonne[colonna]
    distribution = distribution_colonne[colonna]
    print(f"{colonna}: Media={media}, Varianza={varianza}, Mediana={median}, Distribution={distribution}")
    plt.hist(valori_colonne[colonna], bins='auto', alpha=0.7, density=True)

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

###############################
## Step 1): Form Tensor T 5xnxm
###############################
def t0(x, variance):
  return 1
def t1(x, variance):
  return x
def t2(x, variance):
  return x**2 - variance
def t3(x, variance):
  return x**3 - 3*x*variance
def t4(x, variance):
  return  x**4 - 6* (x**2) *variance + 3*(variance**2)
def terror(x, variance):
  assert false, "k must be in range(1:6)"
def T(k, i, l, X, variance):
  switch_cases = {
    1: t0,
    2: t1,
    3: t2,
    4: t3,
    5: t4
  }
  return switch_cases.get(k+1, terror)(X[l,i], variance)

###############################
## Step 2): Define index matrix
###############################
IM = np.array([[1,1],[1,2],[1,3],[2,2],[2,3],[3,3],[1,0],[2,0],[3,0],[0,0]])

#####################################
## Step 3): Form Tensor R 
#####################################
def is_equal(v1, v2):
    if v1 == v2:
        return 1
    else:
        return 0
def R(p,q,i, M):
  return is_equal(M[p,0], i) + is_equal(M[p,1], i)  + is_equal(M[q,0], i)  + is_equal(M[q,1], i) 

#####################################
## Step 4): Compute ni_als
#####################################
def ni_als(p,q,X,variance,M):
   m = X.shape[0]
   n = X.shape[1]
   s = 0
   for l in range(m):
     f = 1
     for i in range(n):
       f = f*T(R(p,q,i+1,M),i,l,X,variance)
     s = s + f
   return s

#####################################
## Step 5): Define index off-diagonal
#####################################
def index_off_diagonal(n=3):
  # only for n=3
  return np.array([1,2,4])


#####################################
## Step 6): form matrix psi_als
#####################################
def make_psi_als(X,variance,M):
  n = X.shape[1]
  n1 = M.shape[0]
  D = index_off_diagonal(n)
  psi_als = np.zeros((n1, n1)) 
  for p in range(n1):
    for q in range(p,n1):
      if p in D and q in D :
        psi_als[p,q] = 4*ni_als(p,q,X,variance, M)
      else:
        if p not in D and q not in D:
          psi_als[p,q] = ni_als(p,q,X,variance, M)
        else:
          psi_als[p,q] = 2*ni_als(p,q,X,variance, M)
      psi_als[q,p] = psi_als[p,q]
  return psi_als

################################################
## Step 7)8): Find eigenvector of min eigenvalue
################################################
## np.linalg.eig non restituisce gli eigenvalues
## in ordine decrescente (ipotesi dell'algoritmo)
## quindi è necessario ordinarli
def eigen(A):
    eigenValues, eigenVectors = np.linalg.eig(A)
    idx = np.argsort(eigenValues)[::-1][:A.shape[0]]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return (eigenValues, eigenVectors)

def calc_eigenvector_psi_als(psi_als):
    eigen_als = eigen(psi_als)
    eigen_values = eigen_als[0]
    eigen_vectors = eigen_als[1]
    idx = np.argmin(np.abs(eigen_values))
    beta = np.matrix(eigen_vectors[:, idx])
    result = beta / np.linalg.norm(beta, ord=2)
    return result

#def calc_eigenvector_psi_als(psi_als):
#  eigen_values, eigen_vectors = np.linalg.eig(psi_als)
#  idx = np.argmin(np.abs(eigen_values))
#  beta = np.delete(eigen_vectors, idx, axis=1)
#  result = beta/np.linalg.norm(beta)
#  return result

################################################
## Step 9): Estimate A, b, c, d
################################################
def estimate_A(beta, n):
    print("beta: ", beta, "n ", n)
    expanded_beta = beta[:(n*(n+1)//2)]
    print("expandend beta: ", expanded_beta)
    expanded_matrix = np.zeros((n, n))

    idx = 0
    for i in range(n):
        for j in range(i, n):
            expanded_matrix[i, j] = expanded_beta[0,idx]
            expanded_matrix[j, i] = expanded_beta[0,idx]
            idx += 1

    return expanded_matrix

def estimate_b(beta, n):
    start_index = int(n*(n+1)/2)
    end_index = beta.shape[1] - 1
    result=beta[0,start_index:end_index]
    return np.transpose(result)

def estimate_d(beta, n):
    return beta[0,-1]

def estimate_c(A, b, n):
    c = -0.5 * np.linalg.solve(A, b)
    return c

def estimate_Ae(A, c, d):
    denominator = np.dot(np.dot(np.transpose(c), A), c) - d
    scaling_factor = 1.0 / denominator
    print("scaling_factor: ", scaling_factor)
    Ae = scaling_factor.item() * A
    return Ae

def to_definite_positive_factor(A):
    eigenvalues = eigen(A)[0]
    result = 0
    if np.prod(eigenvalues >= 0) == 1:
        result = 1
    elif np.prod(eigenvalues < 0) == 1:
        result = -1

    return result

def estimate_all(psi_als, n):
    beta = calc_eigenvector_psi_als(psi_als)
    A = estimate_A(beta, n)
    # forza la matrice ad essere definita positiva se è negativa
    factor = to_definite_positive_factor(A)
    print("factor: ", factor)
    assert factor != 0, "La matrice non può essere indefinita"
    A = factor * A
    
    b = factor * estimate_b(beta, n)
    d = factor * estimate_d(beta, n)
    c = estimate_c(A, b, n)
    Ae = estimate_Ae(A, c, d)
    
    als = {}
    als['A'] = A
    als['b'] = b
    als['d'] = d
    als['c'] = c
    als['Ae'] = Ae
    
    return als

def mag_estimate(als, mag_data):
    mag_model_Q = als['A']
    mag_model_u = als['b']
    mag_model_k = als['d']
    mag_model_b = -0.5 * np.linalg.solve(mag_model_Q, mag_model_u)
    
    # Verifica il risultato
    eigen_result = eigen(mag_model_Q)
    D = np.diag(eigen_result[0])
    V = eigen_result[1]
    TA = np.dot(np.dot(V, D), np.transpose(V))
    assert np.allclose(mag_model_Q, TA, atol=1e-15), "Verifica fallita: mag_model_Q non corrisponde a TA"
    TD = np.dot(np.dot(np.transpose(V), mag_model_Q), V)
    assert np.allclose(D, TD, atol=1e-15), "Verifica fallita: D non corrisponde a TD"
    
    # Radice quadrata matrice Q
    appo = np.dot(mag_model_Q, np.transpose(mag_model_Q))
    appo_eigen = eigen(appo)
    mag_model_V = appo_eigen[1]
    mag_model_D = np.diag(eigen(mag_model_Q)[0])
    print("mag model u", mag_model_u)
    print("mag model k", mag_model_k)
    print("mag model b", mag_model_b)
    print("mag model V", mag_model_V)
    print("eigenvals mag model V", appo_eigen[0])
    print("mag model D", mag_model_D)
    print("mag model Q", mag_model_Q)
    mag_model_magnetic_norm = (np.dot(np.dot(np.transpose(mag_model_b), mag_model_Q), mag_model_b) - mag_model_k).item()
    print("mag_model_magnetic_norm", mag_model_magnetic_norm)
    P1 = np.dot(np.transpose(mag_model_V), mag_model_u)
    print("P1: ", P1)
    P2 = np.dot(np.linalg.inv(mag_model_D), np.dot(np.transpose(mag_model_V), mag_model_u))
    print("P2: ", P2)
    P3 = np.dot(np.transpose(P1), P2).item()
    print("P3: ", P3)
    mag_model_alpha = (-4 * mag_model_magnetic_norm / (4 * mag_model_k - P3))
    print("mag_model_alpha: ", mag_model_alpha)
    mag_model_B = np.dot(np.dot(mag_model_V, np.sqrt(mag_model_alpha * mag_model_D)) , np.transpose(mag_model_V))
    print("mag_model_B: ", mag_model_B)
    mag_model_inv_A = mag_model_B
    mag_model_A = np.linalg.inv(mag_model_inv_A)
    ## FIXME
    #f = 1 / np.sqrt(mag_model_magnetic_norm) # real ..
    f = 1.09 / np.sqrt(mag_model_magnetic_norm) # but works ..
    ##
    mag_model_scale_factors = np.zeros((3, 3))
    np.fill_diagonal(mag_model_scale_factors, f)
    result = {}
    result['Q'] = mag_model_Q
    result['b'] = mag_model_b
    result['k'] = mag_model_k
    result['offset'] = mag_model_b
    result['V'] = mag_model_V
    result['D'] = mag_model_D
    result['B'] = mag_model_B
    result['Hm2'] = mag_model_magnetic_norm
    result['alpha'] = mag_model_alpha
    result['invA'] = mag_model_inv_A
    result['data_source'] = mag_data
    result['scale_factors'] = mag_model_scale_factors
    result['als'] = als
    
    return result

def mag_apply_estimator(mag_model):
    mag_data = mag_model['data_source']
    mag_model_inv_A = mag_model['invA']
    mag_model_b = mag_model['b']
    mag_data_target = np.zeros(mag_data.shape)
    
    # Applica il modello ai dati
    for i in range(mag_data.shape[0]):
        x = np.transpose(mag_data[i, :] - np.transpose(mag_model_b))
        mag_data_target[i, :] = np.transpose(np.dot(mag_model['scale_factors'], np.dot(mag_model_inv_A, x)))
    
    return mag_data_target

def mag_plot_data(mag_data, sphere_radius=-1, title=""):
    fig = plt.figure(figsize=(8,6), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(mag_data[:, 0], mag_data[:, 1], mag_data[:, 2], c=mag_data[:, 2], cmap=None)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim([-1.005, 1.005])
    ax.set_ylim([-1.005, 1.005])
    ax.set_zlim([-1.005, 1.005])

    ax.set_title(title)
    plt.show()

imu_data_original = pd.read_csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-data-mag.csv')
# Axis North-East-Down
imu_data_NED = imu_data_original.rename(columns={'MX': 'MX2', 'MY': 'MX'})
imu_data_NED = imu_data_NED.rename(columns={'MX2': 'MY'})

imu_data_NED['MY'] = -imu_data_NED['MY']
imu_data_NED['MX'] = -imu_data_NED['MX']
imu_data_NED['AY'] = -imu_data_NED['AY']
imu_data_NED['AZ'] = -imu_data_NED['AZ']
imu_data_NED['GY'] = -imu_data_NED['GY']
imu_data_NED['GZ'] = -imu_data_NED['GZ']

imu_data_mag = imu_data_NED[['MX', 'MY', 'MZ']].loc[imu_data_NED['MV'] == 1]

#psi_als_mag = make_psi_als(np.matrix(imu_data_mag), (5.618457 ** 2), IM)
psi_als_mag = make_psi_als(np.matrix(imu_data_mag), (6.054304 ** 2), IM)
print("psi_als_mag: ", psi_als_mag)
als_mag = estimate_all(psi_als_mag, 3)
imu_mag_estimator = mag_estimate(als_mag, np.matrix(imu_data_mag))
imu_data_mag_estimated = mag_apply_estimator(imu_mag_estimator)

centroid = np.array([np.mean(imu_data_mag_estimated[:,0]), np.mean(imu_data_mag_estimated[:,1]),np.mean(imu_data_mag_estimated[:,2])])
imu_data_mag_centrate = imu_data_mag_estimated - centroid
print("means == 0?", np.mean(imu_data_mag_centrate[:,0]), np.mean(imu_data_mag_centrate[:,1]), np.mean(imu_data_mag_centrate[:,2]))
distances=np.sqrt(np.sum(imu_data_mag_centrate*imu_data_mag_centrate, axis=1))

print("distances.shape: ", distances.shape)
dist_media = np.mean(distances)
print("distances.mean: ", dist_media)
dist_varianza=np.var(distances)
print("distances.var: ", dist_varianza)

plt.hist(distances, bins='auto', alpha=0.7, density=True)

xmin, xmax = plt.xlim()
print("distances: xmin, xmax: ", xmin, xmax)
xlen = 100
# Genera un array di valori x per la curva della distribuzione normale
# Calcola i valori della funzione di densità di probabilità della distribuzione normale
x = np.linspace(xmin, xmax, xlen)
pdf = stats.norm.pdf(x, loc=dist_media, scale=np.sqrt(dist_varianza))
# Traccia la curva della distribuzione normale
plt.plot(x, pdf, 'r', label='Distribuzione Normale')

# Aggiunta delle linee verticali per i parametri
plt.axvline(dist_media, color='r', linestyle='dashed', linewidth=2, label='Media')
plt.axvline(dist_varianza, color='g', linestyle='dashed', linewidth=2, label='Varianza')
plt.legend()

# Aggiunta del testo per i parametri
testo = f"distance\nMedia: {dist_media:.2f}\nVarianza: {dist_varianza:.2f}\n"
plt.text(0.98, 0.80, testo, ha='right', va='top', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

# Mostra il grafico
plt.show()


mag_plot_data(imu_data_mag_estimated, title="imu_data_mag_estimated")

