import csv
import numpy as np
from scipy import stats
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

################################################
## Find eigenvector of min eigenvalue
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
   

class ImuEllipsoidEstimator:
   def __init__(self, dati):
     # Data
     self.dati = dati
     self.variance = self.approx_ellipsoid_std() ** 2
     # Index Matrix
     self.IM = np.array([[1,1],[1,2],[1,3],[2,2],[2,3],[3,3],[1,0],[2,0],[3,0],[0,0]])
     self.n = 3
     self.psi_als = self.make_psi_als()
     self.als = self.estimate_all()
     self.model = self.mag_estimate()
     estimated_data = self.mag_apply_estimator()

     # scale factors correction ...
     distances=np.sqrt(np.sum(estimated_data*estimated_data, axis=1))
     dist_media = np.mean(distances)
     self.model['scale_factors'] = self.model['scale_factors']/dist_media
     estimated_data = self.mag_apply_estimator()

     # Calculate mean and var of radius (expected mean close to 1)
     distances=np.sqrt(np.sum(estimated_data*estimated_data, axis=1))
     dist_media = np.mean(distances)
     dist_variance = np.var(distances)
     print("radius mean: ", dist_media)
     print("radius variance: ", dist_variance)
     assert np.isclose(dist_media, 1.0), "Il raggio medio deve essere circa uno"

     self.estimated_data = estimated_data

   ###############################
   ## Estimate ellipsoid std
   ## (very approxymated)
   ###############################
   def approx_ellipsoid_std(self): 
     ## Calcolo approssimativamente la deviazione standard del raggio
     centroid = np.array([np.mean(self.dati[:,0]), np.mean(self.dati[:,1]),np.mean(self.dati[:,2])])
     dati_centered = self.dati - centroid
     distances=np.sqrt(np.sum(dati_centered*dati_centered, axis=1))
     dist_media = np.mean(distances)
     dist_varianza=np.var(distances)
     radius_std = dist_varianza ** (1/3)
     radius_std -= radius_std*7.9/100
     return radius_std

   ###############################
   ## Form Tensor T 5xnxm
   ###############################
   def T(self, k, i, l):
     X = self.dati
     variance = self.variance
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

     switch_cases = {
       1: t0,
       2: t1,
       3: t2,
       4: t3,
       5: t4
     }
     return switch_cases.get(k+1, terror)(X[l,i], variance)

   #####################################
   ## Form Tensor R 
   #####################################
   def R(self, p,q,i):
     M = self.IM
     def is_equal(v1, v2):
       if v1 == v2:
         return 1
       else:
         return 0
     return is_equal(M[p,0], i) + is_equal(M[p,1], i)  + is_equal(M[q,0], i)  + is_equal(M[q,1], i) 

   #####################################
   ## Compute ni_als
   #####################################
   def ni_als(self, p,q):
      X = self.dati
      m = X.shape[0]
      n = X.shape[1]
      s = 0
      for l in range(m):
        f = 1
        for i in range(n):
          f = f*self.T(self.R(p,q,i+1),i,l)
        s = s + f
      return s
   
   #####################################
   ## form matrix psi_als
   #####################################
   def make_psi_als(self):
     X = self.dati
     M = self.IM

     # Define index off-diagonal
     def index_off_diagonal():
       # only for n=3
       return np.array([1,2,4])
   
     n = X.shape[1]
     n1 = M.shape[0]
     D = index_off_diagonal()
     psi_als = np.zeros((n1, n1)) 
     for p in range(n1):
       for q in range(p,n1):
         if p in D and q in D :
           psi_als[p,q] = 4*self.ni_als(p,q)
         else:
           if p not in D and q not in D:
             psi_als[p,q] = self.ni_als(p,q)
           else:
             psi_als[p,q] = 2*self.ni_als(p,q)
         psi_als[q,p] = psi_als[p,q]
     return psi_als
   
   ################################################
   ## Find eigenvector of min eigenvalue
   ################################################
   def calc_eigenvector_psi_als(self):
       psi_als = self.psi_als
       eigen_als = eigen(psi_als)
       eigen_values = eigen_als[0]
       eigen_vectors = eigen_als[1]
       idx = np.argmin(np.abs(eigen_values))
       beta = np.matrix(eigen_vectors[:, idx])
       result = beta / np.linalg.norm(beta, ord=2)
       return result
   
   ################################################
   ## Estimate A, b, c, d
   ################################################
   def estimate_A(self, beta):
       print("beta: ", beta, "n ", self.n)
       expanded_beta = beta[:(self.n*(self.n+1)//2)]
       print("expandend beta: ", expanded_beta)
       expanded_matrix = np.zeros((self.n, self.n))
   
       idx = 0
       for i in range(self.n):
           for j in range(i, self.n):
               expanded_matrix[i, j] = expanded_beta[0,idx]
               expanded_matrix[j, i] = expanded_beta[0,idx]
               idx += 1
   
       return expanded_matrix
   
   def estimate_b(self, beta):
       start_index = int(self.n*(self.n+1)/2)
       end_index = beta.shape[1] - 1
       result=beta[0,start_index:end_index]
       return np.transpose(result)
   
   def estimate_d(self, beta):
       return beta[0,-1]
   
   def estimate_c(self, A, b):
       c = -0.5 * np.linalg.solve(A, b)
       return c
   
   def estimate_Ae(self, A, c, d):
       denominator = np.dot(np.dot(np.transpose(c), A), c) - d
       scaling_factor = 1.0 / denominator
       print("scaling_factor: ", scaling_factor)
       Ae = scaling_factor.item() * A
       return Ae
   
   def to_definite_positive_factor(self, A):
       eigenvalues = eigen(A)[0]
       result = 0
       if np.prod(eigenvalues >= 0) == 1:
           result = 1
       elif np.prod(eigenvalues < 0) == 1:
           result = -1
   
       return result
   
   def estimate_all(self):
       beta = self.calc_eigenvector_psi_als()
       A = self.estimate_A(beta)
       # forza la matrice ad essere definita positiva se è negativa
       factor = self.to_definite_positive_factor(A)
       print("factor: ", factor)
       assert factor != 0, "La matrice non può essere indefinita"
       A = factor * A
       
       b = factor * self.estimate_b(beta)
       d = factor * self.estimate_d(beta)
       c = self.estimate_c(A, b)
       Ae = self.estimate_Ae(A, c, d)
       
       als = {}
       als['A'] = A
       als['b'] = b
       als['d'] = d
       als['c'] = c
       als['Ae'] = Ae
       
       return als
   
   def mag_estimate(self):
       als = self.als
       mag_data = self.dati
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
       ## START FIXME
       #f = 1 / np.sqrt(mag_model_magnetic_norm) # real ..
       f = 1.09 / np.sqrt(mag_model_magnetic_norm) # but works ..
       mag_model_scale_factors = np.zeros((3, 3))
       np.fill_diagonal(mag_model_scale_factors, f)
       ## provo a calcolare i factor in altro modo
       eigen_Q = eigen(mag_model_Q)
       factors=np.sqrt(1/eigen_Q[0])
       factors=factors/np.linalg.norm(factors, ord=2)
       print("factors: ", factors)
       mag_model_scale_factors = np.diag(factors)
       ##
       ## END FIXME
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
   
   def mag_apply_estimator(self):
       mag_model = self.model
       mag_data = mag_model['data_source']
       mag_model_inv_A = mag_model['invA']
       mag_model_b = mag_model['b']
       result = np.zeros(mag_data.shape)
       
       # Applica il modello ai dati
       for i in range(mag_data.shape[0]):
           x = np.transpose(mag_data[i, :] - np.transpose(mag_model_b))
           result[i, :] = np.transpose(np.dot(mag_model_inv_A, np.dot(mag_model['scale_factors'], x)))
       return result


def estimate_mag_acc(file_csv):
  #################################################################################
  ######### Fase 4: Calcolo offset e matrici di correzione
  #########         dal secondo set di dati utilizzando media, varianza e std
  #################################################################################
  dati = pd.read_csv(file_csv)

  ## Magnetometer
  dati_mag = np.array(dati[['MX', 'MY', 'MZ']].loc[dati['MV'] == 1])

  estimator_mag = ImuEllipsoidEstimator(dati_mag)

  print("Magnetometer Results:")
  print("Offset: ", estimator_mag.model['offset'])
  print("Matrix: ", estimator_mag.model['invA'])
  print("Scale Factors: ", estimator_mag.model['scale_factors'])
  print("Scaled Matrix: ", np.dot(estimator_mag.model['invA'], estimator_mag.model['scale_factors']))
  print("eigen(Scaled Matrix)", eigen(np.dot(estimator_mag.model['invA'], estimator_mag.model['scale_factors'])))

  ## Accelerometer
  dati_acc = np.array(dati[['AX', 'AY', 'AZ']])

  estimator_acc = ImuEllipsoidEstimator(dati_acc)

  print("Accelerometer Results:")
  print("Offset: ", estimator_acc.model['offset'])
  print("Matrix: ", estimator_acc.model['invA'])
  print("Scale Factors: ", estimator_acc.model['scale_factors'])
  print("Scaled Matrix: ", np.dot(estimator_acc.model['invA'], estimator_acc.model['scale_factors']))
  print("eigen(Scaled Matrix)", eigen(np.dot(estimator_acc.model['invA'], estimator_acc.model['scale_factors'])))

  return estimator_mag, estimator_acc

def estimate_gyro(file_csv):
  dati = pd.read_csv(file_csv)
  dati_gyro = np.array(dati[['GX', 'GY', 'GZ']])
  centroid = np.array([np.mean(dati_gyro[:,0]), np.mean(dati_gyro[:,1]),np.mean(dati_gyro[:,2])])
  gyro_variances = np.array([np.var(dati_gyro[:,0]), np.var(dati_gyro[:,1]),np.var(dati_gyro[:,2])])
  imu_data_gyro_centered = dati_gyro - centroid
  gyro_means = np.array([np.mean(imu_data_gyro_centered[:,0]), np.mean(imu_data_gyro_centered[:,1]),np.mean(imu_data_gyro_centered[:,2])])
  print("gyro centered means: ", gyro_means)
  assert np.allclose(gyro_means, 0.0), "Means must be close to zero ..."

  return np.identity(3, dtype=float), centroid, gyro_variances,imu_data_gyro_centered

def imu_ellipsoid_estimator_example():
   def prepare_data():
     ##################################################
     ### Data Preparation
     ##################################################
     df = pd.read_csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/02-imu-raw-data.csv')
     return df   

   def estimate_ellipsoid(imu_data_mag):
     centroid = np.array([np.mean(imu_data_mag[:,0]), np.mean(imu_data_mag[:,1]),np.mean(imu_data_mag[:,2])])
     imu_data_mag_centered = imu_data_mag - centroid
     print("post: means == 0?", np.mean(imu_data_mag_centered[:,0]), np.mean(imu_data_mag_centered[:,1]), np.mean(imu_data_mag_centered[:,2]))

     distances=np.sqrt(np.sum(imu_data_mag_centered*imu_data_mag_centered, axis=1))

     print("raw distances.shape: ", distances.shape)
     dist_media = np.mean(distances)
     print("raw distances.mean: ", dist_media)
     dist_varianza=np.var(distances)
     print("raw distances.var: ", dist_varianza)
     raw_std = dist_varianza ** (1/3)
     raw_std -= raw_std * 7.9 / 100.0
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
     plt.legend()

     # Aggiunta del testo per i parametri
     testo = f"distance\nMedia: {dist_media:.2f}\nVarianza: {dist_varianza:.2f}\n"
     plt.text(0.98, 0.80, testo, ha='right', va='top', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

     # Mostra il grafico
     plt.show()

     ##################################################
     ### Estimation
     ##################################################
     #estimator = ImuEllipsoidEstimator(np.matrix(imu_data_mag), (6.054304 ** 2))
     #estimator = ImuEllipsoidEstimator(np.matrix(imu_data_mag), (11.054304 ** 2))
     #estimator = ImuEllipsoidEstimator(np.matrix(imu_data_mag), (11.025304 ** 2))
     estimator = ImuEllipsoidEstimator(np.matrix(imu_data_mag), (raw_std ** 2))
     imu_data_mag_estimated = estimator.estimated_data
     distances=np.sqrt(np.sum(imu_data_mag_estimated*imu_data_mag_estimated, axis=1))
   
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

     ##################################################
     ### Plot
     ##################################################
     def mag_plot_data(mag_data, sphere_radius=-1, title=""):
       # Genera le coordinate (x, y, z) per la sfera di raggio 1
       theta = np.linspace(0, 2 * np.pi, 100)
       phi = np.linspace(0, np.pi, 100)
       x1 = np.outer(np.cos(theta), np.sin(phi))
       y1 = np.outer(np.sin(theta), np.sin(phi))
       z1 = np.outer(np.ones(np.size(theta)), np.cos(phi))

       fig = plt.figure(figsize=(8,6), dpi=300)
       ax = fig.add_subplot(111, projection='3d')

       # Disegna la sfera di raggio 1
       ax.plot_surface(x1, y1, z1, color='b', alpha=0.8)

       ax.scatter(mag_data[:, 0], mag_data[:, 1], mag_data[:, 2], c='g', cmap=None, alpha=0.4)
       ax.set_xlabel('X')
       ax.set_ylabel('Y')
       ax.set_zlabel('Z')

       ax.set_xlim([-1.005, 1.005])
       ax.set_ylim([-1.005, 1.005])
       ax.set_zlim([-1.005, 1.005])

       ax.set_title(title)
       # Imposta la scala degli assi in modo da essere uniforme
       plt.show()

     mag_plot_data(imu_data_mag_estimated, title="imu_data_mag_estimated")
   
     print("Results:")
     print("Offset: ", estimator.model['offset'])
     print("Matrix: ", estimator.model['invA'])
     print("Scale Factors: ", estimator.model['scale_factors'])
     print("Scaled Matrix: ", np.dot(estimator.model['invA'], estimator.model['scale_factors']))
     print("eigen(Scaled Matrix)", eigen(np.dot(estimator.model['invA'], estimator.model['scale_factors'])))
     return estimator

   df = prepare_data()
   imu_data_mag = np.array(df[['MX', 'MY', 'MZ']].loc[df['MV'] == 1])
   mag_estimator = estimate_ellipsoid(imu_data_mag)
   imu_data_acc = np.array(df[['AX', 'AY', 'AZ']])
   acc_estimator = estimate_ellipsoid(imu_data_acc)

