import csv
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import pandas as pd

file_csv = 'examples/01-imu-cal-data.csv'
dati = pd.read_csv(file_csv)
dati_mag = np.array(dati[['MX', 'MY', 'MZ']].loc[dati['MV'] == 1])
dati_gyro = np.array(dati[['GX', 'GY', 'GZ']])
dati_acc = np.array(dati[['AX', 'AY', 'AZ']])

means_gyro = np.array([np.mean(dati_gyro[:,0]), np.mean(dati_gyro[:,1]),np.mean(dati_gyro[:,2])])
variance_gyro = np.array([np.var(dati_gyro[:,0]), np.var(dati_gyro[:,1]),np.var(dati_gyro[:,2])])
print("means_gyro: ", means_gyro)
print("variance_gyro: ", variance_gyro)

means_acc = np.array([np.mean(dati_acc[:,0]), np.mean(dati_acc[:,1]),np.mean(dati_acc[:,2])])
variance_acc = np.array([np.var(dati_acc[:,0]), np.var(dati_acc[:,1]),np.var(dati_acc[:,2])])
print("means_acc: ", means_acc)
print("variance_acc: ", variance_acc)

means_mag = np.array([np.mean(dati_mag[:,0]), np.mean(dati_mag[:,1]),np.mean(dati_mag[:,2])])
variance_mag = np.array([np.var(dati_mag[:,0]), np.var(dati_mag[:,1]),np.var(dati_mag[:,2])])
print("means_mag: ", means_mag)
print("variance_mag: ", variance_mag)

