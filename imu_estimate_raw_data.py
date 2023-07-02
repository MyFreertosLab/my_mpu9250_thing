import numpy as np
import pandas as pd
import struct
import time
import imu_ellipsoid_estimator as iee
from imu_ellipsoid_estimator import estimate_mag_acc,estimate_gyro

################################################
## Gyroscope
################################################
file_csv_static = "examples/dati-reali/20230702-1945-01-imu-raw-data.csv"
gyro_model_matrix, gyro_model_bias, gyro_model_bias_variance,_ = estimate_gyro(file_csv_static)
print("Gyro Matrix: ", gyro_model_matrix)
print("Gyro bias: ", gyro_model_bias)
print("Gyro bias variance: ", gyro_model_bias_variance)

################################################
## Magnetometer
################################################
file_csv = "examples/dati-reali/20230701-1053-02-imu-raw-data.csv"
estimator_mag, estimator_acc = estimate_mag_acc(file_csv)

mag_model_matrix = np.dot(estimator_mag.model['invA'], estimator_mag.model['scale_factors'])
mag_model_bias = estimator_mag.model['b']
print("Mag Matrix: ", mag_model_matrix)
print("Mag bias: ", mag_model_bias)

acc_model_matrix = np.dot(estimator_acc.model['invA'], estimator_acc.model['scale_factors'])
acc_model_bias = estimator_acc.model['b']
print("Acc Matrix: ", acc_model_matrix)
print("Acc bias: ", acc_model_bias)

