import numpy as np
import pandas as pd
import struct
import time
import imu_ellipsoid_estimator as iee
from imu_ellipsoid_estimator import estimate_mag_acc,estimate_static_gyro,estimate_static_acc, estimate_static_mag,ImuOffsetKalmanEstimator

################################################
## Gyroscope
################################################
file_csv_static = "examples/01-imu-raw-data.csv"
gyro_kalman_R, gyro_model_bias, gyro_model_variance,gyro_model_covariance,gyro_kalman_state, _ = estimate_static_gyro(file_csv_static)
print("Gyro bias: ", gyro_model_bias)
print("Gyro variance: ", gyro_model_variance)
print("Gyro covariance: ", gyro_model_covariance)
print("Gyro kalman X: ", gyro_kalman_state)
print("Gyro kalman R: ", gyro_kalman_R)

################################################
## Magnetometer Static
################################################
mag_kalman_R, mag_model_bias, mag_model_variance,mag_model_covariance,mag_kalman_state, _ = estimate_static_mag(file_csv_static)
print("Mag means: ", mag_model_bias)
print("Mag variance: ", mag_model_variance)
print("Mag covariance: ", mag_model_covariance)
print("Mag kalman X: ", mag_kalman_state)
print("Mag kalman R: ", mag_kalman_R)

################################################
## Accelerometer Static
################################################
acc_kalman_R, acc_model_bias, acc_model_variance,acc_model_covariance,acc_kalman_state, _ = estimate_static_acc(file_csv_static)
print("Acc means: ", acc_model_bias)
print("Acc variance: ", acc_model_variance)
print("Acc covariance: ", acc_model_covariance)
print("Acc kalman X: ", acc_kalman_state)
print("Acc kalman R: ", acc_kalman_R)

################################################
## Magnetometer Ellipsoid
################################################

file_csv = "examples/dati-reali/20230701-1053-02-imu-raw-data.csv"
estimator_mag, estimator_acc = estimate_mag_acc(file_csv)

mag_model_matrix = np.dot(estimator_mag.model['invA'], estimator_mag.model['scale_factors'])
mag_model_bias = estimator_mag.model['b']
print("Mag Matrix: ", mag_model_matrix)
print("Mag bias: ", mag_model_bias)

################################################
## Accelerometer Ellipsoid
################################################
acc_model_matrix = np.dot(estimator_acc.model['invA'], estimator_acc.model['scale_factors'])
acc_model_bias = estimator_acc.model['b']
print("Acc Matrix: ", acc_model_matrix)
print("Acc bias: ", acc_model_bias)
