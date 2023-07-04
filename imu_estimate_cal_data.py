import csv
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import pandas as pd
from imu_ellipsoid_estimator import estimate_static_gyro,estimate_static_acc, estimate_static_mag

file_csv_static = 'examples/01-imu-cal-data.csv'

################################################
## Gyroscope static
################################################
gyro_kalman_R, gyro_model_bias, gyro_model_variance, gyro_model_covariance,gyro_kalman_state,_ = estimate_static_gyro(file_csv_static)
print("Gyro bias: ", gyro_model_bias)
print("Gyro variance: ", gyro_model_variance)
print("Gyro kalman R: ", gyro_kalman_R)
print("Gyro kalman X: ", gyro_kalman_state)
print("Gyro covariance: ", gyro_model_covariance)

################################################
## Magnetometer static
################################################
mag_kalman_R, _, mag_model_variance, mag_model_covariance,mag_kalman_state,_ = estimate_static_mag(file_csv_static)
print("Mag variance: ", mag_model_variance)
print("Mag kalman R: ", mag_kalman_R)
print("Mag kalman X: ", mag_kalman_state)
print("Mag covariance: ", mag_model_covariance)

################################################
## Accelerometer static
################################################
acc_kalman_R, _, acc_model_variance, acc_model_covariance,acc_kalman_state,_ = estimate_static_acc(file_csv_static)
print("Acc variance: ", acc_model_variance)
print("Acc kalman R: ", acc_kalman_R)
print("Acc kalman X: ", acc_kalman_state)
print("Acc covariance: ", acc_model_covariance)

