import numpy as np
import pandas as pd
import paho.mqtt.client as mqtt
import struct
import time
import imu_ellipsoid_estimator as iee
from imu_ellipsoid_estimator import estimate_mag_acc,estimate_gyro

# Configurazione connessione MQTT
mqtt_broker = "192.168.1.60"
mqtt_port = 1883
mqtt_topic_mag = "/imu/calibration/mag/model"
mqtt_topic_acc = "/imu/calibration/acc/model"
mqtt_topic_gyro = "/imu/calibration/gyro/model"

# Creazione del client MQTT
client = mqtt.Client()
client.keepalive = 60

# Connessione al broker MQTT
client.connect(mqtt_broker, mqtt_port)

client.loop_start()

################################################
## Gyroscope
################################################
file_csv_static = "examples/01-imu-raw-data.csv"
gyro_model_matrix, gyro_model_bias, gyro_model_variance, gyro_model_covariance,kalman_state,_ = estimate_gyro(file_csv_static)
print("Gyro Matrix: ", gyro_model_matrix)
print("Gyro bias: ", gyro_model_bias)
print("Gyro variance: ", gyro_model_variance)
print("Gyro covariance: ", gyro_model_covariance)

################################################
## Magnetometer
################################################
file_csv = "examples/dati-reali/20230701-1053-02-imu-raw-data.csv"
estimator_mag, estimator_acc = estimate_mag_acc(file_csv)

mag_model_matrix = np.dot(estimator_mag.model['invA'], estimator_mag.model['scale_factors'])
mag_model_bias = estimator_mag.model['b']
print("Mag Matrix: ", mag_model_matrix)
print("Mag bias: ", mag_model_bias)

################################################
## Send Data
################################################
# preparo payload da inviare
blob = b''

for row in gyro_model_matrix:
    for value in row:
        binary_value = struct.pack('f', value.item())  # Converte il float in formato binario
        blob += binary_value
for value in kalman_state:
        print("value bias: ", value.item())
        binary_value = struct.pack('f', value.item())  # Converte il float in formato binario
        blob += binary_value
  
print(blob)
print("len(blob)", len(blob))

# Pubblicazione del blob sul topic MQTT
client.publish(mqtt_topic_gyro, payload=blob)
time.sleep(10)

# preparo payload da inviare
blob = b''

for row in mag_model_matrix:
    for value in row:
        binary_value = struct.pack('f', value.item())  # Converte il float in formato binario
        blob += binary_value
for row in mag_model_bias:
    for value in row:
        print("value bias: ", value.item())
        binary_value = struct.pack('f', value.item())  # Converte il float in formato binario
        blob += binary_value
  
print(blob)
print("len(blob)", len(blob))

# Pubblicazione del blob sul topic MQTT
client.publish(mqtt_topic_mag, payload=blob)
time.sleep(10)

################################################
## Accelerometer
################################################
acc_model_matrix = np.dot(estimator_acc.model['invA'], estimator_acc.model['scale_factors'])
acc_model_bias = estimator_acc.model['b']
print("Acc Matrix: ", acc_model_matrix)
print("Acc bias: ", acc_model_bias)

# preparo payload da inviare
blob = b''

for row in acc_model_matrix:
    for value in row:
        binary_value = struct.pack('f', value.item())  # Converte il float in formato binario
        blob += binary_value
for row in acc_model_bias:
    for value in row:
        print("value bias: ", value.item())
        binary_value = struct.pack('f', value.item())  # Converte il float in formato binario
        blob += binary_value

print(blob)
print("len(blob)", len(blob))

# Pubblicazione del blob sul topic MQTT
client.publish(mqtt_topic_acc, payload=blob)

time.sleep(10)

# Disconnessione dal broker MQTT
client.disconnect()

