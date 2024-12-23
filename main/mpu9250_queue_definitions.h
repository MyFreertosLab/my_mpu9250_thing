/*
 * calibration_queue.h
 *
 *  Created on: 22 dic 2024
 *      Author: andrea
 */
#ifndef MPU9250_QUEUE_DEFINITIONS_H
#define MPU9250_QUEUE_DEFINITIONS_H

#include <freertos/FreeRTOS.h>
#include "freertos/queue.h"
#include <stddef.h>
#include <stdint.h>
#include <mpu9250.h>

typedef enum  {
	MESSAGE_TYPE_DATA_ACC_GYRO = 1,
	MESSAGE_TYPE_DATA_MAG,
	MESSAGE_TYPE_CONFIG,
	MESSAGE_TYPE_COMMAND_START,
	MESSAGE_TYPE_COMMAND_STOP,
	MESSAGE_TYPE_COMMAND_DISABLE_KALMAN,
	MESSAGE_TYPE_COMMAND_ENABLE_KALMAN,
	MESSAGE_TYPE_COMMAND_DISABLE_KALMAN_ACC,
	MESSAGE_TYPE_COMMAND_ENABLE_KALMAN_ACC,
	MESSAGE_TYPE_COMMAND_DISABLE_KALMAN_GYRO,
	MESSAGE_TYPE_COMMAND_ENABLE_KALMAN_GYRO,
	MESSAGE_TYPE_COMMAND_DISABLE_KALMAN_MAG,
	MESSAGE_TYPE_COMMAND_ENABLE_KALMAN_MAG,
	MESSAGE_TYPE_MAG_MODEL,
	MESSAGE_TYPE_ACC_MODEL,
	MESSAGE_TYPE_GYRO_MODEL
} message_type_enum_t;

// Struttura per dati dal sensore
typedef struct {
    message_type_enum_t type; // Tipo di messaggio
    uint32_t timestamp;  // Timestamp
    union {
        struct {
            float accel[3];   // Accelerometro (X, Y, Z)
            float gyro[3];    // Giroscopio (X, Y, Z)
            float temp;       // temperature
        } imu_data;          // Dati IMU
        struct {
            float magneto[3]; // Magnetometro (X, Y, Z)
        } magneto_data;      // Dati Magnetometro
        struct {
        	uint8_t accel_fsr;
        	uint8_t gyro_fsr;
        	uint8_t mag_precision;
        	uint8_t mag_drdy;
        	uint16_t frequenzy_hz;
        } mpu9250_config_data;
		struct {
			uint8_t blob[60]; // 3x3 matrix + 3x1 bias + 3x1 variances = 15 floats = 15 * sizeof(float) = 15x4 = 60 bytes
		} sensor_model;
    } data;                  // Dati unificati
} sensor_message_t;

// Dichiarazione delle code
extern QueueHandle_t mpu9250_output_queue;
extern QueueHandle_t mpu9250_input_queue;
extern uint8_t my_mpu9250_activate_send_message;

#endif // MPU9250_QUEUE_DEFINITIONS_H
