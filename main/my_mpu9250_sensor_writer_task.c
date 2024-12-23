/*
 * my_mpu9250_sensor_writer.c
 *
 *  Created on: 22 dic 2024
 *      Author: andrea
 */
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <string.h>
#include <esp_log.h>
#include <mpu9250_accel.h>
#include <mpu9250_gyro.h>
#include <mpu9250_mag.h>
#include <mpu9250_baro.h>
#include <mpu9250_queue_definitions.h>
#include <my_mpu9250_sensor_writer_task.h>

static const char *TAG = "my_mpu9250_sensor_writer";

// global variable for enable/disable sending of sensor data
uint8_t my_mpu9250_activate_send_message = 1;

void my_mpu9250_sensor_writer_task(void *arg)
{
    sensor_message_t sensor_msg;

    while (true) {
        // Attendi dati dalla coda del sensore ed aggiorna il sensore
        if (xQueueReceive(mpu9250_input_queue, &sensor_msg, portMAX_DELAY) == pdTRUE) {
            ESP_LOGI(TAG, "Ricezione dati: %d", sensor_msg.type);
            switch(sensor_msg.type) {
              case MESSAGE_TYPE_COMMAND_START:
                __atomic_store_n(&my_mpu9250_activate_send_message, 1, __ATOMIC_SEQ_CST);
            	break;
              case MESSAGE_TYPE_COMMAND_STOP:
                __atomic_store_n(&my_mpu9250_activate_send_message, 0, __ATOMIC_SEQ_CST);
             	break;
              case MESSAGE_TYPE_MAG_MODEL:
                  mpu9250_mag_save_calibration_params(mpu9250_handle, (char*)&sensor_msg.data.sensor_model, sizeof(sensor_msg.data.sensor_model));
                  mpu9250_mag_load_calibration_params(mpu9250_handle);
            	  break;
              case MESSAGE_TYPE_ACC_MODEL:
                  mpu9250_acc_save_calibration_params(mpu9250_handle, (char*)&sensor_msg.data.sensor_model, sizeof(sensor_msg.data.sensor_model));
                  mpu9250_acc_load_calibration_params(mpu9250_handle);
            	  break;
              case MESSAGE_TYPE_GYRO_MODEL:
                  mpu9250_gyro_save_calibration_params(mpu9250_handle, (char*)&sensor_msg.data.sensor_model, sizeof(sensor_msg.data.sensor_model));
                  mpu9250_gyro_load_calibration_params(mpu9250_handle);
            	  break;
              default:
                  ESP_LOGW(TAG, "Topic non riconosciuto: %d", sensor_msg.type);
            	  break;
            }
            memset(&sensor_msg, '\0', sizeof(sensor_message_t)); // libero il buffer
        }
    }
}

