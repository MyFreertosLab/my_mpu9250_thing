/*
 * my_mp9250_task.c
 *
 *  Created on: 29 gen 2021
 *      Author: andrea
 */
#include "my_mpu9250_sensor_reader_task.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <esp_system.h>
#include <sys/socket.h>
#include <errno.h>
#include <netdb.h>
#include <arpa/inet.h>
#include "esp_netif.h"
#include "esp_log.h"
#include <nvs_flash.h>
#include <nvs.h>
#include <driver/spi_master.h>
#include <driver/gpio.h>
#include <mpu9250_accel.h>
#include <mpu9250_gyro.h>
#include <mpu9250_mag.h>
#include <mpu9250_baro.h>
#include <mpu9250_queue_definitions.h>

static const char *TAG = "my_mpu9250_sensor_reader";


static message_type_enum_t message_type = MESSAGE_TYPE_ACC_MODEL;

static nvs_handle_t my_nvs_handle;

// Initialize mpu9250_handle
static mpu9250_init_t mpu9250;
mpu9250_handle_t mpu9250_handle = &mpu9250;

esp_err_t my_mpu9250_init_flash(){
    esp_err_t err = nvs_flash_init();
    if (err == ESP_ERR_NVS_NO_FREE_PAGES || err == ESP_ERR_NVS_NEW_VERSION_FOUND) {
        // NVS partition was truncated and needs to be erased
        // Retry nvs_flash_init
        ESP_ERROR_CHECK(nvs_flash_erase());
        err = nvs_flash_init();
    }
	return err;
}

static void my_mpu9250_task_init(mpu9250_handle_t mpu9250) {
	ESP_ERROR_CHECK(my_mpu9250_init_flash());
    ESP_ERROR_CHECK(nvs_open("CAL_DATA", NVS_READWRITE, &my_nvs_handle));
    mpu9250_handle->data.nvs_cal_data = my_nvs_handle;
}

void my_mpu9250_temperature_read_data_cycle(mpu9250_handle_t mpu9250_handle) {
	const TickType_t xMaxBlockTime = pdMS_TO_TICKS( 500 );
	uint32_t counter = 0;

	while (true) {
		counter++;
		if( ulTaskNotifyTake( pdTRUE,xMaxBlockTime ) == 1) {
			ESP_ERROR_CHECK(mpu9250_load_raw_data(mpu9250_handle));
			// TODO: probabilmente offset non è 0. Verificare
			if(counter%100 == 0) {
				float temp_deg = mpu9250_handle->data.raw_data.data_s_xyz.temp_data/333.87 + 21;

				printf("Temperature: [%2.3f][%d]\n", temp_deg, mpu9250_handle->data.raw_data.data_s_xyz.temp_data);
			}
	    } else {
	    	ESP_ERROR_CHECK(mpu9250_test_connection(mpu9250_handle));
			if(counter%100 == 0) {
		    	printf("SORRY!! Interrupt LOST!\n");
			}
	    }
	}
}

static esp_err_t mpu9250_discard_messages(mpu9250_handle_t mpu9250_handle, uint16_t num_msgs) {
	printf("Discarding %d Samples ... \n", num_msgs);
	for(uint16_t i = 0; i < num_msgs; i++) {
		ulTaskNotifyTake( pdTRUE,pdMS_TO_TICKS( 500 ) );
	}
	return ESP_OK;
}

void my_mpu9250_read_data_cycle(mpu9250_handle_t mpu9250_handle) {
	const TickType_t xMaxBlockTime = pdMS_TO_TICKS( 500 );

	ESP_ERROR_CHECK(mpu9250_gyro_set_fsr(mpu9250_handle, INV_FSR_1000DPS));
	ESP_ERROR_CHECK(mpu9250_acc_set_fsr(mpu9250_handle, INV_FSR_8G));
    ESP_ERROR_CHECK(mpu9250_mag_set_mode2_with_precision(mpu9250_handle, INV_MAG_PRECISION_16_BITS));
	ESP_ERROR_CHECK(mpu9250_mag_set_continuous_reading(mpu9250_handle));
	ESP_ERROR_CHECK(mpu9250_discard_messages(mpu9250_handle, 10000));

	uint32_t counter = 0;

	sensor_message_t data_message_cal;
	uint16_t data_size_cal = sizeof(sensor_message_t);


    while(true) {
      counter++;
   	  counter %= 100;
      uint8_t activate_send = __atomic_load_n(&my_mpu9250_activate_send_message, __ATOMIC_SEQ_CST);
      if( ulTaskNotifyTake( pdTRUE,xMaxBlockTime ) == 1) {
        ESP_ERROR_CHECK(mpu9250_load_data(mpu9250_handle));
    	if(counter == 0) {
    	  ESP_LOGD(TAG, "Message Type: %d\n", message_type);
    	  ESP_LOGI(TAG, "Gyro .: [%3.5f][%3.5f][%3.5f]\n", mpu9250_handle->data.cal_data.data_s_xyz.gyro_data_x, mpu9250_handle->data.cal_data.data_s_xyz.gyro_data_y, mpu9250_handle->data.cal_data.data_s_xyz.gyro_data_z);
    	  ESP_LOGI(TAG, "Accel : [%3.5f][%3.5f][%3.5f]\n", mpu9250_handle->data.cal_data.data_s_xyz.accel_data_x, mpu9250_handle->data.cal_data.data_s_xyz.accel_data_y, mpu9250_handle->data.cal_data.data_s_xyz.accel_data_z);
    	  ESP_LOGI(TAG, "Mag ..: [%3.5f][%3.5f][%3.5f]\n", mpu9250_handle->data.cal_data.data_s_xyz.mag_data_x, mpu9250_handle->data.cal_data.data_s_xyz.mag_data_y, mpu9250_handle->data.cal_data.data_s_xyz.mag_data_z);

    	  if(false) {
			//while(data_message_cal.type != 0) {}; // Attendo che un eventuale messaggio precedente sia stato consegnato.
            memset(&data_message_cal, '\0', sizeof(sensor_message_t)); // elimino eventuali valori residui nel buffer
			// preparo il messaggio
			data_message_cal.type = MESSAGE_TYPE_CONFIG;
			data_message_cal.data.mpu9250_config_data.accel_fsr = mpu9250_handle->data.accel.fsr;
			data_message_cal.data.mpu9250_config_data.gyro_fsr = mpu9250_handle->data.gyro.fsr;
			data_message_cal.data.mpu9250_config_data.mag_precision = mpu9250_handle->data.mag.precision;
			data_message_cal.data.mpu9250_config_data.frequenzy_hz = 500;

            // invio il messaggio
            if (xQueueSend(mpu9250_output_queue, &data_message_cal, pdMS_TO_TICKS(100)) != pdPASS) {
                ESP_LOGW(TAG, "Coda piena, impossibile inviare il messaggio");
                memset(&data_message_cal, '\0', sizeof(sensor_message_t)); // libero il buffer
            }
    	  }
        }
    	if(activate_send) {
		  //while(data_message_cal.type != 0) {}; // Attendo che un eventuale messaggio precedente sia stato consegnato.
          memset(&data_message_cal, '\0', sizeof(sensor_message_t)); // elimino eventuali valori residui nel buffer

          // copy data to send
 		  data_message_cal.type = MESSAGE_TYPE_DATA_ACC_GYRO;
          memcpy((char*)&data_message_cal.data.imu_data.accel, (char*)mpu9250_handle->data.cal_data.data_s_vector.accel, sizeof(mpu9250_handle->data.cal_data.data_s_vector.accel));
          memcpy((char*)&data_message_cal.data.imu_data.gyro, (char*)mpu9250_handle->data.cal_data.data_s_vector.gyro, sizeof(mpu9250_handle->data.cal_data.data_s_vector.gyro));
          data_message_cal.data.imu_data.temp = mpu9250_handle->data.cal_data.data_s_vector.temp;
          data_message_cal.timestamp = mpu9250_handle->data.timestamp;

          // Invia i dati
          if (xQueueSend(mpu9250_output_queue, &data_message_cal, pdMS_TO_TICKS(100)) != pdPASS) {
            ESP_LOGW(TAG, "Coda piena, impossibile inviare i dati");
			memset(&data_message_cal, '\0', sizeof(sensor_message_t)); // libero il buffer

          }

		  if(mpu9250_handle->data.mag.drdy) {
		    //while(data_message_cal.type != 0) {}; // Attendo che un eventuale messaggio precedente sia stato consegnato.
            memset(&data_message_cal, '\0', sizeof(sensor_message_t)); // elimino eventuali valori residui nel buffer
 		    data_message_cal.type = MESSAGE_TYPE_DATA_MAG;
            memcpy((char*)&data_message_cal.data.magneto_data, (char*)mpu9250_handle->data.cal_data.data_s_vector.mag, sizeof(mpu9250_handle->data.cal_data.data_s_vector.mag));
            data_message_cal.timestamp = mpu9250_handle->data.timestamp;
            if (xQueueSend(mpu9250_output_queue, &data_message_cal, pdMS_TO_TICKS(100)) != pdPASS) {
              ESP_LOGW(TAG, "Coda piena, impossibile inviare i dati");
			  memset(&data_message_cal, '\0', sizeof(sensor_message_t)); // libero il buffer
            }
		  }
    	}
      }
    }
}

void my_mpu9250_sensor_reader_task(void *arg) {
    // Avvia MQTT client
//    mqtt_app_start();
    ESP_LOGI(TAG, "Task STARTED!");
	my_mpu9250_task_init(mpu9250_handle);
    ESP_LOGI(TAG, "Task INITIALIZED!");

	// Init MPU9250
	ESP_ERROR_CHECK(mpu9250_init(mpu9250_handle));
    ESP_LOGI(TAG, "MPU INITIALIZED!");

	my_mpu9250_read_data_cycle(mpu9250_handle);
}
