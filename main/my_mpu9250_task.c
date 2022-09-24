/*
 * my_mp9250_task.c
 *
 *  Created on: 29 gen 2021
 *      Author: andrea
 */
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

#include <my_mpu9250_task.h>
#include <mpu9250_accel.h>
#include <mpu9250_gyro.h>
#include <mpu9250_mag.h>
#include <mpu9250_baro.h>

static const char *MY_MPU9250_TAG = "my_mpu9250";
#define HOST_IP "192.168.1.64"
#define HOST_PORT 65432


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


void my_mpu9250_task_init(mpu9250_handle_t mpu9250) {
	ESP_ERROR_CHECK(my_mpu9250_init_flash());
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

    int sock =  -1;
    while(true) {
    	while(true) {
    	    char host_ip[] = HOST_IP;
    	    int addr_family = AF_INET;
    	    int ip_protocol = IPPROTO_IP;
    	    struct sockaddr_in dest_addr;
    	    inet_pton(AF_INET, host_ip, &dest_addr.sin_addr);
    	    dest_addr.sin_family = AF_INET;
    	    dest_addr.sin_port = htons(HOST_PORT);

    		sock = socket(addr_family, SOCK_STREAM, ip_protocol);
    		if (sock < 0) {
    			ESP_LOGE(MY_MPU9250_TAG, "Unable to create socket: errno %d",
    					errno);
    			break;
    		}
    	    ESP_LOGI(MY_MPU9250_TAG, "Socket created, connecting to %s:%d", host_ip, HOST_PORT);

    		int err = connect(sock, (struct sockaddr*) &dest_addr,
    				sizeof(dest_addr));
    		if (err != 0) {
    			ESP_LOGE(MY_MPU9250_TAG, "Socket unable to connect: errno %d",
    					errno);
    			break;
    		}
    	    ESP_LOGI(MY_MPU9250_TAG, "Successfully connected");

    		uint32_t counter = 0;
    		while (true) {
    			counter++;
    			if( ulTaskNotifyTake( pdTRUE,xMaxBlockTime ) == 1) {
    				counter %= 100;
    				ESP_ERROR_CHECK(mpu9250_load_data(mpu9250_handle));
    				if(counter == 0) {
    					printf("Gyro .: [%d][%d][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.gyro_data_x, mpu9250_handle->data.raw_data.data_s_xyz.gyro_data_y, mpu9250_handle->data.raw_data.data_s_xyz.gyro_data_z);
    					printf("Accel : [%d][%d][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.accel_data_x, mpu9250_handle->data.raw_data.data_s_xyz.accel_data_y, mpu9250_handle->data.raw_data.data_s_xyz.accel_data_z);
    					printf("Mag ..: [%d][%d][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.mag_data_x, mpu9250_handle->data.raw_data.data_s_xyz.mag_data_y, mpu9250_handle->data.raw_data.data_s_xyz.mag_data_z);
    				}
    				int err = send(sock, &mpu9250_handle->data.raw_data.data_s_xyz,
    						sizeof(mpu9250_raw_data_t), 0);
    				if (err < 0) {
    					ESP_LOGE(MY_MPU9250_TAG,
    							"Error occurred during sending: errno %d", errno);
    					break;
    				}
    		    } else {
    		    	ESP_ERROR_CHECK(mpu9250_test_connection(mpu9250_handle));
    				if(counter == 0) {
    			    	printf("SORRY!! Interrupt LOST!\n");
    				}
    		    }
    		}

    		break;
    	}

        if (sock != -1) {
            ESP_LOGE(MY_MPU9250_TAG, "Shutting down socket and restarting...");
            shutdown(sock, 0);
            close(sock);
        }
        vTaskDelay(pdMS_TO_TICKS(5000));
    }
}

void my_mpu9250_task(void *arg) {
	const TickType_t xMaxBlockTime = pdMS_TO_TICKS( 500 );

	// MPU9250 Handle
	mpu9250_init_t mpu9250;
	mpu9250_handle_t mpu9250_handle = &mpu9250;
	my_mpu9250_task_init(mpu9250_handle);

	// Init MPU9250
	ESP_ERROR_CHECK(mpu9250_init(mpu9250_handle));

	// load circular buffer
	for(uint8_t i = 0; i < CIRCULAR_BUFFER_SIZE; i++) {
		if( ulTaskNotifyTake( pdTRUE,xMaxBlockTime ) == 1) {
			ESP_ERROR_CHECK(mpu9250_load_raw_data(mpu9250_handle));
		}
	}
	my_mpu9250_read_data_cycle(mpu9250_handle);
}
