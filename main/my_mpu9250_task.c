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
#include <mqtt_client.h>

#include <my_mpu9250_task.h>
#include <mpu9250_accel.h>
#include <mpu9250_gyro.h>
#include <mpu9250_mag.h>
#include <mpu9250_baro.h>

static const char *TAG = "my_mpu9250";
#define HOST_IP "192.168.1.60"
#define HOST_PORT 65432
#define BUFF_LEN 128

#define MY_MPU9250_SENSORS_MSG_CONFIG_START "<CN>"
#define MY_MPU9250_SENSORS_MSG_CONFIG_END "</CN>"
#define MY_MPU9250_SENSORS_MSG_RAW_DATA_START "<RW>"
#define MY_MPU9250_SENSORS_MSG_RAW_DATA_END "</RW>"

typedef struct my_raw_data_s {
	mpu9250_raw_data_t data;
	uint8_t mag_data_drdy;
	uint32_t timestamp;
} mpu9250_raw_message_t;

typedef struct my_cal_data_s {
	mpu9250_calibrated_data_t data;
	uint8_t mag_data_drdy;
	uint32_t timestamp;
} mpu9250_cal_message_t;

typedef struct {
	uint8_t accel_fsr;
	uint8_t gyro_fsr;
	uint8_t mag_precision;
	uint8_t mag_drdy;
	uint16_t frequenzy_hz;
} mpu9250_config_data_t;


static esp_err_t mqtt_event_handler_cb(esp_mqtt_event_handle_t event)
{
    esp_mqtt_client_handle_t client = event->client;
    int msg_id;

    switch (event->event_id) {
        case MQTT_EVENT_CONNECTED:
            ESP_LOGI(TAG, "Connesso al broker MQTT");
            msg_id = esp_mqtt_client_subscribe(client, "/imu/calibration/commands", 0);
            ESP_LOGI(TAG, "Sottoscrizione al topic esempio, msg_id=%d", msg_id);
            break;

        case MQTT_EVENT_DISCONNECTED:
            ESP_LOGI(TAG, "Disconnesso dal broker MQTT");
            break;

        case MQTT_EVENT_DATA:
            ESP_LOGI(TAG, "Dati ricevuti sul topic %.*s : %.*s",
                     event->topic_len, event->topic, event->data_len, event->data);

            break;

        default:
            break;
    }
    return ESP_OK;
}

void mqtt_app_start(void)
{
    esp_mqtt_client_config_t mqtt_cfg = {
        .uri = "192.168.1.60", // Indirizzo del broker MQTT
        .event_handle = mqtt_event_handler_cb,
    };

    esp_mqtt_client_handle_t client = esp_mqtt_client_init(&mqtt_cfg);
    esp_mqtt_client_start(client);
}

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

static esp_err_t mpu9250_send_message(int sock, char* data, char* buff, uint8_t data_len, const char* tag_start, const char* tag_end) {
	uint8_t tag_start_len = strlen(tag_start);
	uint8_t tag_end_len = strlen(tag_end);

	// clean buffer
	memset(buff, 0, BUFF_LEN);
	// set start tag for data
	strcpy(buff, tag_start);
	// copy data to send
	memcpy(buff+tag_start_len, data, data_len);
	// set end tag for data
	strcpy(buff+tag_start_len+data_len, tag_end);
	// send data
	int err = send(sock, buff,
			data_len+tag_start_len+tag_end_len, 0);

	if (err < 0) {
		ESP_LOGE(TAG,
				"Error occurred during sending: errno %d", errno);
		return ESP_FAIL;
	}
	return ESP_OK;
}

uint8_t my_mpu9250_activate_send_message = 1;

void my_mpu9250_read_data_cycle(mpu9250_handle_t mpu9250_handle) {
	const TickType_t xMaxBlockTime = pdMS_TO_TICKS( 500 );

	ESP_ERROR_CHECK(mpu9250_gyro_set_fsr(mpu9250_handle, INV_FSR_1000DPS));
	ESP_ERROR_CHECK(mpu9250_acc_set_fsr(mpu9250_handle, INV_FSR_8G));
    ESP_ERROR_CHECK(mpu9250_mag_set_mode2_with_precision(mpu9250_handle, INV_MAG_PRECISION_16_BITS));
	ESP_ERROR_CHECK(mpu9250_mag_set_continuous_reading(mpu9250_handle));
	ESP_ERROR_CHECK(mpu9250_discard_messages(mpu9250_handle, 10000));

    int sock =  -1;
    mpu9250_config_data_t config_data = {
    		.accel_fsr = mpu9250_handle->data.accel.fsr,
    		.gyro_fsr = mpu9250_handle->data.gyro.fsr,
    		.mag_precision = mpu9250_handle->data.mag.precision,
    		.mag_drdy = mpu9250_handle->data.mag.drdy,
			.frequenzy_hz = 500
    };

	char buff[BUFF_LEN];
    while(true) {
    	while(true) { // TCP connection
    	    char host_ip[] = HOST_IP;
    	    int addr_family = AF_INET;
    	    int ip_protocol = IPPROTO_IP;
    	    struct sockaddr_in dest_addr;
    	    inet_pton(AF_INET, host_ip, &dest_addr.sin_addr);
    	    dest_addr.sin_family = AF_INET;
    	    dest_addr.sin_port = htons(HOST_PORT);

    		sock = socket(addr_family, SOCK_STREAM, ip_protocol);
    		if (sock < 0) {
    			ESP_LOGE(TAG, "Unable to create socket: errno %d",
    					errno);
    			break;
    		}
    	    ESP_LOGI(TAG, "Socket created, connecting to %s:%d", host_ip, HOST_PORT);

    		int err = connect(sock, (struct sockaddr*) &dest_addr,
    				sizeof(dest_addr));
    		if (err != 0) {
    			ESP_LOGE(TAG, "Socket unable to connect: errno %d",
    					errno);
    			break;
    		}
    	    ESP_LOGI(TAG, "Successfully connected");

    		uint32_t counter = 0;
#ifdef CONFIG_ESP_DATA_CAL
    		mpu9250_cal_message_t data_message;
    		uint16_t data_size = sizeof(mpu9250_cal_message_t);
    		char* data_pointer = (char*)&mpu9250_handle->data.cal_data.data_s_xyz;
#else
    		mpu9250_raw_message_t data_message;
    		uint16_t data_size = sizeof(mpu9250_raw_message_t);
    		char* data_pointer = (char*)&mpu9250_handle->data.raw_data.data_s_xyz;
#endif
			printf("sizeof data: %d\n", data_size);
			printf("sizeof float: %d\n", sizeof(float));
    		while (true) { // read/send cycle
        		memset(&data_message, 0, data_size);
    			counter++;
   				counter %= 100;
    			if( ulTaskNotifyTake( pdTRUE,xMaxBlockTime ) == 1) {
     				ESP_ERROR_CHECK(mpu9250_load_data(mpu9250_handle));
    				if(counter == 0) {
#ifdef CONFIG_ESP_DATA_CAL
    					printf("Gyro .: [%d][%d][%d]\n", mpu9250_handle->data.cal_data.data_s_xyz.gyro_data_x, mpu9250_handle->data.cal_data.data_s_xyz.gyro_data_y, mpu9250_handle->data.cal_data.data_s_xyz.gyro_data_z);
    					printf("Accel : [%5.3f][%5.3f][%5.3f]\n", mpu9250_handle->data.cal_data.data_s_xyz.accel_data_x, mpu9250_handle->data.cal_data.data_s_xyz.accel_data_y, mpu9250_handle->data.cal_data.data_s_xyz.accel_data_z);
    					printf("Mag ..: [%d][%d][%d][%5.3f][%5.3f][%5.3f][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.mag_data_x, mpu9250_handle->data.raw_data.data_s_xyz.mag_data_y, mpu9250_handle->data.raw_data.data_s_xyz.mag_data_z, mpu9250_handle->data.cal_data.data_s_xyz.mag_data_x, mpu9250_handle->data.cal_data.data_s_xyz.mag_data_y, mpu9250_handle->data.cal_data.data_s_xyz.mag_data_z, mpu9250_handle->data.mag.drdy);
#else
    					printf("Gyro .: [%d][%d][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.gyro_data_x, mpu9250_handle->data.raw_data.data_s_xyz.gyro_data_y, mpu9250_handle->data.raw_data.data_s_xyz.gyro_data_z);
    					printf("Accel : [%d][%d][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.accel_data_x, mpu9250_handle->data.raw_data.data_s_xyz.accel_data_y, mpu9250_handle->data.raw_data.data_s_xyz.accel_data_z);
    					printf("Mag ..: [%d][%d][%d]\n", mpu9250_handle->data.raw_data.data_s_xyz.mag_data_x, mpu9250_handle->data.raw_data.data_s_xyz.mag_data_y, mpu9250_handle->data.raw_data.data_s_xyz.mag_data_z);
#endif
    					if(my_mpu9250_activate_send_message) {
            				esp_err_t res = mpu9250_send_message(sock, (char*)&config_data, buff, sizeof(mpu9250_config_data_t), MY_MPU9250_SENSORS_MSG_CONFIG_START, MY_MPU9250_SENSORS_MSG_CONFIG_END);
            				if(res != ESP_OK) {
            					break;
            				}
    					}

    				}
    				if(my_mpu9250_activate_send_message) {
        				// copy data to send
        				memcpy((char*)&data_message.data, data_pointer, data_size);
        				data_message.mag_data_drdy = mpu9250_handle->data.mag.drdy;
        				data_message.timestamp = mpu9250_handle->data.timestamp;
        				esp_err_t res = mpu9250_send_message(sock, (char*)&data_message, buff, data_size, MY_MPU9250_SENSORS_MSG_RAW_DATA_START, MY_MPU9250_SENSORS_MSG_RAW_DATA_END);
        				if(res != ESP_OK) {
        					break;
        				}
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
            ESP_LOGE(TAG, "Shutting down socket and restarting...");
            shutdown(sock, 0);
            close(sock);
        }
        vTaskDelay(pdMS_TO_TICKS(5000));
    }
}

void my_mpu9250_task(void *arg) {
    // Avvia MQTT client
    mqtt_app_start();

	// MPU9250 Handle
	mpu9250_init_t mpu9250;
	mpu9250_handle_t mpu9250_handle = &mpu9250;
	my_mpu9250_task_init(mpu9250_handle);

	// Init MPU9250
	ESP_ERROR_CHECK(mpu9250_init(mpu9250_handle));

	my_mpu9250_read_data_cycle(mpu9250_handle);
}
