/*
 * my_mpu9250_mqtt_task.c
 *
 *  Created on: 22 dic 2024
 *      Author: andrea
 */
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <esp_err.h>
#include <esp_log.h>
#include <mpu9250_queue_definitions.h>
#include <my_mpu9250_mqtt_definitions.h>

#define MQTT_TOPIC_COMMANDS "/imu/calibration/commands"
#define MQTT_TOPIC_MAG_MODEL "/imu/calibration/mag/model"
#define MQTT_TOPIC_ACC_MODEL "/imu/calibration/acc/model"
#define MQTT_TOPIC_GYRO_MODEL "/imu/calibration/gyro/model"

static const char *TAG = "my_mpu9250_mqtt_receiver";

// NOTE:
// la frequenza di ricezione è bassa.
// non serve allocare dinamicamente il messaggio
// per cui lo alloco qui staticamente
static sensor_message_t msg = {0};

esp_err_t mqtt_event_handler_cb(esp_mqtt_event_handle_t event)
{
    esp_mqtt_client_handle_t client = event->client;
    int msg_id;
    ESP_LOGI(TAG, "Received MQTT event %d", event->event_id);
    uint32_t timestamp = (uint32_t)(esp_timer_get_time()/1000);
    switch (event->event_id) {
        case MQTT_EVENT_CONNECTED:
            ESP_LOGI(TAG, "Connesso al broker MQTT");
            msg_id = esp_mqtt_client_subscribe(client, MQTT_TOPIC_COMMANDS, 0);
            ESP_LOGI(TAG, "Sottoscrizione al topic %s. msg_id=%d",MQTT_TOPIC_COMMANDS, msg_id);
            msg_id = esp_mqtt_client_subscribe(client, MQTT_TOPIC_MAG_MODEL, 0);
            ESP_LOGI(TAG, "Sottoscrizione al topic %s. msg_id=%d",MQTT_TOPIC_MAG_MODEL, msg_id);
            msg_id = esp_mqtt_client_subscribe(client, MQTT_TOPIC_ACC_MODEL, 0);
            ESP_LOGI(TAG, "Sottoscrizione al topic %s. msg_id=%d",MQTT_TOPIC_ACC_MODEL, msg_id);
            msg_id = esp_mqtt_client_subscribe(client, MQTT_TOPIC_GYRO_MODEL, 0);
            ESP_LOGI(TAG, "Sottoscrizione al topic %s. msg_id=%d",MQTT_TOPIC_GYRO_MODEL, msg_id);
            break;

        case MQTT_EVENT_DISCONNECTED:
            ESP_LOGI(TAG, "Disconnesso dal broker MQTT");
            break;

        case MQTT_EVENT_DATA:
            ESP_LOGI(TAG, "Dati ricevuti sul topic %.*s : %.*s",
                     event->topic_len, event->topic, event->data_len, event->data);

             //while(msg.type != 0) {} // attendo che il precedente messaggio sia consumato (type == 0)
             memset(&msg, '\0', sizeof(sensor_message_t)); // elimino eventuali valori residui nel buffer

             if(strncmp(event->topic, MQTT_TOPIC_COMMANDS, event->topic_len) == 0) {
               if (strncmp(event->data, "start", strlen("start")) == 0) {
                 msg.type = MESSAGE_TYPE_COMMAND_START;
                 msg.timestamp = timestamp;
               } else if (strncmp(event->data, "stop", strlen("stop")) == 0)  {
                 msg.type = MESSAGE_TYPE_COMMAND_STOP;
                 msg.timestamp = timestamp;
               }
             } else if(strncmp(event->topic, MQTT_TOPIC_MAG_MODEL, event->topic_len) == 0) {
            	msg.type = MESSAGE_TYPE_MAG_MODEL;
                msg.timestamp = timestamp;
                memcpy((void*)&msg.data.sensor_model, event->data, event->data_len);
            } else if(strncmp(event->topic, MQTT_TOPIC_ACC_MODEL, event->topic_len) == 0) {
            	msg.type = MESSAGE_TYPE_ACC_MODEL;
                msg.timestamp = timestamp;
                memcpy((void*)&msg.data.sensor_model, event->data, event->data_len);
            } else if(strncmp(event->topic, MQTT_TOPIC_GYRO_MODEL, event->topic_len) == 0) {
            	msg.type = MESSAGE_TYPE_GYRO_MODEL;
                msg.timestamp = timestamp;
                memcpy((void*)&msg.data.sensor_model, event->data, event->data_len);
            }

            // Invia il messaggio mqtt ricevuto alla coda
            if (xQueueSend(mpu9250_input_queue, &msg, pdMS_TO_TICKS(100)) != pdPASS) {
                ESP_LOGW(TAG, "Coda piena, impossibile inviare il messaggio");
                memset(&msg, '\0', sizeof(sensor_message_t)); // libero il buffer
            }
            break;

//        case MQTT_EVENT_DATA:
//            ESP_LOGI(TAG, "Dati ricevuti sul topic %.*s : %.*s",
//                     event->topic_len, event->topic, event->data_len, event->data);
//            if(strncmp(event->topic, MQTT_TOPIC_COMMANDS, event->topic_len) == 0) {
//                if (strncmp(event->data, "start", strlen("start")) == 0) {
////                	my_mpu9250_activate_send_message = 1;
//                } else if (strncmp(event->data, "stop", strlen("stop")) == 0)  {
////                	my_mpu9250_activate_send_message = 0;
//                }
//                if (strncmp(event->data, "raw_start", strlen("raw_start")) == 0) {
////                	my_mpu9250_activate_send_message = 1;
////                	message_type = MESSAGE_TYPE_RAW;
//                    ESP_LOGI(TAG, "Message type RAW selected");
//                } else if (strncmp(event->data, "cal_start", strlen("cal_start")) == 0)  {
////                	my_mpu9250_activate_send_message = 1;
//                    ESP_LOGI(TAG, "Message type CAL selected");
////                	message_type = MESSAGE_TYPE_CAL;
//                }
//            } else if(strncmp(event->topic, MQTT_TOPIC_MAG_MODEL, event->topic_len) == 0) {
////            	mpu9250_mag_save_calibration_params(mpu9250_handle, event->data, event->data_len);
////            	mpu9250_mag_load_calibration_params(mpu9250_handle);
//
//            } else if(strncmp(event->topic, MQTT_TOPIC_ACC_MODEL, event->topic_len) == 0) {
////               	mpu9250_acc_save_calibration_params(mpu9250_handle, event->data, event->data_len);
////            	mpu9250_acc_load_calibration_params(mpu9250_handle);
//            } else if(strncmp(event->topic, MQTT_TOPIC_GYRO_MODEL, event->topic_len) == 0) {
////               	mpu9250_gyro_save_calibration_params(mpu9250_handle, event->data, event->data_len);
////            	mpu9250_gyro_load_calibration_params(mpu9250_handle);
//            }
//            break;

        default:
            break;
    }
    return ESP_OK;
}

void my_mpu9250_mqtt_receiver_task(void *arg)
{

    ESP_LOGI(TAG, "Task STARTED!");

    // Creazione code per invio/ricezione dati
    mpu9250_output_queue = xQueueCreate(10, sizeof(sensor_message_t));
    mpu9250_input_queue = xQueueCreate(10, sizeof(sensor_message_t));

    if (mpu9250_output_queue == NULL || mpu9250_input_queue == NULL) {
        ESP_LOGE("APP", "Errore nella creazione delle code");
        return;
    }
    ESP_LOGI(TAG, "QUEUES CREATED!");
    // TODO: configurare in flash memory
    esp_mqtt_client_config_t mqtt_cfg = {
        .uri = "mqtt://192.168.1.61:1883", // Indirizzo del broker MQTT
        .event_handle = mqtt_event_handler_cb,
        .buffer_size = 2048,
    };

    ESP_LOGI(TAG, "STARTING MQTT CLIENT ....");
    mqtt_client = esp_mqtt_client_init(&mqtt_cfg);
    ESP_ERROR_CHECK(esp_mqtt_client_start(mqtt_client));
    ESP_LOGI(TAG, "MQTT CLIENT STARTED!!");

    while(true) {
    	vTaskDelay(pdMS_TO_TICKS(20000));
    }

}
