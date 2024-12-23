/*
 * my_mpu9250_mqtt_sender_task.c
 *
 *  Created on: 22 dic 2024
 *      Author: andrea
 */


#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <esp_log.h>
#include <mpu9250_queue_definitions.h>
#include <my_mpu9250_mqtt_definitions.h>

static const char *TAG = "my_mpu9250_mqtt_sender";
#define MQTT_TOPIC_SENSOR_DATA "/imu/calibration/data"

void my_mpu9250_mqtt_sender_task(void *arg)
{
    sensor_message_t sensor_msg_in;
    sensor_message_t sensor_msg_out;

    while (true) {
        // Attendi dati dalla coda del sensore
        if (xQueueReceive(mpu9250_output_queue, &sensor_msg_in, portMAX_DELAY) == pdTRUE) {
            //ESP_LOGI(TAG, "Preparazione invio dati: %d", sensor_msg_in.type);
            memset(&sensor_msg_out, '\0', sizeof(sensor_message_t)); // preparo il buffer in uscita
            memcpy((char*)&sensor_msg_out, (char*)&sensor_msg_in, sizeof(sensor_message_t)); // copio i dat da inviare
            memset(&sensor_msg_in, '\0', sizeof(sensor_message_t)); // libero il buffer di ingresso
            
             // Invia i dati tramite MQTT
            // il messaggio viene bufferizzato
            int msg_id = esp_mqtt_client_publish(mqtt_client, MQTT_TOPIC_SENSOR_DATA, (char*)&sensor_msg_out, sizeof(sensor_message_t), 0, 0);
            if (msg_id == -1) {
              ESP_LOGE(TAG, "Errore nell'invio del messaggio mqtt");
            }
        }
    }
}

