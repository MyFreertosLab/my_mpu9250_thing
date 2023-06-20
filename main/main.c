#include "freertos/FreeRTOS.h"
#include "esp_wifi.h"
#include "esp_system.h"
#include "esp_event.h"
#include "nvs_flash.h"
#include "esp_log.h"
#include "my_mpu9250_task.h"

static const char *TAG = "my_mpu9250_thing_main";
static uint8_t wifi_ready = 0;


static void got_ip_event_handler(void* arg, esp_event_base_t event_base, int32_t event_id, void* event_data)
{
    ip_event_got_ip_t* event = (ip_event_got_ip_t*)event_data;
    ESP_LOGI(TAG, "Indirizzo IP ottenuto: " IPSTR, IP2STR(&event->ip_info.ip));
    ESP_LOGI(TAG, "WIFI READY!");
    wifi_ready = 1;
}


static void wifi_event_handler(void* arg, esp_event_base_t event_base,
                                    int32_t event_id, void* event_data)
{
	switch(event_id) {
	case WIFI_EVENT_STA_CONNECTED: {
//        wifi_event_sta_connected_t* event = (wifi_event_sta_connected_t*) event_data;
        ESP_LOGI(TAG, "CONNECTED");
		break;
	}
	case WIFI_EVENT_STA_DISCONNECTED: {
//        wifi_event_sta_disconnected_t* event = (wifi_event_sta_disconnected_t*) event_data;
        ESP_LOGI(TAG, "DISCONNECTED");
        esp_err_t ret = esp_wifi_connect();
        if(ret != ESP_OK) {
        	ESP_LOGE(TAG, "Sorry!! not connected to wifi ...");
        } else {
        	ESP_LOGI(TAG, "Connected!");
        }
		break;
	}
	case WIFI_EVENT_STA_START: {
        ESP_LOGI(TAG, "START");
        esp_err_t ret = esp_wifi_connect();
        if(ret != ESP_OK) {
        	ESP_LOGE(TAG, "Sorry!! not connected to wifi ...");
        } else {
        	ESP_LOGI(TAG, "Connected!");
        }
		break;
	}
	case WIFI_EVENT_STA_STOP: {
        ESP_LOGI(TAG, "STOP");
		break;
	}
	case SYSTEM_EVENT_WIFI_READY: {
        ESP_LOGI(TAG, "WIFI READY!");
        wifi_ready = 1;
		break;
	}
	default: {
        ESP_LOGI(TAG, "WIFI received event [%d]", event_id);
		break;
	}
	}
}


void app_main(void)
{
    ESP_LOGI(TAG, "[APP] Startup..");
    ESP_LOGI(TAG, "[APP] Free memory: %d bytes", esp_get_free_heap_size());
    ESP_LOGI(TAG, "[APP] IDF version: %s", esp_get_idf_version());
    esp_log_level_set("*", ESP_LOG_INFO);
    esp_log_level_set("TRANS_TCP", ESP_LOG_DEBUG);

    nvs_flash_init();
    tcpip_adapter_init();

    esp_netif_init();
    esp_event_loop_create_default();
    esp_event_handler_register(IP_EVENT, IP_EVENT_STA_GOT_IP, &got_ip_event_handler, NULL);
//    ESP_ERROR_CHECK( esp_event_loop_init(event_handler, NULL) );

    wifi_init_config_t cfg = WIFI_INIT_CONFIG_DEFAULT();
    ESP_ERROR_CHECK( esp_wifi_init(&cfg) );

    ESP_ERROR_CHECK( esp_wifi_set_storage(WIFI_STORAGE_RAM) );
    ESP_ERROR_CHECK( esp_wifi_set_mode(WIFI_MODE_STA) );

    ESP_ERROR_CHECK(esp_event_handler_instance_register(WIFI_EVENT,
                                                        ESP_EVENT_ANY_ID,
                                                        &wifi_event_handler,
														NULL,
                                                        NULL));

    wifi_config_t sta_config = {
        .sta = {
        	.channel = 12,
			.pmf_cfg.capable = false,
            .ssid = CONFIG_ESP_WIFI_SSID,
			.scan_method = WIFI_ALL_CHANNEL_SCAN,
            .password = CONFIG_ESP_WIFI_PASSWORD,
            .bssid_set = false
        }
    };
    ESP_ERROR_CHECK( esp_wifi_set_config(WIFI_IF_STA, &sta_config) );
    ESP_ERROR_CHECK( esp_wifi_start() );
    ESP_LOGI(TAG, "wifi_init_softap finished. SSID:%s password:%s",
    		CONFIG_ESP_WIFI_SSID, CONFIG_ESP_WIFI_PASSWORD);

    while(wifi_ready == 0) {
    	vTaskDelay(pdMS_TO_TICKS(2000));
    }
    xTaskCreate(my_mpu9250_task, "my_mpu9250_task", 4096, NULL, 5, NULL);

}

