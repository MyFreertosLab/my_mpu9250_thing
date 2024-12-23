/*
 * mpu9250_queue_definitions.c
 *
 *  Created on: 22 dic 2024
 *      Author: andrea
 */

#include <mpu9250_queue_definitions.h>

// Code globali
QueueHandle_t mpu9250_output_queue = NULL;
QueueHandle_t mpu9250_input_queue = NULL;

