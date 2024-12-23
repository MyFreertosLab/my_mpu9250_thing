# Edit following two lines to set component requirements (see docs)
set(COMPONENT_REQUIRES )
set(COMPONENT_PRIV_REQUIRES )

set(COMPONENT_SRCS "main.c" "mpu9250_queue_definitions.c" "my_mpu9250_mqtt_definitions.c" "my_mpu9250_sensor_reader_task.c" "my_mpu9250_sensor_writer_task.c" "my_mpu9250_mqtt_receiver_task.c" "my_mpu9250_mqtt_sender_task.c")
set(COMPONENT_ADD_INCLUDEDIRS ".")

register_component()
