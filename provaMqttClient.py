import struct
import paho.mqtt.client as mqtt
import time
import csv
import threading

# Enumerazione dei tipi di messaggio
class MessageType:
    IMU_DATA = 1
    MAGNETO_DATA = 2
    MPU9250_CONFIG_DATA = 3
    CALIBRATION_COMMAND = 4 
    SENSOR_MODEL = 5

# Formati binari
FORMATS = {
    MessageType.IMU_DATA: "I3f3ff32s",  # uint32_t timestamp, float[3] accel, float[3] gyro, float temp
    MessageType.MAGNETO_DATA: "I3f48s", # uint32_t timestamp, float[3] magneto
    MessageType.MPU9250_CONFIG_DATA: "IBBBH",  # uint32_t timestamp, 3x uint8, 1x uint16
    MessageType.SENSOR_MODEL: "I60s", # uint32_t timestamp, 60 bytes blob
}

# File CSV
imu_data_file = "imu-data.csv"
mag_data_file = "mag-data.csv"
MAX_IMU_DATA_MSGS = 100000
imu_data_count = 0
process_active = True

# Double buffers for reading and writing
imu_data_buffers = [[], []]
mag_data_buffers = [[], []]
imu_current_buffer = 0
mag_current_buffer = 0
BUFFER_FLUSH_THRESHOLD = 1000

# Lock and condition for thread synchronization
write_condition = threading.Condition()

# Function to flush a buffer to a file
def flush_buffer(buffer, file_path):
    with open(file_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(buffer)
    buffer.clear()

# Writer thread function
def writer_thread():
    global imu_data_buffers, mag_data_buffers, process_active
    while True:
        with write_condition:
            write_condition.wait()  # Wait for the signal from the reader thread

            if not process_active and not imu_data_buffers[1 - imu_current_buffer] and not mag_data_buffers[1 - mag_current_buffer]:
                break

            # Get the buffers to flush
            imu_write_buffer = imu_data_buffers[1 - imu_current_buffer]
            mag_write_buffer = mag_data_buffers[1 - mag_current_buffer]

            # Flush the buffers
            if imu_write_buffer:
                flush_buffer(imu_write_buffer, imu_data_file)
            if mag_write_buffer:
                flush_buffer(mag_write_buffer, mag_data_file)

# Reader functions to handle data
def write_imu_data(data):
    global imu_data_count, imu_current_buffer, process_active

    imu_data_buffers[imu_current_buffer].append([data[0], *data[1:4], *data[4:7], data[7]])
    imu_data_count += 1

    if imu_data_count >= MAX_IMU_DATA_MSGS:
        switch_buffers()
        return True

    if len(imu_data_buffers[imu_current_buffer]) >= BUFFER_FLUSH_THRESHOLD:
        switch_buffers()

    return False

def write_mag_data(data):
    global mag_current_buffer

    mag_data_buffers[mag_current_buffer].append([data[0], *data[1:4]])

    if len(mag_data_buffers[mag_current_buffer]) >= BUFFER_FLUSH_THRESHOLD:
        switch_buffers()

def switch_buffers():
    global imu_current_buffer, mag_current_buffer

    with write_condition:
        # Switch the current buffers
        imu_current_buffer = 1 - imu_current_buffer
        mag_current_buffer = 1 - mag_current_buffer
        # Notify the writer thread
        write_condition.notify()

# Callback for receiving MQTT messages
def on_message(client, userdata, message):
    global process_active

    if not process_active:
      return
    # print("on message started")
    try:
        # Decode the message type
        payload = message.payload
        msg_type, = struct.unpack("I", payload[:4])

        if msg_type in FORMATS:
            # Decode the message data
            data = struct.unpack(FORMATS[msg_type], payload[4:])
            #print(f"Messaggio ricevuto su {message.topic}:")

            if msg_type == MessageType.IMU_DATA:
                if write_imu_data(data):
                    print("Raggiunto il limite massimo di messaggi IMU_DATA, terminazione dello script.")
                    client.loop_stop()
                    send_calibration_command(client, "stop")
                    process_active = False
 
            elif msg_type == MessageType.MAGNETO_DATA:
                write_mag_data(data)
                print(f"Timestamp: {data[0]}, Magneto: {data[1:4]}")

        else:
            print(f"Tipo di messaggio sconosciuto: {msg_type}")
    except Exception as e:
        print(f"Errore nel parsing del messaggio: {e}")

# MQTT connection callback
def on_connect(client, userdata, flags, rc):
    print("Connesso al broker MQTT con codice di risultato: " + str(rc))

# Function to send calibration commands
# Funzione per inviare il messaggio di calibrazione
def send_calibration_command(client, msg):
    topic = "/imu/calibration/commands"
    timestamp = int(time.time())  # Timestamp corrente (uint32_t)
    msg_type = MessageType.CALIBRATION_COMMAND
    command = msg.encode("utf-8")  # Campo 'data' come stringa

    # Formatta il messaggio binario
    #payload = struct.pack("I", msg_type) + struct.pack("I", timestamp) + command
    payload = command
    # Pubblica il messaggio sul topic
    client.publish(topic, payload)
    print(f"Messaggio di calibrazione inviato su {topic}: timestamp={timestamp}, type={msg_type}, data={msg}")

# Configure MQTT client
broker = "localhost"
port = 1883

client = mqtt.Client()

client.on_message = on_message
client.on_connect = on_connect

client.connect(broker, port, 60)

# Subscribe to topics
topics = ["/imu/calibration/data"]
for topic in topics:
    client.subscribe(topic)

# Start writer thread
writer = threading.Thread(target=writer_thread)
writer.start()

# Start MQTT loop
print("In attesa di messaggi MQTT...")
client.loop_start()

# Keep process active
try:
    # s = True
    while process_active:
        time.sleep(10)
        # if s:
        #     msg = "start"
        # else:
        #     msg = "stop"
        # s = not s
        # send_calibration_command(client, msg)

except KeyboardInterrupt:
    print("Interruzione da tastiera, arresto del client...")
finally:
    with write_condition:
        write_condition.notify()
    writer.join()
    client.loop_stop()
    client.disconnect()
