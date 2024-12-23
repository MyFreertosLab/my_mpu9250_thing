import struct
import paho.mqtt.client as mqtt
import time

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

# Callback per la ricezione di messaggi MQTT
def on_message(client, userdata, message):
    print("on message started")
    try:
        # Decodifica il tipo di messaggio dal payload
        payload = message.payload
        print(len(payload))
        msg_type, = struct.unpack("I", payload[:4])  # Legge il tipo di messaggio (4 byte iniziali)
        if msg_type in FORMATS:
            # Decodifica il resto del messaggio in base al tipo
            data = struct.unpack(FORMATS[msg_type], payload[4:])
            print(f"Messaggio ricevuto su {message.topic}:")
            if msg_type == MessageType.IMU_DATA:
                print(f"Timestamp: {data[0]}, Accel: {data[1:4]}, Gyro: {data[4:7]}, Temp: {data[7]}")
            elif msg_type == MessageType.MAGNETO_DATA:
                print(f"Timestamp: {data[0]}, Magneto: {data[1:4]}")
            elif msg_type == MessageType.MPU9250_CONFIG_DATA:
                print(f"Timestamp: {data[0]}, Config: {data[1:]}")
            elif msg_type == MessageType.SENSOR_MODEL:
                print(f"Timestamp: {data[0]}, Sensor Model Blob: {data[1]}")
        else:
            print(f"Tipo di messaggio sconosciuto: {msg_type}")
    except Exception as e:
        print(f"Errore nel parsing del messaggio: {e}")


# Funzione di callback chiamata quando il client si connette al broker MQTT
def on_connect(client, userdata, flags, rc):
    print("Connesso al broker MQTT con codice di risultato: " + str(rc))


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

# Configura il client MQTT
broker = "localhost"  # Sostituisci con il tuo broker MQTT
port = 1883

client = mqtt.Client()

# Configura i callback
client.on_message = on_message

# Imposta le funzioni di callback
client.on_connect = on_connect

# Connessione al broker
client.connect(broker, port, 60)

# Sottoscrizione ai topic
topics = ["/imu/calibration/data"]
for topic in topics:
    client.subscribe(topic)

# Invia il messaggio di calibrazione
#send_calibration_command(client)

# Avvia il loop
print("In attesa di messaggi MQTT...")
client.loop_start()

# Mantieni il processo attivo
try:
    s = True
    while True:
        time.sleep(10)  # Attendi per evitare che il thread principale consumi CPU
        if s:
            msg = "start"
        else:
            msg = "stop"
        s = not s
        send_calibration_command(client, msg)

except KeyboardInterrupt:
    print("Interruzione da tastiera, arresto del client...")
    client.loop_stop()
    client.disconnect()