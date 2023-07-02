import ctypes
import struct
import socket
import select
import threading
import queue
from queue import Queue
import numpy as np
import csv
from numpy import dtype
import paho.mqtt.client as mqtt
from scipy import stats
import pandas as pd
import imu_ellipsoid_estimator as iee
from imu_ellipsoid_estimator import estimate_mag_acc

cal_data_queue = Queue()
end_queue = Queue(maxsize=1)

class CalDataFields(ctypes.Structure):
    _fields_ = [
        ("ax",   ctypes.c_float),
        ("ay",   ctypes.c_float),
        ("az",   ctypes.c_float),
        ("temp_data",   ctypes.c_float),
        ("gx",   ctypes.c_float),
        ("gy",   ctypes.c_float),
        ("gz",   ctypes.c_float),
        ("mx",   ctypes.c_float),
        ("my",   ctypes.c_float),
        ("mz",   ctypes.c_float),
        ("mag_drdy",   ctypes.c_uint8),
        ("dummy",   ctypes.c_uint8),
        ("dummy1",   ctypes.c_uint16),
        ("ts",   ctypes.c_uint32),
    ]


def parse_cal_data_payload(cal_data):
  fmt = "<ffffffffffBBhI"
  fmt_size = struct.calcsize(fmt)
  k = CalDataFields();
  k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.mag_drdy, k.dummy, k.dummy1, k.ts = struct.unpack(fmt, cal_data[:fmt_size])
  return k
  
def cal_data_renderer_function(file_path, in_queue):
    max_qsize = 0
    data = np.zeros((1000,10), dtype=np.double)
    
    with open(file_path, "w") as file:
      writer = csv.writer(file, delimiter = ',')
      record = ['AX', 'AY', 'AZ', 'TEMP', 'GX', 'GY', 'GZ', 'MX', 'MY', 'MZ', 'MV', 'TS']
      writer.writerow(record)
      while True:
        payload = in_queue.get(timeout=10)
        len_payload = len(payload)
        if len_payload != 48:
          print(f"len(cal_payoad) by {len_payload}")
          continue
        k = parse_cal_data_payload(payload)
        record = np.array([k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.mag_drdy, k.ts])
        writer.writerow(record)


def read_records(sock, out_queue, end_queue, num_records):
    buffer = b""
    tot_records = 0
    while True:
      ready_to_read, _, _ = select.select([sock], [], [], 5)
      if sock in ready_to_read:
        try:

          data = sock.recv(1024)  # Ricevi i dati dal socket
          if not data:
              break

          buffer += data

          # Verifica se ci sono record completi nel buffer
          start_index = buffer.find(b"<RW>")
          end_index = buffer.find(b"</RW>")
          while start_index != -1 and end_index != -1:
            record = buffer[start_index + len(b"<RW>"):end_index]
            if tot_records < num_records:
              out_queue.put(record)
              tot_records +=1
            else:
              if end_queue.empty():
                end_queue.put(f"{tot_records}")   
            buffer = buffer[end_index + len(b"</RW>"):]  # Rimuovi il record dal buffer
            start_index = buffer.find(b"<RW>")
            end_index = buffer.find(b"</RW>")

        except socket.timeout:
            print("Timeout scaduto durante la ricezione")
            break
        except socket.error as e:
            print("Errore durante la ricezione:", e)
            break
      else:
          print("Timeout scaduto prima che il socket fosse pronto per la lettura")
          break

    print("read_records terminato")

# Funzione di callback chiamata quando il client si connette al broker MQTT
def on_connect(client, userdata, flags, rc):
    print("Connesso al broker MQTT con codice di risultato: " + str(rc))

# Lettura files csv
def leggi_dati_da_csv(file_path, delimiter=','):
    dati = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=delimiter)
        for riga in reader:
            dati.append([float(valore) for valore in riga])
    return dati

if __name__ == "__main__":

    # Crea un'istanza del client MQTT
    client = mqtt.Client()

    # Imposta le funzioni di callback
    client.on_connect = on_connect
 
    # Specifica l'indirizzo IP del broker MQTT
    broker_address = "localhost"

    # Connettiti al broker MQTT
    client.connect(broker_address, 1883, 60)

    # Avvia il loop del client MQTT
    client.loop_start()

    # Creazione del socket TCP per il server
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    # Configurazione del socket
    server_address = ('192.168.1.60', 65432)  # IP del server (vuoto per accettare connessioni da qualsiasi interfaccia di rete), porta del server
    server_socket.bind(server_address)
    server_socket.listen(1)  # Limita il numero massimo di connessioni pendenti

    print("Server in ascolto su {}:{}".format(*server_address))

    print("In attesa di connessioni...")
    client_socket, client_address = server_socket.accept()  # Accetta la connessione del client
    print("Connessione accettata da {}:{}".format(*client_address))

    print("please do not move the sensor ...")
    # Due thread, uno scrive su file e l'altro legge dal socket
    rd = threading.Thread(target=cal_data_renderer_function, args=("examples/01-imu-cal-data.csv",cal_data_queue))
    rd.start()
    sr = threading.Thread(target=read_records, args=(client_socket,cal_data_queue, end_queue, 100000 ))
    # avvio la produzione di messaggi
    topic = "/imu/calibration/commands"  # Il topic su cui pubblicare il messaggio
    message = "cal_start"  # comando l'avvio di produzione dei messaggi
    client.publish(topic, message)
    sr.start()
    # Attendo il segnale di fine elaborazione del thread
    print("attendo il termine di lettura del socket")
    end_queue.get()
    # fermo la produzione di messaggi
    print("invio messaggio di stop")
    message = "stop"  # # comando l'interruzione di produzione dei messaggi
    client.publish(topic, message)
    print("attendo la chiusura del thread di lettura del socket")
    sr.join()
    print("attendo la chiusura del thread di scrittura su file")
    rd.join()
    end_queue.queue.clear()    
    #################################################################################
    ######### Fase 7: Verifico calibrazione e corretto allineamento assi accel e mag
    #################################################################################
    ## TODO: utilizzare script R
    
    #################################################################################
    ######### Fase 8: Invio correzioni e trasformazioni per allineamento assi
    #########         il sensore memorizza il tutto e li adotta per
    #########         produrre attitude con fusion dei sensori
    #################################################################################
    ## TODO: utilizzare script R
    
    #################################################################################
    #################################### F I N E ####################################
    #################################################################################
    
    
    # Fermo mqtt
    client.loop_stop()
    client.disconnect()
    
    # Chiusura del socket del client
    client_socket.close()

