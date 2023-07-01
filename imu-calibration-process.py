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
import imuEllipsoidEstimator as iee
from imuEllipsoidEstimator import estimate_mag_acc

raw_data_queue = Queue()
end_queue = Queue(maxsize=1)

class RawDataFields(ctypes.Structure):
    _fields_ = [
        ("ax",   ctypes.c_short),
        ("ay",   ctypes.c_short),
        ("az",   ctypes.c_short),
        ("temp_data",   ctypes.c_short),
        ("gx",   ctypes.c_short),
        ("gy",   ctypes.c_short),
        ("gz",   ctypes.c_short),
        ("mx",   ctypes.c_short),
        ("my",   ctypes.c_short),
        ("mz",   ctypes.c_short),
        ("ext_1",   ctypes.c_short),
        ("ext_2",   ctypes.c_short),
        ("ext_3",   ctypes.c_short),
        ("ext_4",   ctypes.c_short),
        ("ext_5",   ctypes.c_short),
        ("ext_6",   ctypes.c_short),
        ("ext_7",   ctypes.c_short),
        ("ext_8",   ctypes.c_short),
        ("ext_9",   ctypes.c_short),
        ("mag_drdy",   ctypes.c_uint8),
        ("dummy",   ctypes.c_uint8),
        ("ts",   ctypes.c_uint32),
    ]

def parse_raw_data_payload(raw_data):
  fmt = "<hhhhhhhhhhhhhhhhhhhBBI"
  fmt_size = struct.calcsize(fmt)
  k = RawDataFields();
  k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.ext_1, k.ext_2, k.ext_3, k.ext_4, k.ext_5, k.ext_6, k.ext_7, k.ext_8, k.ext_9, k.mag_drdy, k.dummy,k.ts = struct.unpack(fmt, raw_data[:fmt_size])
  return k
  
def raw_data_renderer_function(file_path, in_queue):
    max_qsize = 0
    data = np.zeros((1000,10), dtype=np.int16)
    
    with open(file_path, "w") as file:
      writer = csv.writer(file, delimiter = ',')
      record = ['AX', 'AY', 'AZ', 'TEMP', 'GX', 'GY', 'GZ', 'MX', 'MY', 'MZ', 'MV', 'TS']
      writer.writerow(record)
      while True:
        try:
          payload = in_queue.get(timeout=10)
          if len(payload) != 44:
            continue
          k = parse_raw_data_payload(payload)
          record = np.array([k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.mag_drdy, k.ts])
          writer.writerow(record)
        except queue.Empty:
          print("Timeout scaduto senza ricevere un messaggio")
          break
    print("raw_data_renderer_function terminato")

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

#def calcola_media_varianza_colonne(dati):
#    num_colonne = len(dati[0])
#    media_colonne = []
#    varianza_colonne = []
#    
#    for colonna in range(num_colonne):
#        valori_colonna = [riga[colonna] for riga in dati]
#        media_colonna = statistics.mean(valori_colonna)
#        varianza_colonna = statistics.variance(valori_colonna)
#        
#        media_colonne.append(media_colonna)
#        varianza_colonne.append(varianza_colonna)
#    
#    return media_colonne, varianza_colonne

def calcola_media_varianza_colonne(dati):
    media_colonne = {}
    varianza_colonne = {}
    median_colonne = {}
    distribution_colonne = {}
    valori_colonne = {}

    colonne_speciali = ['MX', 'MY', 'MZ']
    colonne_normali = [colonna for colonna in dati[0].keys() if colonna not in colonne_speciali]

    for colonna in colonne_normali:
        valori_colonna = [riga[colonna] for riga in dati]
        media_colonna = np.mean(valori_colonna)
        varianza_colonna = np.var(valori_colonna)
        median_colonna = np.median(valori_colonna)
        distribution_colonna = stats.describe(valori_colonna)

        media_colonne[colonna] = media_colonna
        varianza_colonne[colonna] = varianza_colonna
        median_colonne[colonna] = median_colonna
        distribution_colonne[colonna] = distribution_colonna
        valori_colonne[colonna] = valori_colonna

    for colonna in colonne_speciali:
        valori_colonna = [riga[colonna] for riga in dati if riga['MV'] == 1]
        media_colonna = np.mean(valori_colonna)
        varianza_colonna = np.var(valori_colonna)
        median_colonna = np.median(valori_colonna)
        distribution_colonna = stats.describe(valori_colonna)

        media_colonne[colonna] = media_colonna
        varianza_colonne[colonna] = varianza_colonna
        median_colonne[colonna] = median_colonna
        distribution_colonne[colonna] = distribution_colonna
        valori_colonne[colonna] = valori_colonna

    return media_colonne, varianza_colonne, median_colonne, distribution_colonne, valori_colonne


if __name__ == "__main__":

#    # Crea un'istanza del client MQTT
#    client = mqtt.Client()
#
#    # Imposta le funzioni di callback
#    client.on_connect = on_connect
# 
#    # Specifica l'indirizzo IP del broker MQTT
#    broker_address = "localhost"
#
#    # Connettiti al broker MQTT
#    client.connect(broker_address, 1883, 60)
#
#    # Avvia il loop del client MQTT
#    client.loop_start()
#
#    # Creazione del socket TCP per il server
#    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#
#    # Configurazione del socket
#    server_address = ('192.168.1.60', 65432)  # IP del server (vuoto per accettare connessioni da qualsiasi interfaccia di rete), porta del server
#    server_socket.bind(server_address)
#    server_socket.listen(1)  # Limita il numero massimo di connessioni pendenti
#
#    print("Server in ascolto su {}:{}".format(*server_address))
#
#    print("In attesa di connessioni...")
#    client_socket, client_address = server_socket.accept()  # Accetta la connessione del client
#    print("Connessione accettata da {}:{}".format(*client_address))
#
##################################################################################
########## Fase 1: produco dati con sensore fermo
##################################################################################
#    print("please do not move the sensor ...")
#    # Due thread, uno scrive su file e l'altro legge dal socket
#    rd = threading.Thread(target=raw_data_renderer_function, args=("examples/01-imu-raw-data.csv",raw_data_queue))
#    rd.start()
#    sr = threading.Thread(target=read_records, args=(client_socket,raw_data_queue, end_queue, 20000 ))
#    # avvio la produzione di messaggi
#    topic = "/imu/calibration/commands"  # Il topic su cui pubblicare il messaggio
#    message = "start"  # comando l'avvio di produzione dei messaggi
#    client.publish(topic, message)
#    sr.start()
#    # Attendo il segnale di fine elaborazione del thread
#    print("attendo il termine di lettura del socket")
#    end_queue.get()
#    # fermo la produzione di messaggi
#    print("invio messaggio di stop")
#    message = "stop"  # # comando l'interruzione di produzione dei messaggi
#    client.publish(topic, message)
#    print("attendo la chiusura del thread di lettura del socket")
#    sr.join()
#    print("attendo la chiusura del thread di scrittura su file")
#    rd.join()
#    end_queue.queue.clear()
#
##################################################################################
########## Fase 2: produco dati con sensore in lenta rotazione su tutti gli assi
##################################################################################
#    print("please slowly rotate the sensor in all directions ...")
#    # Due thread, uno scrive su file e l'altro legge dal socket
#    rd = threading.Thread(target=raw_data_renderer_function, args=("examples/02-imu-raw-data.csv",raw_data_queue))
#    rd.start()
#    sr = threading.Thread(target=read_records, args=(client_socket,raw_data_queue, end_queue, 100000 ))
#    # avvio la produzione di messaggi
#    topic = "/imu/calibration/commands"  # Il topic su cui pubblicare il messaggio
#    message = "start"  # comando l'avvio di produzione dei messaggi
#    client.publish(topic, message)
#    sr.start()
#    # Attendo il segnale di fine elaborazione del thread
#    print("attendo il termine di lettura del socket")
#    end_queue.get()
#    # fermo la produzione di messaggi
#    print("invio messaggio di stop")
#    message = "stop"  # # comando l'interruzione di produzione dei messaggi
#    client.publish(topic, message)
#    print("attendo la chiusura del thread di lettura del socket")
#    sr.join()
#    print("attendo la chiusura del thread di scrittura su file")
#    rd.join()
#    end_queue.queue.clear()
#
    #################################################################################
    ######### Fase 3: Calcolo media, varianza e deviazione standard 
    #########         dal primo set di dati
    #################################################################################
    file_csv = "examples/01-imu-raw-data.csv"
    dati = pd.read_csv(file_csv)

# FIXME: da implementare partendo da dataframe pandas
#    media_colonne, varianza_colonne, median_colonne, distribution_colonne, valori_colonne = calcola_media_varianza_colonne(dati)
#    print("Media per ogni colonna:")
#    for i, media in enumerate(media_colonne):
#        print(f"Colonna {i+1}: {media}")
#    print("\nVarianza per ogni colonna:")
#    for i, varianza in enumerate(varianza_colonne):
#        print(f"Colonna {i+1}: {varianza}")

    #################################################################################
    ######### Fase 4: Calcolo offset e matrici di correzione 
    #########         dal secondo set di dati utilizzando media, varianza e std
    #################################################################################
    file_csv = "examples/02-imu-raw-data.csv"
    estimator_mag, estimator_acc = estimate_mag_acc(file_csv)
    
    #################################################################################
    ######### Fase 5: Invio offset e matrici di correzione al sensore
    #########         il sensore memorizza il tutto e li adotta per
    #########         produrre dati calibrati
    #################################################################################
    ## TODO: inviare via mqtt
    
    #################################################################################
    ######### Fase 6: Produco dati calibrati con sensore in lenta rotazione
    #################################################################################
    ## TODO: inviare richiesta via mqtt per produrre dati calibrati 
    ##       e switchare parser
    
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

