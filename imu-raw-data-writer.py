# socket-server.py

import ctypes
import struct
import socket
import threading
from queue import Queue
import numpy as np
import csv
from numpy import dtype

HOST = "192.168.1.61"  
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)
raw_data_queue = Queue(maxsize = 1000)

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
  
def raw_data_renderer_function(name, in_queue):
    max_qsize = 0
    data = np.zeros((1000,10), dtype=np.int16)
    
    with open("imu-raw-data.csv", "w") as file:
      writer = csv.writer(file, delimiter = ',')
      record = ['AX', 'AY', 'AZ', 'TEMP', 'GX', 'GY', 'GZ', 'MX', 'MY', 'MZ', 'MV', 'TS']
      writer.writerow(record)
      while True:
        payload = in_queue.get()
        if len(payload) != 44:
          print("len(payload)", len(payload))
          continue
        k = parse_raw_data_payload(payload)
        record = np.array([k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.mag_drdy, k.ts])
        writer.writerow(record)

def socket_receiver_function(name, out_queue):
  with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.bind((HOST, PORT))
    s.listen()
    conn, addr = s.accept()
    with conn:
        print(f"Connected by {addr}")
        data_prev = ''.encode()
        while True:
            data_curr = conn.recv(1024)
            if not data_curr:
                break
            data = data_prev + data_curr
            rl = len(data)
            start_pos = 0
            end_pos = 0
            while(start_pos < rl):
              end_pos = start_pos+4
              prefix = data[start_pos:end_pos]
              if prefix == "<RW>".encode():
                 suffix = '</RW>'.encode()
              elif prefix == "<CN>".encode():
                 suffix = '</CN>'.encode()
              else:
                 print("Wrong Prefix! {} start: {}, end: {}, data size: {}".format(prefix,start_pos,end_pos,len(data)))
                 data_prev = prefix
                 break

              payload_start_pos = end_pos
              end_pos=data[payload_start_pos:].find(suffix)+payload_start_pos
              if end_pos < payload_start_pos:
                 data_prev=data[start_pos:]
                 break
              payload = data[payload_start_pos:end_pos]
              
              suffix=data[end_pos:end_pos+5]

              suffix_start_pos = end_pos
              end_pos = suffix_start_pos + 5
              suffix=data[suffix_start_pos:end_pos]
              rwl=len(payload)
              out_queue.put(payload)
              start_pos = end_pos
              data_prev = ''.encode()


if __name__ == "__main__":
    rd = threading.Thread(target=raw_data_renderer_function, args=("raw_data_renderer",raw_data_queue, ))
    rd.start()
    sr = threading.Thread(target=socket_receiver_function, args=("receiver",raw_data_queue, ))
    sr.start()
    sr.join()
    
