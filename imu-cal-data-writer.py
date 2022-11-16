# socket-server.py

import ctypes
import struct
import socket
import threading
from queue import Queue
import numpy as np
import csv
from numpy import dtype

HOST = "192.168.1.64"  
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)
cal_data_queue = Queue(maxsize = 1000)

class CalDataFields(ctypes.Structure):
    _fields_ = [
        ("ax",   ctypes.c_float),
        ("ay",   ctypes.c_float),
        ("az",   ctypes.c_float),
        ("temp_data",   ctypes.c_short),
        ("gx",   ctypes.c_short),
        ("gy",   ctypes.c_short),
        ("gz",   ctypes.c_short),
        ("mx",   ctypes.c_float),
        ("my",   ctypes.c_float),
        ("mz",   ctypes.c_float),
        ("mag_drdy",   ctypes.c_uint8),
        ("dummy",   ctypes.c_uint8),
    ]


def parse_cal_data_payload(cal_data):
  fmt = "<fffhhhhfffBB"
  fmt_size = struct.calcsize(fmt)
  k = CalDataFields();
  k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.mag_drdy, k.dummy = struct.unpack(fmt, cal_data[:fmt_size])
  return k
  
def cal_data_renderer_function(name, in_queue):
    max_qsize = 0
    data = np.zeros((1000,10), dtype=np.double)
    
    with open("imu-cal-data.csv", "w") as file:
      writer = csv.writer(file, delimiter = ',')
      record = ['AX', 'AY', 'AZ', 'TEMP', 'GX', 'GY', 'GZ', 'MX', 'MY', 'MZ', 'MV']
      writer.writerow(record)
      while True:
        payload = in_queue.get()
        len_payload = len(payload)
        if len_payload != 36:
          print(f"len(payoad) by {len_payload}")
          continue
        k = parse_cal_data_payload(payload)
        record = np.array([k.ax, k.ay, k.az, k.temp_data, k.gx, k.gy, k.gz, k.mx, k.my, k.mz, k.mag_drdy], dtype=np.double)
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
                 print("Wrong Prefix! {}".format(prefix))
                 data_prev = ''.encode()
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
    rd = threading.Thread(target=cal_data_renderer_function, args=("cal_data_renderer",cal_data_queue, ))
    rd.start()
    sr = threading.Thread(target=socket_receiver_function, args=("receiver",cal_data_queue, ))
    sr.start()
    sr.join()
    
