# echo-server.py

import socket

HOST = "192.168.1.64"  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.bind((HOST, PORT))
    s.listen()
    conn, addr = s.accept()
    with conn:
        print(f"Connected by {addr}")
        while True:
            data = ''.encode()
            data1 = conn.recv(1024)
            if not data1:
                break
            data = data + data1
            rl = len(data)
            start_pos = 0
            end_pos = 0
            while(start_pos < rl):
              end_pos = start_pos+4
              prefix = data[start_pos:end_pos].decode()
              if prefix == "<RW>":
                 suffix = '</RW>'
              elif prefix == "<CN>":
                 suffix = '</CN>'
              print("1. start={}, end={}, prefix={}, suffix={}".format(start_pos, end_pos, prefix, suffix))

              start_pos = end_pos
              end_pos=data[start_pos:].find(suffix.encode())+start_pos
              if end_pos < start_pos:
                 data=data[start_pos:]
                 break

              rw_record = data[start_pos:end_pos]
              suffix=data[end_pos:end_pos+5].decode()
              print("2. start={}, end={}, prefix={}, suffix={}".format(start_pos, end_pos, prefix, suffix))

              start_pos = end_pos
              end_pos = start_pos + 5
              suffix=data[start_pos:end_pos].decode()
              rwl=len(rw_record)
              print("3. start={}, end={}, prefix={}, suffix={}, rwlen={}, len={}".format(start_pos, end_pos, prefix, suffix, rwl, rl))
              start_pos = end_pos+1
