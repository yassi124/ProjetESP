import socket as sock
import matplotlib.pyplot as p

def receve_DATA():
    
    socket = sock.socket(sock.AF_INET,sock.SOCK_STREAM)
    socket.bind(('192.168.4.2',1555))
    
    while True:
        socket.listen(0)
        client, adress = socket.accept()
        print("{} Connected ".format(adress))
        rec = client.recv(9999999)
        print(rec)
        if "END" in rec:
            break
    data = rec.split(';')
    del data[-1]
    
    print("[DEBUG]Close Connection")
    client.close()
    socket.close()
    return data

print(receve_DATA())
