
# Autor : Yassine OUAISSA
# Date  : 07 / 11 / 2018 
# Description : receve & visulaze function 
#

import socket as sock
import matplotlib.pyplot as plt
import matplotlib.animation as ani

Data = []

def strToFloat(str):
    for i in range(len(str)):
        str[i] = float(str[i])

def receve_DATA():
    
    socket = sock.socket(sock.AF_INET,sock.SOCK_STREAM)
    socket.bind(('192.168.4.2',1555))
    socket.listen(0)    
    client, adress = socket.accept()
    print("{} Connected ".format(adress))
    rec = client.recv(99999999,sock.MSG_WAITALL)
        #print(rec)
        #if "END" in rec:
         #   break
    data = rec.split(';')
    del data[-1]
    print("[DEBUG]Close Connection")
    client.close()
    socket.close()
    return data
plt.style.use('dark_background')
fig=plt.figure()
ax=fig.add_subplot(1,1,1)


def update_plot(yData):
    try:
        yData = receve_DATA()
        strToFloat(yData)
 	ax.clear()
        ax.plot(yData,'w',linestyle='dotted',marker='o',markerfacecolor='red',markersize=7)
        m  = max(yData)
        idx = yData.index(m)
        plt.annotate('Freq Max = {}'.format(m),xy=(idx,m)
                     ,xytext=(idx+20,m)
                     ,arrowprops = dict(facecolor="blue",shrink=0.05))
    except Exception as e :
        print(e)
       
        
plot_ani = ani.FuncAnimation(fig,update_plot,fargs=(Data),interval = 1)
plt.show()
