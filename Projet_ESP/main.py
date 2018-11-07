
# Autor : Yassine OUAISSA
# Date  : 07 / 11 / 2018 
# Description : Main function for nodeMcu & Micropython
#

import webrepl as web
import machine as m 
import FFT_ESP as fft
import time as t
import socket as sock


def adc_read():
	
	tab_DATA = []
	adc = m.ADC(0)
	for i in range(256):
		tab_DATA.append(adc.read()*(3.3/1023.0))
		#t.sleep_us(1)
	return tab_DATA

# Description : Send all datas via socket for visu
# Params : 
#          rData : Data to send 

def Send_FFT_Values(rData):
	
	#rData = list(range(256))
	#print(rData)
	socket = sock.socket(sock.AF_INET,sock.SOCK_STREAM)
	socket.connect(("192.168.4.2",1555))
	
	dataToSend =''.join("{};".format(e) for e in rData)
	socket.sendall("{}END".format(dataToSend))
	#t.sleep_us(20)
	#print("[DEBUG] Closing socket")
	socket.close()

if __name__=='__main__':

	while True:
		Vr = adc_read()
		fft.FFT_Pack(Vr)
		
		try:
			Send_FFT_Values(Vr)
		except Exception as e:
			print("[ERROR] {}".format(e))

		finally:
			t.sleep(2)
