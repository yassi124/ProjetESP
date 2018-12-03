
# Autor : Yassine OUAISSA
# Date  : 07 / 11 / 2018 
# Description : Main function for nodeMcu & Micropython
#

import webrepl as web
from machine import I2C, Pin
import ads1x15
import FFT_ESP as fft
import time as t
import socket as sock
import math as m 

PI = 3.14

addr = 72
gain = 1
ADC_RATE = 1
_bufferSize = 1024

i2c = I2C(scl=Pin(5), sda=Pin(4), freq=400000)
ads = ads1x15.ADS1115(i2c, addr, gain)

#ads.set_conv(7, 0)
#ads.read_rev()
ads.alert_start(7,0)
t.sleep_ms(2)

def arange(x):
	T = []
	for i in range(x):
		T.append(i*1.0)
	return T

def adc_read(adc = ads.read):
	tab_DATA = [0.0]*_bufferSize
	for i in range(_bufferSize):
		#t.sleep_ms(ADC_RATE)
		tab_DATA[i] = ads.raw_to_v(adc(7,0))*1.0
		#print("Alert_read function return ",data)
		print("read function return en V  ",tab_DATA[i])
		#print("Alert_rev function return ",ads.read_rev())
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

#Vr = [0.0]*_bufferSize


# Description : Read file data and check if value exist in file_Datas
# Params : nameFile : name of file to check 
#	   value    : value to search 
# 

def ifExist(nameFile,value):
	with open(nameFile) as f: # open file as default mode 'read'
		alldata = f.readlines() # read all fileData
		for v in alldata:	
			data = v.split() # Split all data and make it as array 
	if value in data:
		return True
	return False

######  learning mode  ########
# Description : Write data to file with mode append(write value at the end of the file)
# Params : nameFile : name of file to write in 
#	   value    : value to write 
# 

def writeValue(nameFile,value):
	with open(nameFile,'a+') as f: # write value with append mode 
		f.write(value)
		

if __name__=='__main__':

	while True:
		Vr = adc_read()
# Test FFT Function  #####################################
		#Vr = [m.sin(2*PI*50*(i/_bufferSize) for i in range(x))]
		#for i in range(_bufferSize):
		#	Vr[i] = m.sin(2*PI*100*(i/_bufferSize))
		#print(Vr)
#####################################
		fft.FFT_Pack(Vr)
		print(Vr)
		print("MAX peak is {};time : {}".format(max(Vr),t.time()))
		try:
			Send_FFT_Values(Vr)
		except Exception as e:
			print("[ERROR] {}".format(e))

		finally:
			t.sleep(2)
