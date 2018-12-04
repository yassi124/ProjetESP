# Autor :Yassine OUAISSA 
# Date  : 05/11/2018
#

####################################################################################################################################
#Warning : Befor execute this Makefile, you will install those depencies : 							   #
#********	Esptool (tool for flashing the ESP-Micro-Controller ), picocom(interfacing with the UART bus)			   #
#																   #
#      *** Installation Linux ***												   #
#    Python Pip : sudo apt-get install python-setuptools python-pip / python3-XXXX(for python3)					   #
#    ESPTOOL : with pip : sudo pip install esptool / pip3 install esptool (python3-pip)						   #
#    PICOCOM : sudo apt-get install picocom
# 
#	Installation with one command : 
#   					sudo apt-get install python-pip python-setuptools picocom; sudo pip install esptool        #
####################################################################################################################################

PORT ?=/dev/ttyUSB0
BAUD ?=115200
ESPTOOL ?=esptool.py

MICRO_BIN_DIR ?=./Micropython_bin
MICRO_BIN_NAME ?=esp8266-20180511
MICRO_BIN_VERSION ?=1.9.4
MICRO_BIN=$(MICRO_BIN_DIR)/$(MICRO_BIN_NAME)-v$(MICRO_BIN_VERSION).bin
BIN =''
clear_flash:
	 $(ESPTOOL) --baud $(BAUD) --port $(PORT) erase_flash
flash_microPython:
	$(ESPTOOL) --baud $(BAUD) --port $(PORT) write_flash --flash_size=detect 0 $(MICRO_BIN)
flash:
	$(ESPTOOL) --baud $(BAUD) --port $(PORT) write_flash --flash_size=detect 0 $(BIN)

connect:
	picocom $(PORT) -b $(BAUD)
