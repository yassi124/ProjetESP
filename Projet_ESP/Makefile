# Autor :Yassine OUAISSA 
# Date  : 05/11/2018
#

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
