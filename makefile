CC := gcc
EXECUTABLE := image_compressor.exe

.PHONY: all executable clean

all: executable

executable:
	$(CC) image_compressor.c -o $(EXECUTABLE)

clean:
	del $(wildcard *_compressed.bmp)
	del $(EXECUTABLE)