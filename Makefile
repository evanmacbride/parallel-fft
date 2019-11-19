SOURCES = src/serial-fft.c
TARGET = serial
CFLAGS = -std=c99 -lm -O2
CC = gcc

serial:
	$(CC) $(SOURCES) $(CFLAGS) -o $(TARGET)

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(TARGET)

run:
	./$(TARGET)

clean:
	rm $(TARGET)
