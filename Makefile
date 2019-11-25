SOURCES = src/main.c src/fft.c
SERIAL_TARGET = serial-fft
OMP_TARGET = omp-fft
ACC_TARGET = acc-fft
CFLAGS = -std=c99 -lm -O2 -g
CC = gcc

serial:
	$(CC) $(SOURCES) $(CFLAGS) -o $(SERIAL_TARGET)

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(OMP_TARGET)

acc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(ACC_TARGET)

all:
	$(CC) $(SOURCES) $(CFLAGS) -o $(SERIAL_TARGET)
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(OMP_TARGET)
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o $(ACC_TARGET)

run:
	./$(SERIAL_TARGET)

clean:
	rm $(SERIAL_TARGET)
	rm $(OMP_TARGET)
	rm $(ACC_TARGET)
