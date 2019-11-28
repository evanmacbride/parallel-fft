SOURCES = src/main.cpp src/fft.cpp
ACC_SOURCES = src/main.cpp src/acc-fft.cpp
SERIAL_TARGET = serial-fft
OMP_TARGET = omp-fft
ACC_TARGET = acc-fft
CFLAGS = -lm -O2 -g -pg
CC = g++

serial:
	$(CC) $(SOURCES) $(CFLAGS) -o $(SERIAL_TARGET)

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(OMP_TARGET)

acc:
	$(CC) $(ACC_SOURCES) $(CFLAGS) -fopenacc -o $(ACC_TARGET)

all:
	$(CC) $(SOURCES) $(CFLAGS) -o $(SERIAL_TARGET)
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(OMP_TARGET)
	$(CC) $(ACC_SOURCES) $(CFLAGS) -fopenacc -o $(ACC_TARGET)

run:
	./$(SERIAL_TARGET)

clean:
	rm $(SERIAL_TARGET)
	rm $(OMP_TARGET)
	rm $(ACC_TARGET)
