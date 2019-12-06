SOURCES = src/main.cpp src/fft.cpp
ACC_SOURCES = src/main.cpp src/acc-fft.cpp
SERIAL_TARGET = serial-fft
OMP_TARGET = omp-fft
ACC_GPU_TARGET = acc-gpu-fft
ACC_MC_TARGET = acc-mc-fft
CFLAGS = -lm -O2 -g -pg
CC = g++

serial:
	$(CC) $(SOURCES) $(CFLAGS) -o $(SERIAL_TARGET)

omp:
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(OMP_TARGET)

accgpu:
	pgc++ -ta=tesla:managed -fast -O3 src/main.cpp src/acc-fft.cpp -g -pg -o acc-gpu-fft

accmc:
	$(CC) $(SOURCES) $(CFLAGS) -fopenacc -o acc-mc-fft

all:
	$(CC) $(SOURCES) $(CFLAGS) -o $(SERIAL_TARGET)
	$(CC) $(SOURCES) $(CFLAGS) -fopenmp -o $(OMP_TARGET)
	pgc++ -ta=tesla:managed -fast -O3 src/main.cpp src/acc-fft.cpp -g -pg -o acc-gpu-fft

run:
	./$(SERIAL_TARGET)

clean:
	rm $(SERIAL_TARGET)
	rm $(OMP_TARGET)
	rm $(ACC_GPU_TARGET)
	rm $(ACC_MC_TARGET)
