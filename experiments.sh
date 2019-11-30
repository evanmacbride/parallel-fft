#!/bin/bash

echo "Compiling..."

make all

echo "Begin serial experiment"

./serial-fft -noisy 0.05 440 > serial-results.txt
./serial-fft -noisy 0.5 440 >> serial-results.txt
./serial-fft -noisy 5 440 >> serial-results.txt
./serial-fft -noisy 50 440 >> serial-results.txt
./serial-fft -noisy 500 440 >> serial-results.txt

echo "Begin OpenMP experiment"

srun -c 4 ./omp-fft -noisy 0.05 440 > omp-4cores-results.txt
srun -c 4 ./omp-fft -noisy 0.5 440 >> omp-4cores-results.txt
srun -c 4 ./omp-fft -noisy 5 440 >> omp-4cores-results.txt
srun -c 4 ./omp-fft -noisy 50 440 >> omp-4cores-results.txt
srun -c 4 ./omp-fft -noisy 500 440 >> omp-4cores-results.txt

srun -c 8 ./omp-fft -noisy 0.05 440 > omp-8cores-results.txt
srun -c 8 ./omp-fft -noisy 0.5 440 >> omp-8cores-results.txt
srun -c 8 ./omp-fft -noisy 5 440 >> omp-8cores-results.txt
srun -c 8 ./omp-fft -noisy 50 440 >> omp-8cores-results.txt
srun -c 8 ./omp-fft -noisy 500 440 >> omp-8cores-results.txt

srun -c 12 ./omp-fft -noisy 0.05 440 > omp-12cores-results.txt
srun -c 12 ./omp-fft -noisy 0.5 440 >> omp-12cores-results.txt
srun -c 12 ./omp-fft -noisy 5 440 >> omp-12cores-results.txt
srun -c 12 ./omp-fft -noisy 50 440 >> omp-12cores-results.txt
srun -c 12 ./omp-fft -noisy 500 440 >> omp-12cores-results.txt

srun -c 16 ./omp-fft -noisy 0.05 440 > omp-16cores-results.txt
srun -c 16 ./omp-fft -noisy 0.5 440 >> omp-16cores-results.txt
srun -c 16 ./omp-fft -noisy 5 440 >> omp-16cores-results.txt
srun -c 16 ./omp-fft -noisy 50 440 >> omp-16cores-results.txt
srun -c 16 ./omp-fft -noisy 500 440 >> omp-16cores-results.txt

echo "Begin OpenACC experiment"

srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 0.05 440 > acc-results.txt
srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 0.5 440 >> acc-results.txt
srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 5 440 >> acc-results.txt
srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 50 440 >> acc-results.txt
srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 500 440 >> acc-results.txt
