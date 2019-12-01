#!/bin/bash

echo "Compiling..."

make all

for (( i=0; i<5; i++ ))
do
  echo "Begin serial experiment $i"

  ./serial-fft -noisy 0.05 440 > serial-results$i.txt
  ./serial-fft -noisy 0.5 440 >> serial-results$i.txt
  ./serial-fft -noisy 5 440 >> serial-results$i.txt
  ./serial-fft -noisy 50 440 >> serial-results$i.txt
  ./serial-fft -noisy 500 440 >> serial-results$i.txt

  echo "Begin OpenMP experiment $i"

  srun -c 2 ./omp-fft -noisy 0.05 440 > omp-2cores-results$i.txt
  srun -c 2 ./omp-fft -noisy 0.5 440 >> omp-2cores-results$i.txt
  srun -c 2 ./omp-fft -noisy 5 440 >> omp-2cores-results$i.txt
  srun -c 2 ./omp-fft -noisy 50 440 >> omp-2cores-results$i.txt
  srun -c 2 ./omp-fft -noisy 500 440 >> omp-2cores-results$i.txt

  srun -c 4 ./omp-fft -noisy 0.05 440 > omp-4cores-results$i.txt
  srun -c 4 ./omp-fft -noisy 0.5 440 >> omp-4cores-results$i.txt
  srun -c 4 ./omp-fft -noisy 5 440 >> omp-4cores-results$i.txt
  srun -c 4 ./omp-fft -noisy 50 440 >> omp-4cores-results$i.txt
  srun -c 4 ./omp-fft -noisy 500 440 >> omp-4cores-results$i.txt

  srun -c 8 ./omp-fft -noisy 0.05 440 > omp-8cores-results$i.txt
  srun -c 8 ./omp-fft -noisy 0.5 440 >> omp-8cores-results$i.txt
  srun -c 8 ./omp-fft -noisy 5 440 >> omp-8cores-results$i.txt
  srun -c 8 ./omp-fft -noisy 50 440 >> omp-8cores-results$i.txt
  srun -c 8 ./omp-fft -noisy 500 440 >> omp-8cores-results$i.txt

  srun -c 16 ./omp-fft -noisy 0.05 440 > omp-16cores-results$i.txt
  srun -c 16 ./omp-fft -noisy 0.5 440 >> omp-16cores-results$i.txt
  srun -c 16 ./omp-fft -noisy 5 440 >> omp-16cores-results$i.txt
  srun -c 16 ./omp-fft -noisy 50 440 >> omp-16cores-results$i.txt
  srun -c 16 ./omp-fft -noisy 500 440 >> omp-16cores-results$i.txt

  echo "Begin OpenACC experiment $i"

  srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 0.05 440 > acc-results$i.txt
  srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 0.5 440 >> acc-results$i.txt
  srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 5 440 >> acc-results$i.txt
  srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 50 440 >> acc-results$i.txt
  srun -p cisc372 --gres=gpu:1 ./acc-fft -noisy 500 440 >> acc-results$i.txt
done
