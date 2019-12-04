#!/bin/bash

TIMESTAMP=$(date +%s)
DIR=results$TIMESTAMP
mkdir $DIR
WAVEFORM=-noisy
FREQUENCY=440
DURATION_START=0.05
DURATION_STEP_FACTOR=10

echo "Compiling..."

make all

for (( i=0; i<5; i++ ))
do
  echo "Begin serial experiment $i"
  echo "Running serial implementation"

  for (( j=0; j<5; j++ ))
  do
    srun ./serial-fft $WAVEFORM (( $DURATION_START*$DURATION_STEP_FACTOR**$j )) $FREQUENCY > $DIR/serial-results$i.txt
  done

  echo "Begin OpenMP experiment $i"
  echo "Testing OpenMP implementation on 2 cores"

  for (( j=0; j<5; j++ ))
  do
    srun -c 2 ./omp-fft $WAVEFORM (( $DURATION_START*$DURATION_STEP_FACTOR**$j )) $FREQUENCY > $DIR/omp-2cores-results$i.txt
  done

  echo "Testing OpenMP implementation on 4 cores"

  for (( j=0; j<5; j++ ))
  do
    srun -c 4 ./omp-fft $WAVEFORM (( $DURATION_START*$DURATION_STEP_FACTOR**$j )) $FREQUENCY > $DIR/omp-4cores-results$i.txt
  done

  echo "Testing OpenMP implementation on 8 cores"

  for (( j=0; j<5; j++ ))
  do
    srun -c 8 ./omp-fft $WAVEFORM (( $DURATION_START*$DURATION_STEP_FACTOR**$j )) $FREQUENCY > $DIR/omp-8cores-results$i.txt
  done

  echo "Testing OpenMP implementation on 16 cores"

  for (( j=0; j<5; j++ ))
  do
    srun -c 16 ./omp-fft $WAVEFORM (( $DURATION_START*$DURATION_STEP_FACTOR**$j )) $FREQUENCY > $DIR/omp-16cores-results$i.txt
  done

  echo "Begin OpenACC (GPU) experiment $i"

  for (( j=0; j<5; j++ ))
  do
    srun -p cisc372 --gres=gpu:1 ./acc-fft $WAVEFORM (( $DURATION_START*$DURATION_STEP_FACTOR**$j )) $FREQUENCY > $DIR/acc-GPU-results$i.txt
  done
done
