#!/bin/bash

TIMESTAMP=$(date +%s)
DIR=results$TIMESTAMP
mkdir $DIR
WAVEFORM=-noisy
FREQUENCY=440
DURATION_START=0.05
DURATION_STEP_FACTOR=10

echo "Compiling..."
module load pgi/19.4
make all

for (( i=0; i<5; i++ ))
do
  echo "Begin serial experiment $i"
  echo "Running serial implementation"
  > $DIR/serial-results$i.txt

  for (( j=0; j<5; j++ ))
  do
    MAG=$(($DURATION_STEP_FACTOR**$j))
    #DURATION=$(( $DURATION_START*$MAG ))
    DURATION=$(echo $DURATION_START $MAG | awk '{printf "%.6f\n",$1*$2}')
    srun ./serial-fft $WAVEFORM $DURATION $FREQUENCY >> $DIR/serial-results$i.txt
  done

  echo "Begin OpenMP experiment $i"
  echo "Testing OpenMP implementation on 2 cores"
  > $DIR/omp-2cores-results$i.txt

  for (( j=0; j<5; j++ ))
  do
    MAG=$(($DURATION_STEP_FACTOR**$j))
    DURATION=$(echo $DURATION_START $MAG | awk '{printf "%.6f\n",$1*$2}')
    srun -c 2 ./omp-fft $WAVEFORM $DURATION $FREQUENCY >> $DIR/omp-2cores-results$i.txt
  done

  echo "Testing OpenMP implementation on 4 cores"
  > $DIR/omp-4cores-results$i.txt

  for (( j=0; j<5; j++ ))
  do
    MAG=$(($DURATION_STEP_FACTOR**$j))
    DURATION=$(echo $DURATION_START $MAG | awk '{printf "%.6f\n",$1*$2}')
    srun -c 4 ./omp-fft $WAVEFORM $DURATION $FREQUENCY >> $DIR/omp-4cores-results$i.txt
  done

  echo "Testing OpenMP implementation on 8 cores"
  > $DIR/omp-8cores-results$i.txt

  for (( j=0; j<5; j++ ))
  do
    MAG=$(($DURATION_STEP_FACTOR**$j))
    DURATION=$(echo $DURATION_START $MAG | awk '{printf "%.6f\n",$1*$2}')
    srun -c 8 ./omp-fft $WAVEFORM $DURATION $FREQUENCY >> $DIR/omp-8cores-results$i.txt
  done

  echo "Testing OpenMP implementation on 16 cores"
  > $DIR/omp-16cores-results$i.txt

  for (( j=0; j<5; j++ ))
  do
    MAG=$(($DURATION_STEP_FACTOR**$j))
    DURATION=$(echo $DURATION_START $MAG | awk '{printf "%.6f\n",$1*$2}')
    srun -c 16 ./omp-fft $WAVEFORM $DURATION $FREQUENCY >> $DIR/omp-16cores-results$i.txt
  done

  echo "Begin OpenACC (GPU) experiment $i"
  > $DIR/acc-GPU-results$i.txt

  for (( j=0; j<5; j++ ))
  do
    MAG=$(($DURATION_STEP_FACTOR**$j))
    DURATION=$(echo $DURATION_START $MAG | awk '{printf "%.6f\n",$1*$2}')
    srun -p cisc372 --gres=gpu:1 ./acc-gpu-fft $WAVEFORM $DURATION $FREQUENCY >> $DIR/acc-GPU-results$i.txt
  done
done
