# parallel-fft

## Parallel implementations of Fast Fourier Transform with OpenMP and OpenACC

This project uses the [OpenMP](https://www.openmp.org/) and
[OpenACC](https://www.openacc.org/) parallel programming frameworks to implement
the
[Cooley-Tukey](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)
fast Fourier transform (FFT) algorithm. It was created by Evan MacBride, Seth
Prentice, and Ryan Barber for the University of Delaware's CISC372 Parallel
Computing course.

### How To Use

#### Compiling

A makefile is provided. The options are
```
make serial
make omp
make accgpu
make all
```

```make serial``` compiles the original, unparallelized code, ```make omp``` compiles an OpenMP parallelization, ```make accgpu``` compiles an OpenACC (GPU) parallelization, and ```make all``` compiles all three.

There are also commands for ```make clean``` to remove all executables and
```make run``` to run ```serial-fft``` with default settings.

#### Running the program

Command line arguments can be given as follows:
```
./executable -waveform duration frequency
```

where ```duration``` is a decimal number in seconds, ```frequency``` is an
integer greater than zero, the options for ```waveform``` are

```
-sine
-octaves
-noise
-noisy
```

and the options for ```executable``` are

```
serial-fft
omp-fft
acc-gpu-fft
```

```serial-fft``` is the original, unparallelized code, ```omp-fft``` is an
OpenMP parallelization, and ```acc-gpu-fft``` is an OpenACC (GPU) parallelization.


```-sine``` generates a sine wave at the given test frequency. ```-octaves```
generates a sine wave at the given test frequency with two more sine waves at
one and two octaves added above. ```-noise``` generates random (unseeded) noise
(and so in this case the frequency flag is irrelevant). ```-noisy``` combines
the ```-octaves``` and ```-noise``` waveforms.

If no command line arguments are given, the defaults are ```-sine 0.5 60```.

Input duration and runtime are printed as output. Additionally, two .TXT files
(```wave.txt``` and ```fft.txt```) are generated, both of
which can be rendered in [gnuplot](http://www.gnuplot.info/). ```wave.txt```
plots the test waveform, and ```fft.txt``` plots the component frequencies of
the test waveform.

### Experiments

```experiments.sh``` is included to help run the program and organize output
data into a format that can be plotted in [gnuplot](http://www.gnuplot.info/).
Input duration and runtime are given as output, which is saved in .txt
files and organized into timestamped folders. ```wave.txt``` and ```fft.txt```
are generated with each iteration of each implementation, and all but the last
instances of these files will be overwritten. Input durations and runtimes for
individual iterations of the different implementations will be saved into their
own appropriately named .txt files.

### Sources

Some code was adapted from the following sources:
* [Recursive serial code snippet](https://rosettacode.org/wiki/Fast_Fourier_transform#C)
* [Iterative serial code](https://www.nayuki.io/page/free-small-fft-in-multiple-languages)
* [Rounding to the next power of 2](https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2)
* [Modeling a sine wave](https://stackoverflow.com/questions/203890/creating-sine-or-square-wave-in-c-sharp)
* [Plotting component frequencies](https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/)
* [Timer code](https://www.appentra.com/download/8899/)
* [Using managed memory with C++ vectors and OpenACC](https://stackoverflow.com/a/54027319)
