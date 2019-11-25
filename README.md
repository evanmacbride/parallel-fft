# parallel-fft

## Parallel implementations of Fast Fourier Transform with OpenMP and OpenACC

This project uses the [OpenMP](https://www.openmp.org/) and
[OpenACC](https://www.openacc.org/) parallel programming frameworks to implement
the
[Cooley-Tukey](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)
fast Fourier transform (FFT) algorithm.

### How To Use

A makefile is provided. The options are
```
make serial
make omp
make acc
make all
```

Command line arguments can be given as follows:
```
./serial -waveform duration frequency
```

where ```duration``` is a decimal number in seconds, ```frequency``` is an
integer greater than zero, and the options for ```waveform``` are

```
-sine
-octaves
-noise
-noisy
```

```-sine``` generates a sine wave at the given test frequency. ```-octaves```
generates a sine wave at the given test frequency with two more sine waves at
one and two octaves added above. ```-noise``` generates random (unseeded) noise
(and so in this case the frequency flag is irrelevant). ```-noisy``` combines
the ```-octaves``` and ```-noise``` waveforms.

If no command line arguments are given, the defaults are ```-sine 0.5 60```.

Output is written to two .TXT files (```wave.txt``` and ```fft.txt```), both of
which can be rendered in [gnuplot](http://www.gnuplot.info/). ```wave.txt```
plots the test waveform, and ```fft.txt``` plots the component frequencies of
the test waveform.

### Sources

Several functions were adapted from the following sources:

[Original serial code snippet](https://rosettacode.org/wiki/Fast_Fourier_transform#C)
[Rounding to the next power of 2](https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2)
[Modeling a sine wave](https://stackoverflow.com/questions/203890/creating-sine-or-square-wave-in-c-sharp)
[Plotting component frequencies](https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/)
