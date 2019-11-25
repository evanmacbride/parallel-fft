# parallel-fft

Parallel implementations of Fast Fourier Transform with OpenMP and OpenACC.

Command line arguments can be given as follows:
```
./serial -waveform duration frequency
```

Where duration is a decimal number in seconds, frequency is an integer greater
than zero, and the options for test waveforms are

```-sine
-octaves
-noise
-noisy
```

```-sine``` generates a sine wave at the given test frequency. ```-octaves```
generates a sine wave at the given test frequency with two more sine waves at
one and two octaves above added. ```-noise``` generates random (unseeded) noise
(and so the frequency flag is irrelevant). ```-noisy``` combines the
```-octaves``` and ```-noise``` waveforms.

If no command line arguments are given, the defaults are ```-sine 0.5 60```.

Output is written to two .TXT files (```wave.txt``` and ```fft.txt```), both of
which can be rendered in gnuplot. ```wave.txt``` plots the test waveform, and
```fft.txt``` plots the component frequencies of the test waveform.
