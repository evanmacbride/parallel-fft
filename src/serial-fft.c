// Retrieved from https://rosettacode.org/wiki/Fast_Fourier_transform#C

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

double PI;
typedef double complex cplx;

void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);

		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}

void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++) out[i] = buf[i];

	_fft(buf, out, n, 1);
}

// Unnecessary if we're generating data files for gnuplot
void show(const char * s, cplx buf[], int size) {
	printf("%s", s);
	for (int i = 0; i < size; i++)
		if (!cimag(buf[i]))
      // If there's no imaginary component, print 0.
			printf("(%g, 0) ", creal(buf[i]));
		else
			printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}

// Get next highest power of 2. Adapted from
// https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2
unsigned long next_2_power(unsigned long v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

int main(int argc, char *argv[])
{
  unsigned long num;
  unsigned long size;
  // TODO: Generate different test waves based on command line arguments. Accept
	// -sine, -square (or -saw), and -noise. Other command line arguments should
	// be frequency and duration. Duration should be adjusted so that, when
	// multiplied by a sampling frequency (such as the standard 44.1k), the number
	// of samples will equal the next highest power of 2 (or else the simple
	// version of the FFT algorithm in our function won't work).
  //
  // Generate a test array based on command line input parameters. Make an array
  // at least as big as a number passed after a -g flag.
  if (argc > 2) {
    if (strcmp(argv[1],"-g") == 0) {
      num = atol(argv[2]);
      // Size of buf must be a power of 2 for algorithm to work.
      size = next_2_power(num);
      // Initialize the buffer
      cplx buf[size];
      memset(buf,0,size*sizeof(cplx));

      //int sampleRate = 44100;

      // Fill the buffer with a size number of samples spaced a step apart.
			// Adapted from
      // https://stackoverflow.com/questions/203890/creating-sine-or-square-wave-in-c-sharp
      PI = atan2(1, 1) * 4;
      int freq = 60;	// A low bass note with a long enough wavelength to be be
											// easily visible. TODO: Accept command line args for freq.
			double duration = 0.25;	// A long enough time period to show several
															// cycles of freq. TODO: Accept command line
															// args for duration, then adjust so duration
															// multiplied by sampling frequency equals the
															// next highest power of 2.
			double step = 0.0;
			// Write test wave input to file for gnuplot to render
			FILE *wavePlotFile = NULL;
			wavePlotFile = fopen("wave.txt","w");
      for (int i = 0; i < size; i++) {
				step += duration / size;
				buf[i] = (cplx)(sin(2 * PI * step * freq));
				fprintf(wavePlotFile, "%f\t%g\n",step,creal(buf[i]));
      }
			fclose(wavePlotFile);

			// Run the fft
			int bufsize = sizeof(buf)/sizeof(cplx);
			fft(buf, bufsize);

			// Write fft output to file for gnuplot to render
			// TODO: Plot amplitudes against frequencies instead of the raw FFT. I
			// want to do what's done here
			// https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/
			// with Python.
			FILE *fftPlotFile = NULL;
			fftPlotFile = fopen("fft.txt","w");
			for (int i = 0; i < size; i++) {
				fprintf(fftPlotFile, "%i\t%g\n",i,creal(buf[i]));
			}
			fclose(fftPlotFile);
    }
  }
	return 0;
}
