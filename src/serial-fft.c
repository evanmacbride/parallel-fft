// Original serial code retrieved from:
// https://rosettacode.org/wiki/Fast_Fourier_transform#C

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

double PI;
typedef double complex cplx;

// TODO: See if we can find an iterative version of the algorithm to get better
// parallelization potential.
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
  v |= v >> 32;
  v++;
  return v;
}

int main(int argc, char *argv[])
{
	double duration = 0.5;
	int freq = 60;
	unsigned int sampling_frequency = 44100;
  unsigned long num_samples;
  // TODO: Generate square or saw wave. (Use -square or -saw flag.)
	// Fix seg fault caused by too long durations (on my machine, anything over
	// ~3s seg faults).
  //
  // Generate a test array based on command line input parameters. Make an array
  // at least as big as duration times sampling_frequency.
	if (argc == 3) {
		duration = atof(argv[2]);
	}
  if (argc > 3) {
		duration = atof(argv[2]);
		freq = atoi(argv[3]);
	}
  // Size of buf must be a power of 2 for algorithm to work.
	num_samples = next_2_power((unsigned long)(duration * sampling_frequency));
  // Initialize the buffer
  cplx buf[num_samples];
  memset(buf,0,num_samples*sizeof(cplx));

  // Fill the buffer with a num_samples number of samples spaced a step apart.
	// Adapted from
  // https://stackoverflow.com/questions/203890/creating-sine-or-square-wave-in-c-sharp
  PI = atan2(1, 1) * 4;
	double step = 0.0;
	// Write test wave input to file for gnuplot to render
	FILE *wavePlotFile = NULL;
	wavePlotFile = fopen("wave.txt","w");
	// Default to sine wave
	if (!argv[1] || strcmp(argv[1],"-sine") == 0) {
    for (int i = 0; i < num_samples; i++) {
			step += duration / num_samples;
			buf[i] = (cplx)(sin(2 * PI * step * freq));
			fprintf(wavePlotFile, "%f\t%g\n",step,creal(buf[i]));
    }
	}
	// Generate noise by using random values
	else if (strcmp(argv[1],"-noise") == 0) {
		for (int i = 0; i < num_samples; i++) {
			step += duration / num_samples;
			int r = rand();
			double randD = (double)r/INT_MAX;
			buf[i] = (cplx)(randD + (r % 3) - 1);
			fprintf(wavePlotFile, "%f\t%g\n",step,creal(buf[i]));
    }
	}
	fclose(wavePlotFile);

	// Initialize the out buffer. TODO: Allocate memory on the heap... If duration
	// is too long (more than a few seconds) and buf and out are both on the
	// stack, by the time out is initialized there will be stack overflow.
	cplx out[num_samples];
	for (int i = 0; i < num_samples; i++) {
		out[i] = buf[i];
	}
	// Run the FFT
	_fft(buf, out, num_samples, 1);

	// Plot the component frequencies represented by the FFT. Adapted from
	// Python code found here:
	// https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/
	FILE *fftPlotFile = NULL;
	fftPlotFile = fopen("fft.txt","w");
	double freq_step = 0.0;
	int end = num_samples / 2;
	double T = duration / num_samples;
	cplx val;
	for (int i = 0; i < end; i++) {
		val = (cplx)(fabs(buf[i]) / num_samples);
		fprintf(fftPlotFile, "%f\t%g\n",freq_step,val);
		freq_step += (1/T) / num_samples;
	}
	fclose(fftPlotFile);
	return 0;
}
