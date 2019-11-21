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

// Fill an array with generated values.
//cplx* generate(int size) {
//  cplx gen[size];
//  for (int i = 0; i < size; i++) {
//    gen[i] = 1;
//  }
//  printf("%g\n", creal(gen[0]));
//  return gen;
//}

int main(int argc, char *argv[])
{

  //FILE *fptr;
  unsigned long num;
  unsigned long size;
  // TODO: Get filename from command line arguments. Check that argv[1] equals
  // some flag, such as -f .
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

			//int amp = 0.25 * INT_MAX;
      //int sampleRate = 44100;

      // Fill the buffer with a size number of samples spaced a step apart.
			// Adapted from
      // https://stackoverflow.com/questions/203890/creating-sine-or-square-wave-in-c-sharp
      PI = atan2(1, 1) * 4;
      int freq = 60;	// A low bass note with a long enough wavelength to be be
											// easily visible
			double timeInterval = 0.25;	// A long enough time period to show several
																	// cycles of freq
			double step = 0.0;
			// Write test wave input to file for gnuplot to render
			FILE *wavePlotFile = NULL;
			wavePlotFile = fopen("wave.txt","w");
      for (int i = 0; i < size; i++) {
				step += timeInterval / size;
				buf[i] = (cplx)(sin(2 * PI * step * freq));
				fprintf(wavePlotFile, "%f\t%g\n",step,creal(buf[i]));
        //buf[i] = (cplx)(amp * sin((2 * PI * i * freq) / sampleRate));
				//fprintf(wavePlotFile, "%i\t%g\n",i,creal(buf[i]));
      }
			fclose(wavePlotFile);
      //cplx buf[] = {1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0};
      //int bufsize = sizeof(buf)/sizeof(cplx);
    	//show("Data: ", buf, bufsize);
    	//fft(buf, bufsize);
    	//show("\nFFT : ", buf, bufsize);
      //printf("\n");

			int bufsize = sizeof(buf)/sizeof(cplx);
			fft(buf, bufsize);

			// Write fft output to file for gnuplot to render
			FILE *fftPlotFile = NULL;
			fftPlotFile = fopen("fft.txt","w");
			for (int i = 0; i < size; i++) {
				fprintf(fftPlotFile, "%i\t%g\n",i,creal(buf[i]));
			}
			fclose(fftPlotFile);
    }
  //  fptr = fopen(argv[2],"r");
  }

  //if (fptr != NULL) {
  //  fclose(fptr);
  //}
	return 0;
}
