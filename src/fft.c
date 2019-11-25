#include <complex.h>
#include <math.h>
#include "fft.h"

double PI = atan2(1, 1) * 4;

// TODO: See if we can find an iterative version of the algorithm to get better
// parallelization potential.
void _fft(_Complex double buf[], _Complex double out[], int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);

		for (int i = 0; i < n; i += 2 * step) {
			double complex t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
