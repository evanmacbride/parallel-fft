#ifndef COMPUTEFOURIER_H
#define COMPUTEFOURIER_H

#include "fft.h"

#include <complex.h>
//#include <map>
#include "fftw.h"
#include "filters.h"

#define OPTIMIZE_FFTW 0
//#define  WITH_COMB 0 

extern bool WITH_COMB;
extern bool ALGORITHM1;
extern bool VERBOSE;
extern bool TIMING;

//Comments located in the cc file.
  Node *outer_loop(complex_t *origx, int n, Filter *filter, int B2, int num, 
      int B, int W_Comb, int Comb_loops, int loop_threshold, int location_loops, int loops, int *hits_found, int *score, int*hits);

#endif
