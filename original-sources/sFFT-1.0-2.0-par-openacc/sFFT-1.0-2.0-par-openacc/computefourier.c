#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "utils.h"
#include "timer.h"
#include "computefourier.h"
#include "filters.h"
#include "fft.h"
#include "wrapper.h"

#ifdef _OPENACC
#include "cufft.h"
#include <openacc.h>
#endif

bool WITH_COMB = false;
bool ALGORITHM1   = true;
bool VERBOSE      = false;
bool TIMING       = true;

#define vprintf(...) if (VERBOSE){ printf(__VA_ARGS__);}

#if 1 
#define COMPLEX_MULT(index,off,k) \
   ax = d_origx[index];\
   ay = d_origx[index+1];\
   bx = d_filter[off+k];\
   by = d_filter[off+k+1];\
   real+=  ax * bx - ay * by;\
   imag += ax * by + ay * bx
#endif

/*use binary search to find the upper bound*/
inline int find_upper_bound(Pair *permuted_approved, int num_Comb, int val)
{ 
  int mid, low = 0, high = num_Comb - 1;
  
  if(val >= permuted_approved[high].first){
    return high;
  }
  
  mid = (low + high) / 2;
  while (high > low) {
    if (permuted_approved[mid].first >=val)
      high  = mid;
    else
      low = mid + 1;
    
    mid = (low + high) / 2;
  }
  return mid;
}

inline int timesmod(const int x, const int a, const int n) {
  return (int)((((long long int)(x))*a)&(n-1));
}


inline int Comb_Filt(complex_t *origx, int n, int num, int W_Comb, 
                      int* Comb_Approved, complex_t *x_sampt, real_t *samples, unsigned int offset, int sigma, fftw_plan p_comb)
{
  
  //Set y_i = x_[i(n/M)+sigma]
  for(int i = 0; i < W_Comb; i++){
    x_sampt[i] = origx[offset + i*sigma];
  }

  //Compute fftw of y_i
  fftw_execute_dft(p_comb, x_sampt, x_sampt); 

  for(int i = 0; i < W_Comb; i ++)
    samples[i] = cabs2(x_sampt[i]);
  
  //Find 2k largest elements of y_i in freq domain
  find_largest_indices(Comb_Approved, num, samples, W_Comb);
  
  return 0;
}

/*
  Inner loop of the algorithm, part one.

  n-dimensional origx
  permute the fourier spectrum, take the first w coordinates.
  dot with the filter
  B-dimensional FFT
  return the top num samples.

  Output to 
 */
/*
int inner_loop_locate(real_t *d_origx, int n, real_t *d_filter,
                      int B, int ai, real_t *d_x_sampt, int block_size, int filter_size){
  double DDD2 = get_time();

  int nr_blocks = B/block_size;
  int round = filter_size/B;
  unsigned index;
  int tmp, tmp2, tmp3,tmp4,tmp5,tmp6, off, off2, off3, off4;
  int n2_m_1 = n*2-1 ;
  int block_size2 = block_size*2;
  #pragma omp parallel
  {
    #pragma omp for schedule(static) private(tmp2) nowait
    for(int i=0; i<nr_blocks; i++)
    {
      tmp2 = i*block_size2;
      for(int j=0; j< block_size2; j+=2)
      {
        d_x_sampt[tmp2+j] = 0.0;
        d_x_sampt[tmp2+j+1] = 0.0;
      }
    }
    
    for(int i=0; i<round; i++)
    {
      tmp = i*2*B;
      #pragma omp for schedule(static) private(index,off,off2,off3,off4,tmp3,tmp4,tmp5,tmp6) nowait
      for(int j=0; j<nr_blocks; j+=4)
      {

        off = j*block_size2 + tmp;
        off2 = off + block_size2;
        off3 = off2 + block_size2;
        off4 = off3 + block_size2;

        tmp3 = j*block_size2;
        tmp4 = tmp3 + block_size2;
        tmp5 = tmp4 + block_size2;
        tmp6 = tmp5 + block_size2;

        for(int k=0; k<block_size2; k+=8)
        {
          index = ((unsigned)((off+ k)*ai)) & (n2_m_1);
          d_x_sampt[tmp3+k] +=   REAL(index,off,k);
          d_x_sampt[tmp3+k+1] += IMAG(index,off,k);

          index = ((unsigned)((off+ k+2)*ai)) & (n2_m_1);
          d_x_sampt[tmp3+k+2] += REAL(index, off,k+2);
          d_x_sampt[tmp3+k+3] += IMAG(index,off,k+2);
          
          index = ((unsigned)((off+ k+4)*ai)) & (n2_m_1);
          d_x_sampt[tmp3+k+4] += REAL(index, off,k+4);
          d_x_sampt[tmp3+k+5] += IMAG(index,off,k+4);

          index = ((unsigned)((off+ k+6)*ai)) & (n2_m_1);
          d_x_sampt[tmp3+k+6] += REAL(index, off,k+6);
          d_x_sampt[tmp3+k+7] += IMAG(index,off,k+6);
          

          index = ((unsigned)((off2+ k)*ai)) & (n2_m_1);
          d_x_sampt[tmp4+k] += REAL(index,off2,k);
          d_x_sampt[tmp4+k+1] +=IMAG(index,off2,k);

          index = ((unsigned)((off2+ k+2)*ai)) & (n2_m_1);
          d_x_sampt[tmp4+k+2] += REAL(index, off2,k+2);
          d_x_sampt[tmp4+k+3] += IMAG(index,off2,k+2);
        
          index = ((unsigned)((off2+ k+4)*ai)) & (n2_m_1);
          d_x_sampt[tmp4+k+4] += REAL(index, off2,k+4);
          d_x_sampt[tmp4+k+5] += IMAG(index,off2,k+4);
          
          index = ((unsigned)((off2+ k+6)*ai)) & (n2_m_1);
          d_x_sampt[tmp4+k+6] += REAL(index, off2,k+6);
          d_x_sampt[tmp4+k+7] += IMAG(index,off2,k+6);
        }
        
        for(int k=0; k<block_size2; k+=8)
        {
          index = ((unsigned)((off3+ k)*ai)) & (n2_m_1);
          d_x_sampt[tmp5+k] +=  REAL(index,off3,k);
          d_x_sampt[tmp5+k+1] +=IMAG(index,off3,k);

          index = ((unsigned)((off3+ k+2)*ai)) & (n2_m_1);
          d_x_sampt[tmp5+k+2] += REAL(index,off3,k+2);
          d_x_sampt[tmp5+k+3] += IMAG(index,off3,k+2);
        
          index = ((unsigned)((off3+ k+4)*ai)) & (n2_m_1);
          d_x_sampt[tmp5+k+4] += REAL(index,off3,k+4);
          d_x_sampt[tmp5+k+5] += IMAG(index,off3,k+4);
          
          index = ((unsigned)((off3+ k+6)*ai)) & (n2_m_1);
          d_x_sampt[tmp5+k+6] += REAL(index,off3,k+6);
          d_x_sampt[tmp5+k+7] += IMAG(index,off3,k+6);


          index = ((unsigned)((off4+ k)*ai)) & (n2_m_1);
          d_x_sampt[tmp6+k] +=  REAL(index,off4,k);
          d_x_sampt[tmp6+k+1] +=IMAG(index,off4,k);

          index = ((unsigned)((off4+ k+2)*ai)) & (n2_m_1);
          d_x_sampt[tmp6+k+2] +=REAL(index,off4,k+2);
          d_x_sampt[tmp6+k+3] +=IMAG(index,off4,k+2);
        
          index = ((unsigned)((off4+ k+4)*ai)) & (n2_m_1);
          d_x_sampt[tmp6+k+4] += REAL(index,off4,k+4);
          d_x_sampt[tmp6+k+5] += IMAG(index,off4,k+4);
          
          index = ((unsigned)((off4+ k+6)*ai)) & (n2_m_1);
          d_x_sampt[tmp6+k+6] += REAL(index,off4,k+6);
          d_x_sampt[tmp6+k+7] += IMAG(index,off4,k+6);
        }
      }
    }
  }
  
  if (TIMING) {
    double PF_T = get_time()-DDD2;
    vprintf("Step 1.A (PERM + FILTER):------------------------ %lf\n", PF_T);
  }
  return 0;
}

*/

int inner_loop_step_c(int num, int B, int *J, complex_t *x_sampt, real_t *samples)
{
  
  double DDD=get_time(), DDD2;
  
  //fftw_execute_dft(p, x_sampt, x_sampt);

  //if (TIMING) {
  //  vprintf("Step 1.B (FFTW)---------: %lf\n", get_time()-DDD);
  //  DDD2 = get_time();
  //}


  for(int i = 0; i < B; i++){
    samples[i] = cabs2(x_sampt[i]);
  }

  
  //Find the indices corresponding to the num largest element of samples in freq domain.
  //J : Output. J is sorted. J stored the indices of largest elements of smaples in freq domain.
  //num: cutoff number. We only want num largest elements
  //samples: Input
  //B: size  

  find_largest_indices(J, num, samples, B);

  if (TIMING) {
    vprintf("Step 1.C (LARGEST BUCKS): %lf\n", get_time() - DDD2);
  }

  return 0;
}

/*
 Find indices that map to J , i.e., lie within n/(2B) of (J * n/B) after permutation.

 For each such i, increment score[i] and append to hits if score[i] reaches loop_threshold.

*/ 
int inner_loop_filter_regular(int *J, int n, int num, int B, int a, int ai, int b, int loop_threshold,
                              int *score, int *hits, int *hits_found){
  double DDD = get_time();

  // Given the set of large samples, find the locations in [n] that map there
  // and output them
  int my_hits_found= *hits_found;
  //#pragma omp parallel for
  for(int i = 0; i < num; i++){
    int low = ((int)(ceil((J[i] - 0.5) * n / B)) + n)&(n-1);
    int high = ((int)(ceil((J[i] + 0.5) * n / B)) + n)&(n-1);
    int loc = timesmod(low, a, n);
    for(int j = low; j != high; j = (j + 1)&(n-1)) {
      //#pragma omp atomic
      score[loc]++;
      if(score[loc]==loop_threshold){
        //#pragma omp critical
        hits[my_hits_found++]=loc;
      }
      loc = (loc + a)&(n-1);
    }
  }

  *hits_found = my_hits_found;

  if (TIMING) {
    vprintf("Step 1.D (GROUPING):----------------------------- %lf\n\n", get_time()-DDD);
    vprintf("#####################################################################\n\n");
  }
  
  return 0;
}


/*
 Find indices that (1) map to J under the permutation and (2) lie in Comb_Approved mod W_Comb.

 For each such i, increment hits[i] and append to hits_found if hits[i] reaches loop_threshold.
*/
int inner_loop_filter_Comb(int *J, int n, int num, int B, int a, int ai, int b, int loop_threshold,
                      int *score, int *hits, int *hits_found,
                      int *Comb_Approved,int num_Comb, int W_Comb, Pair *permuted_approved){
  
  
  double DDD = get_time();
  
  
  int my_hits= *hits_found; 
  #pragma omp parallel for
  for(int m=0; m< num_Comb; m++){
    int prev = timesmod(Comb_Approved[m], ai, W_Comb);
    permuted_approved[m].first  = prev;
    permuted_approved[m].second = timesmod(prev, a, n);
  }
  
  qsort(permuted_approved, num_Comb, sizeof(permuted_approved[0]), comp_struct3);

  int nthreads = omp_get_max_threads();

  // compute intersection of permuted_approved and indices close to J * n/B, then invert to get true locations.
  if(nthreads >4)
  {
    #pragma omp parallel for
    for(int i = 0; i < num; i++){
      int low = ((int)(ceil((J[i] - 0.5) * n / B)) + n) & (n-1);
      int high = ((int)(ceil((J[i] + 0.5) * n / B)) + n) & (n-1);
      
      int key = low % W_Comb;
      
      int index = find_upper_bound(permuted_approved, num_Comb, key);
      int location = low - (low % W_Comb);
      int locinv = timesmod(location, a, n);
      for(int j = index; ; j++){
        if (j == num_Comb){
          j -= num_Comb;
          location = (location + W_Comb) & (n-1);
          //locinv = (unsigned long)(location * a) & (n-1);
          locinv = timesmod(location, a, n);
        }
        int approved_loc = location + permuted_approved[j].first;
        if((low < high && (approved_loc >= high || approved_loc < low)) ||
           (low > high && (approved_loc >= high && approved_loc < low)))
          break;
        int loc = (locinv + permuted_approved[j].second) & (n-1);
        #pragma omp atomic
        score[loc]++;
        if(score[loc]==loop_threshold){
          #pragma omp critical 
          hits[my_hits++]=loc;
        }
      }
    }
  }else
  {
    for(int i = 0; i < num; i++){
      int low = ((int)(ceil((J[i] - 0.5) * n / B)) + n) & (n-1);
      int high = ((int)(ceil((J[i] + 0.5) * n / B)) + n) & (n-1);
      
      int key = low % W_Comb;
      
      int index = find_upper_bound(permuted_approved, num_Comb, key);
      int location = low - (low % W_Comb);
      int locinv = timesmod(location, a, n);
      for(int j = index; ; j++){
        if (j == num_Comb){
          j -= num_Comb;
          location = (location + W_Comb) & (n-1);
          locinv = timesmod(location, a, n);
        }
        int approved_loc = location + permuted_approved[j].first;
        if((low < high && (approved_loc >= high || approved_loc < low)) ||
           (low > high && (approved_loc >= high && approved_loc < low)))
          break;
        int loc = (locinv + permuted_approved[j].second) & (n-1);
        score[loc]++;
        if(score[loc]==loop_threshold){
          hits[my_hits++]=loc;
        }
      }
    }
    
  }
  
  *hits_found = my_hits;
  
  if (TIMING) {
    vprintf("Step 1.D (GROUPING):----------------------------- %lf\n\n", get_time()-DDD);
    vprintf("#####################################################################\n\n");
  }

  return 0;
}

/*
  hits contains the indices that we want to estimate.

  x_samp contains a B-dimensional array for each of the `loops`
  iterations of the outer loop.  Every coordinate i of x "hashes to" a
  corresponding coordinate (permute[j] * i) mod B of x_samp[j], which
  gives an estimate of x[i].

  We estimate each coordinate as the median (independently in real and
  imaginary axes) of its images in the rows of x_samp.
 */
//Node *estimate_values(const int *hits, int *hits_found,
void estimate_values(const int *hits, int *hits_found,
                complex_t *x_samp,  const int loops, int n,
                const int *permute,
                const int B, const Filter *filter, 
                int location_loops, Node *ans)
{
  int my_hits_found = *hits_found;
  int n_div_B = n/B;
  int n_div_B_div_2 = n/B/2;
  #pragma omp parallel
  {
    real_t *values[2];
    values[0] = (real_t *) malloc(loops * sizeof(*values[0]));
    values[1] = (real_t *) malloc(loops * sizeof(*values[1]));
    #pragma omp for nowait
    for(int i = 0; i < my_hits_found; i++){
      for(int j = 0; j < loops; j++){
        int permuted_index= timesmod(permute[j], hits[i],  n);
        int hashed_to = permuted_index / n_div_B;
        int dist = permuted_index % n_div_B;
        if (dist > n_div_B_div_2) {
          hashed_to = (hashed_to + 1)%B;
          dist -= n_div_B;
        }
        dist = (n - dist) & (n-1); 
        complex_t filter_value = filter->freq[dist];
        complex_t tmp = x_samp[j*B+hashed_to] / filter_value;
        values[0][j] = creal(tmp);
        values[1][j] = cimag(tmp);
      }
    
      int location = (loops - 1) / 2;
       
      
      #ifndef WRAPPER 
      //Use qsort to find the value in location
      qsort(values[0], loops, sizeof(*values[0]), compare);
      qsort(values[1], loops, sizeof(*values[1]), compare);
      #else
      //Use wrapper to find the nth element.
      nthelement(values[0], location, loops);
      nthelement(values[1], location, loops);
      #endif
      
      real_t realv = values[0][location];
      real_t imagv = values[1][location];
      
      ans[i].key = hits[i];
      ans[i].value = realv + I*imagv;
    }
    
    free(values[0]);
    free(values[1]);
  }
}


/*
  Outer loop of the algorithm.

  If we are performing the Comb heuristic, first we do so.

  Then, `loops` times:
    choose a random permutation
    run inner_loop_locate
    if in the first location_loops loops, also run inner_loop_filter

  at the end, `hits` contains the coordinates that appear at least
  loop_threshold of location_loops times.  We estimate the values at
  these coordinates as the median of the images x_samp[loops].

  Returns a map from coordinates to estimates.
 */
Node *outer_loop(complex_t *origx, int n, Filter *filter, int B2,
  	   int num, int B, int W_Comb, int Comb_loops, int loop_threshold, int location_loops,
	   int loops, int *hits_found, int *score, int *hits){
  
  
  if(n%B){
    fprintf(stderr, "Warning: n is not divisible by B, which algorithm expects.\n");
    exit(-1);
  }
  
  //Variables used for timing
  double SCORE_T = 0;
  double PF_ALL =0;
  double G_ALL =0;
  double B_ALL = 0;
  double C_ALL =0;

  //Prepare for location+estimation loops
  int permutea[loops];
  int permute[loops];
  int permuteb[loops];
   
  real_t *restrict d_x_sampt= (real_t *restrict) malloc(B* loops * 2 * sizeof(real_t));
  real_t *samples[loops];
  int *J[loops];

  for(int i = 0; i < loops; i++){
    samples[i] = (real_t*)malloc(B* sizeof(*samples[i]));
    J[i] = (int *)malloc(num * sizeof(*J));
  } 
  
  //Prepare for Comb Filt
  complex_t *x_sampt_comb[Comb_loops];
  real_t *samples_comb[Comb_loops];
  unsigned int offset[Comb_loops]; 
  int sigma = n/W_Comb;
  int *Comb_Approved=NULL;
  if(WITH_COMB)
  {
  
    if(n%W_Comb) {
      fprintf(stderr, "Warning: W_Comb is not divisible by N, which algorithm expects.\n");
      assert(n%W_Comb == 0);
    }

    for(int i=0; i<Comb_loops; i++){
      x_sampt_comb[i] = (complex_t*) malloc(W_Comb*sizeof(*x_sampt_comb[i]));
      samples_comb[i] = (real_t*)malloc(W_Comb * sizeof(*samples_comb[i]));
      offset[i] = (unsigned) floor(drand48() * sigma);
    }
    
    Comb_Approved = (int*) malloc(Comb_loops*num*sizeof(*Comb_Approved));
  }
  
  //Create FFTW plan
  //fftw_plan p1 = fftw_plan_dft_1d(B, (complex_t*)(d_x_sampt), (complex_t *)(d_x_sampt), FFTW_FORWARD, FFTW_ESTIMATE); 
  fftw_plan p_comb = fftw_plan_dft_1d(W_Comb, x_sampt_comb[0], x_sampt_comb[0], FFTW_FORWARD, FFTW_ESTIMATE);
  *hits_found = 0;
  
  //Create cuFFT plan
  cufftHandle plan;
  cufftResult err;
  err = cufftPlan1d(&plan, B, CUFFT_Z2Z, loops); 
  if (err != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: Plan creation failed");
    exit(-1); 
  }


  reset_timer();

  //BEGIN Comb
  double DDD = get_time();
  
  int num_Comb = num;
  
  if(WITH_COMB){
    #pragma omp parallel for
    for(int i=0; i< Comb_loops; i++){
      Comb_Filt(origx, n, num, W_Comb, Comb_Approved+i*num, x_sampt_comb[i], samples_comb[i], offset[i], sigma, p_comb);
    }
    
    if(Comb_loops > 1){
      radix_sort(Comb_Approved, Comb_loops * num);
      int Last =0;
      for(int i = 1; i < Comb_loops * num; i++){
        if(Comb_Approved[i]!=Comb_Approved[Last])
           Comb_Approved[++Last]=Comb_Approved[i];
      }
      num_Comb=Last+1;
      vprintf("Comb:%d----->%d\n\n",num*Comb_loops,num_Comb);
    }
  } 
  

  if(!ALGORITHM1){
    *hits_found = num_Comb * (n/W_Comb);
    for(int j=0; j < n/W_Comb; j++)
      for (int i=0; i<num_Comb; i++)
        hits[j*num_Comb + i] = j*W_Comb + Comb_Approved[i];
  }

  double Comb_time = get_time() - DDD;
  Pair *permuted_approved = (Pair *)malloc(num_Comb * sizeof(*permuted_approved));
  //END Comb
  
  
  //BEGIN INNER LOOPS
  //Prepare for the stride factor  
  for(int i=0; i<loops; i++){
    int a=0;
    permuteb[i] = 0;
    while(gcd(a,n) != 1){
       a = (int)(random() & (n-1));
    }
    permutea[i] = a;
    permute[i] = mod_inverse(a, n);
  }
  DDD=get_time();
  
  
  //cast the complex_t* pointers to real_t* pointers and use interleaved indexing
  real_t *d_origx = (real_t *)(origx);
  real_t *d_filter = (real_t *)(filter->time); 
  int filter_size = filter->sizet;
  int round = filter_size/B;
  

#if 0 
  //simple GPU implementation 
  //int off, off2;
  //int ai;
  int B_2 = B * 2;
  int n2_m_1 = 2*n-1;
  
  #pragma acc data copyin(d_origx[0:2*n], d_filter[0:2*filter_size], permute[0:loops]) \
                   copyout(d_x_sampt[0:loops*B_2])
  {
    #pragma acc kernels loop independent
    for (int ii=0; ii<loops; ii++){
      #pragma acc loop independent
      for(int j=0; j<B_2; j+=4){
        int off = ii * B_2;
        int ai = permute[ii];
        d_x_sampt[off+j] = 0.0;
        d_x_sampt[off+j+1] = 0.0;
        d_x_sampt[off+j+2] = 0.0;
        d_x_sampt[off+j+3] = 0.0;
        for(int i=0; i<round; i+=2){
          int off2 = i * B_2;
          unsigned index = ((unsigned)((off2+j)*ai)) & n2_m_1;
          d_x_sampt[off+j] += REAL(index,off2,j);
          d_x_sampt[off+j+1] += IMAG(index,off2,j);
          
          index = ((unsigned)((off2+j+2)*ai)) & n2_m_1;
          d_x_sampt[off+j+2] += REAL(index,off2,j+2);
          d_x_sampt[off+j+3] += IMAG(index,off2,j+2);
          
          off2 = (i+1) * B_2;
          index = ((unsigned)((off2+j)*ai)) & n2_m_1;
          d_x_sampt[off+j] += REAL(index,off2,j);
          d_x_sampt[off+j+1] += IMAG(index,off2,j);
          
          index = ((unsigned)((off2+j+2)*ai)) & n2_m_1;
          d_x_sampt[off+j+2] += REAL(index,off2,j+2);
          d_x_sampt[off+j+3] += IMAG(index,off2,j+2);
        }
      }
    }
  }
#endif


#if 0
  //pre-ordered the access to filter.time 
  int B_2 = B * 2;
  int n2_m_1 = n*2-1;
  #pragma acc data copyin(d_origx[0:2*n], d_filter[0:2*filter_size], permute[0:loops]) \
                   copyout(d_x_sampt[0:loops*B_2])
  {
    #pragma acc kernels loop independent
    for (int ii=0; ii<loops; ii++){
      #pragma acc loop independent
      for(int i=0; i<B; i++){
        int i_2 = i * 2; 
        int off3 = i_2* round;
        int off = ii * B_2;
        int ai = permute[ii];

        d_x_sampt[off+i_2] = 0.0;
        d_x_sampt[off+i_2+1] = 0.0; 
        
        for(int j=0; j<round*2; j+=4){
          unsigned index = ((unsigned)((i_2+j*B)*ai)) & n2_m_1;
          d_x_sampt[off+i_2] += REAL(index,off3,j);
          d_x_sampt[off+i_2+1] += IMAG(index,off3,j);
          
          index = ((unsigned)((i_2+(j+2)*B)*ai)) & n2_m_1;
          d_x_sampt[off+i_2] += REAL(index,off3,j+2);
          d_x_sampt[off+i_2+1] += IMAG(index,off3,j+2);
        }
      }
    }
  }
#endif

#if 0 
  //pre-order the origx[]
  int off, off2;
  int ai;
  int B_2 = B * 2;
  int n2_m_1 = 2*n-1;
  
  real_t *d_origx_sort = (real_t *)acc_malloc(n*2*sizeof(real_t));

  #pragma acc data copyin(d_origx[0:2*n], d_filter[0:2*filter_size], permute[0:loops])\
                   copy(d_x_sampt[0:loops*2*B])\
                   deviceptr(d_origx_sort)
  {
    //#pragma acc kernels loop independent private(off, ai)
    for (int ii=0; ii<loops; ii++)
    {
      off = ii * B_2;
      ai = permute[ii];
      
      #pragma acc kernels loop vector(256) independent
      for(int i=0; i<filter_size*2; i+=2){
        unsigned index = ((unsigned)(i*ai)) & n2_m_1;
        d_origx_sort[i] = d_origx[index];
        d_origx_sort[i+1] = d_origx[index+1];
      }

      #pragma acc kernels loop vector(256) independent
      for(int i=0; i<B_2; i+=2){
        d_x_sampt[off+i] = 0.0;
        d_x_sampt[off+i+1] = 0.0; 
      }
      
      #pragma acc kernels loop vector(256) independent
      for(int j=0; j<B_2; j+=4){
        for(int i=0; i<round; i++){
          off2 = i * B_2;
          //unsigned index = ((unsigned)((off2+j)*ai)) & n2_m_1;
          d_x_sampt[off+j] += REAL(off2,j);
          d_x_sampt[off+j+1] += IMAG(off2,j);
          
          //index = ((unsigned)((off2+j+2)*ai)) & n2_m_1;
          d_x_sampt[off+j+2] += REAL(off2,j+2);
          d_x_sampt[off+j+3] += IMAG(off2,j+2);
        }
      }
    }
  }  
  acc_free(d_origx_sort);
#endif

#if 0 
  //Blocking + Unrolling Optimizations
  int nr_blocks = B/block_size;
  unsigned index;
  int tmp, tmp2, tmp3,tmp4,tmp5,tmp6, off, off2, off3, off4, off5;
  int n2_m_1 = n*2-1 ;
  int block_size2 = block_size*2;
  int ai;
  #pragma acc data copyin(d_origx[0:2*n], d_filter[0:2*filter_size], permute[0:loops])\
                     copy(d_x_sampt[0:loops*2*B])
  {
    for(int ii=0; ii<loops; ii++)
    {
      ai = permute[ii];
      off5 = ii * 2 * B;
      #pragma acc kernels loop independent
       for(int i=0; i<nr_blocks; i++){
         for(int j=0; j< block_size2; j+=2){
            d_x_sampt[off5+i*block_size2+j] = 0.0;
            d_x_sampt[off5+i*block_size2+j+1] = 0.0;
         }
       }

      for(int i=0; i<round; i++)
      {
        tmp = i*2*B;
        #pragma acc kernels loop vector(256) private(off,off2,off3,off4,tmp3,tmp4,tmp5,tmp6) independent
        for(int j=0; j<nr_blocks; j+=2)
        {
          off = j*block_size2 + tmp;
          off2 = off + block_size2;
          off3 = off2 + block_size2;
          off4 = off3 + block_size2;

          tmp3 = j*block_size2;
          tmp4 = tmp3 + block_size2;
          tmp5 = tmp4 + block_size2;
          tmp6 = tmp5 + block_size2;

          #pragma acc loop private(index) vector(256) independent
          for(int k=0; k<block_size2; k+=4)
          {
            index = ((unsigned)((off+ k)*ai)) & (n2_m_1);
            d_x_sampt[off5+tmp3+k] +=   REAL(index,off,k);
            d_x_sampt[off5+tmp3+k+1] += IMAG(index,off,k);

            index = ((unsigned)((off+ k+2)*ai)) & (n2_m_1);
            d_x_sampt[off5+tmp3+k+2] += REAL(index, off,k+2);
            d_x_sampt[off5+tmp3+k+3] += IMAG(index,off,k+2);
            

            index = ((unsigned)((off2+ k)*ai)) & (n2_m_1);
            d_x_sampt[off5+tmp4+k] += REAL(index,off2,k);
            d_x_sampt[off5+tmp4+k+1] +=IMAG(index,off2,k);

            index = ((unsigned)((off2+ k+2)*ai)) & (n2_m_1);
            d_x_sampt[off5+tmp4+k+2] += REAL(index, off2,k+2);
            d_x_sampt[off5+tmp4+k+3] += IMAG(index,off2,k+2);
          }
        }
      }
    }
  }
#endif

#if 1
  //pre-ordered the access to filter.time 
  int B_2 = B * 2;
  int n2_m_1 = n*2-1;
  int round_2 = round * 2;
  #pragma acc data copyin(d_origx[0:2*n], d_filter[0:2*filter_size], permute[0:loops]) \
                   copyout(d_x_sampt[0:loops*B_2])
  {
    #pragma acc kernels loop gang vector(8) independent
    for (int ii=0; ii<loops; ii++){
      #pragma acc loop gang vector(64) independent
      for(int i=0; i<B; i++){
        int i_2 = i * 2; 
        int off3 = i_2* round;
        int off = ii * B_2;
        int ai = permute[ii];
        real_t real=0.0;
        real_t imag = 0.0;
        real_t ax, ay, bx, by;
        unsigned tmp, index;
        for(int j=0; j<round_2; j+=4){
          tmp = ((unsigned)((i_2+j*B)*ai)); 
          index = tmp & n2_m_1;
          COMPLEX_MULT(index,off3,j);
          
          index = (unsigned)(tmp + B*2*ai) & n2_m_1;
          COMPLEX_MULT(index,off3,j+2);
        }
        d_x_sampt[off+i_2] =real;
        d_x_sampt[off+i_2+1] =imag; 
      }
    }
    PF_ALL = get_time() - DDD;
    
    
    DDD = get_time();
    //Step B -- cuFFT of B-dimensional FFT
    #pragma acc host_data use_device(d_x_sampt)
    {
      err = cufftExecZ2Z(plan, (cufftDoubleComplex *)d_x_sampt, (cufftDoubleComplex *)d_x_sampt, CUFFT_FORWARD);
      if (err != CUFFT_SUCCESS){
        fprintf(stderr, "CUFFT error: Execution failed, error code is %d\n", err);
        exit(-1); 
      }
    }
  }/*End of ACC data region*/
#endif
    B_ALL = get_time() - DDD;

  //cast the real_t* back to complex_t*
  complex_t *x_sampt = (complex_t *)(d_x_sampt);
  
  //Step C -- Find largest indices
  DDD=get_time();
  #pragma omp parallel for
  for(int i=0; i< loops; i++)
  {
      inner_loop_step_c(num, B, J[i], x_sampt+i*B,samples[i]);
  } 
  C_ALL = get_time() - DDD;

  //Step D -- Grouping 
  DDD=get_time();
  if(!WITH_COMB)
    for (int i=0; i<location_loops; i++){
      inner_loop_filter_regular(J[i], n, num, B,
                                permutea[i], permute[i], permuteb[i], loop_threshold,
                                score, hits, hits_found);
  }else
    for (int i=0; i<location_loops; i++){
      inner_loop_filter_Comb(J[i], n, num, B,
                           permutea[i], permute[i], permuteb[i], loop_threshold,
                           score, hits, hits_found,
                           Comb_Approved, num_Comb, W_Comb, permuted_approved);
    }
  
  G_ALL=get_time() - DDD;
  //END INNER LOOPS
  
  //BEGIN ESTIMATATION LOOPS
  DDD = get_time();
  Node *ans = (Node *) malloc(*hits_found*sizeof(*ans));

  estimate_values(hits, hits_found, x_sampt,  loops, n, permute, B, filter, location_loops, ans);

  //END ESTIMATION
  double E_T = get_time() - DDD;

  DDD = get_time();
  
  for(int i = 0; i < loops; i++)
  {
    free(samples[i]);
    free(J[i]);
  }
  free(d_x_sampt);
  if(WITH_COMB)
  {
    for (int i=0; i<Comb_loops; i++)
    {
      free(x_sampt_comb[i]);
      free(samples_comb[i]);
    }
    free(Comb_Approved);
    fftw_destroy_plan(p_comb);
  }
  
  free(permuted_approved);
  
  //destroy plan
  cufftDestroy(plan);


  if(TIMING){
	  printf("Total sFFT time: %lf\n", DDD); 
	  printf("Time distribution: scoretable  Comb     perm+filter grouping estimation  stepB    stepC    other    total\n");
	  printf("                     %lf %lf    %lf %lf   %lf %lf %lf  %lf %lf\n",
			 SCORE_T, Comb_time, PF_ALL, G_ALL, E_T, B_ALL, C_ALL, DDD-PF_ALL-G_ALL-E_T-Comb_time-B_ALL-C_ALL-SCORE_T, DDD);
	  double tott = (DDD)/100;
	  printf("                        %4.1lf%%    %4.1lf%%       %4.1lf%%    %4.1lf%%      %4.1lf%%    %4.1lf%%  %4.1lf%%     %4.1lf%%   %5.1lf%%\n",
			 SCORE_T / tott, Comb_time/tott, PF_ALL/tott, G_ALL/tott, E_T/tott, B_ALL/tott, C_ALL/tott, (DDD-PF_ALL-G_ALL-E_T-Comb_time-B_ALL-C_ALL-SCORE_T)/tott, (DDD)/tott);

	  printf("\n");
  }

  return ans;

}  
