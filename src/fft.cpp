// Copyright notice from original iterative serial code follows:
/*
 * Free FFT and convolution (C++)
 *
 * Copyright (c) 2019 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <omp.h>
#include "fft.hpp"

using std::complex;
using std::size_t;
using std::vector;

// Private function prototypes
static size_t reverseBits(size_t x, int n);

void Fft::transformRadix2(vector<complex<double> > &vec) {
	// Length variables
	size_t n = vec.size();
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if (static_cast<size_t>(1U) << levels != n)
		throw std::domain_error("Length is not a power of 2");

	// Trignometric table
	vector<complex<double> > expTable(n / 2);
	size_t i;
#pragma omp parallel for private(i) shared(expTable,n)
	for (i = 0; i < n / 2; i++)
		expTable[i] = std::polar(1.0, -2 * M_PI * i / n);

	// Bit-reversed addressing permutation
#pragma omp parallel for private(i) shared(vec,n)
	for (i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i)
			std::swap(vec[i], vec[j]);
	}

	size_t halfsize;
	size_t tablestep;
	complex<double> temp;
	size_t size, j, k;
	// Cooley-Tukey decimation-in-time radix-2 FFT
#pragma omp parallel private(size,halfsize,tablestep,i,j,k,temp) shared(vec)
{
	for (size = 2; size <= n; size *= 2) {
		halfsize = size / 2;
		tablestep = n / size;
#pragma omp for
		for (i = 0; i < n; i += size) {
			for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				temp = vec[j + halfsize] * expTable[k];
				vec[j + halfsize] = vec[j] - temp;
				vec[j] += temp;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}
}

static size_t reverseBits(size_t x, int n) {
	size_t result = 0;
	for (int i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1U);
	return result;
}
