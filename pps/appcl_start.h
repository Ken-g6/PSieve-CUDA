/* appcu.h -- (C) Ken Brazier August 2010.

   Proth Prime Search sieve OpenCL portion (for many K and many N).
   appcl_start.h is combined with appcl.cl, using Perl, to produce appcl.h.
   You should edit appcl_start.h, not appcl.h.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
*/
#ifndef _APPCL_H
#define _APPCL_H 1

// BLOCKSIZE should be a power of two for greatest efficiency.
#define BLOCKSIZE 64
// Maximum iterations to be run by any one kernel.  Breaking up the kernel this way should improve display effects.
#define ITERATIONS_PER_KERNEL 512

#define VECSIZE 2

// Maximum time to sleep per K set, per N, in nanoseconds.
// For PPSE sieve (7/2010), 3*2,000,000 N's = about 6 seconds.
#define MAX_NS_DELAY_PER_N 3

#ifdef __cplusplus
#define KERNELS_FILENAME "appcl.cl"
#define KERNELS_FILEPATH "../../appcl.cl" // for Linux and Mac

// Kernel argument positions:
// Common.  Must always line up!:
#define ARGNO_d_P 0
#define ARGNO_d_Ps 1
#define ARGNO_d_K 2
#define ARGNO_d_factor_found 3
// check_more_ns_kernel only:
#define ARGNO_N 4

// Constants:
// start_ns_kernel
#define ARGNO_d_r0 4
#define KERNELOF_d_r0 start_ns_kernel

// check_more_ns_kernel
#define ARGNO_shift 5
#define KERNELOF_shift check_more_ns_kernel
#define ARGNO_d_kmin 6
#define KERNELOF_d_kmin check_more_ns_kernel
#define ARGNO_d_kmax 7
#define KERNELOF_d_kmax check_more_ns_kernel

// Sorry, I went a little nuts here.
// The follwing URL may help decyphering this:
// http://www.cakoose.com/wiki/c_preprocessor_abuse
#define GET_KERNELOF(SYMBOL) KERNELOF_##SYMBOL
#define RELAY_KERNELOF(SYMBOL) GET_KERNELOF(SYMBOL)
#define GET_ARGNO(SYMBOL) ARGNO_##SYMBOL
#define RELAY_ARGNO(SYMBOL) GET_ARGNO(SYMBOL)
#define STROF(S) #S
#define RELAY_STROF(S) STROF(S)
// Now, a macro to set the constants.
#define CL_MEMCPY_TO_SYMBOL(SYMBOL, SOURCE, SIZE) \
	status = clSetKernelArg(RELAY_KERNELOF(SYMBOL), RELAY_ARGNO(SYMBOL), SIZE, (void *) SOURCE); \
	if (status != CL_SUCCESS) { \
		fprintf(stderr, "Error setting %s constant %s\n", RELAY_STROF(RELAY_KERNELOF(SYMBOL)), #SYMBOL); \
		checkCUDAErr(status, "setting constant"); \
	}

// And likewise the standard arguments.
#define CL_SET_ARG(KERNEL, SYMBOL, SIZE) \
	status = clSetKernelArg(KERNEL, RELAY_ARGNO(SYMBOL), SIZE, (void *)& SYMBOL); \
	if (status != CL_SUCCESS) { \
		fprintf(stderr, "Error setting %s argument %s\n", #KERNEL, #SYMBOL); \
		checkCUDAErr(status, "setting argument"); \
	}
#define CL_SET_BUF_ARG(KERNEL, SYMBOL) CL_SET_ARG(KERNEL, SYMBOL, sizeof(cl_mem))

extern "C" {
#endif

unsigned int cuda_app_init(int gpuno, unsigned int cthread_count, int use_nvidia);
//void setup_ps(const uint64_t *P, unsigned int cthread_count);
void check_ns(const uint64_t *P, const uint64_t *Ps, const uint64_t *k0, const unsigned int cthread_count);
void get_factors_found(unsigned int *factor_found, const unsigned int cthread_count, const uint64_t start_t, int *check_ns_delay);
void cuda_finalize(void);
int get_n_subsection_start(int index);

extern unsigned int ld_nstep;
extern int ld_bbits;
extern uint64_t ld_r0;
extern unsigned int vecsize;
#ifdef __cplusplus
}
#endif

#endif
