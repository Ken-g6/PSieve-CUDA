/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.cu -- (C) Ken Brazier February 2010.

   Proth Prime Search sieve CUDA portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "main.h"
#include "util.h"
#include "app.h"
#include "appcu.h"
#define INLINE static inline
#define BITSATATIME 4
#define BITSMASK ((1<<BITSATATIME)-1)
// BLOCKSIZE should be a power of two for greatest efficiency.
#define BLOCKSIZE 128
#if(BITSATATIME == 3)
  #define SHIFT_CAST unsigned int
#elif(BITSATATIME == 4)
  #define SHIFT_CAST uint64_t
#else
  #error "Invalid BITSATATIME."
#endif
// Extern vars in appcu.h:
unsigned int ld_nstep;

// Device constants
__constant__ unsigned int d_bitsatatime;
__constant__ unsigned int d_len;//=(1<<bitsatatime); 
__constant__ unsigned int d_halflen;//=(1<<bitsatatime)/2; 
__constant__ uint64_t d_kmax;
__constant__ unsigned int d_bpernstep;
__constant__ unsigned int d_nmin;
__constant__ unsigned int d_nmax;
__constant__ unsigned int d_nstep;
__constant__ unsigned int d_search_proth;
// Device arrays
uint64_t *d_P;
uint64_t *d_K;
uint64_t *d_bitsskip;
unsigned char *d_factor_found;

/* This function is called once before any threads are started.
 */
unsigned int cuda_app_init(int gpuno)
{
  unsigned int i;
  struct cudaDeviceProp gpuprop;
  unsigned int ld_bitsatatime = 0;
  //unsigned int ld_halflen=(1<<bitsatatime)/2; 
  unsigned int ld_bitsmask;
  unsigned int ld_bpernstep;
  unsigned int cthread_count;

  // Find the GPU's properties.
  if(cudaGetDeviceProperties(&gpuprop, gpuno) != cudaSuccess) {
    fprintf(stderr, "GPU %d not compute-capable.\n", gpuno);
    return 0;
  }
  cudaSetDevice(gpuno);
  fprintf(stderr, "Detected GPU %d: %s\n", gpuno, gpuprop.name);
  fprintf(stderr, "Detected compute capability: %d.%d\n", gpuprop.major, gpuprop.minor);
  fprintf(stderr, "Detected %d multiprocessors.\n", gpuprop.multiProcessorCount);
  fprintf(stderr, "Detected %lu bytes of device memory.\n", gpuprop.totalGlobalMem);

  // Use them to set cthread_count.
  // First, threads per multiprocessor, based on compute capability.
  cthread_count = (gpuprop.major == 1 && gpuprop.minor < 2)?384:768;
  cthread_count *= gpuprop.multiProcessorCount;

  if(gpuprop.totalGlobalMem < cthread_count*48) {
    fprintf(stderr, "Insufficient GPU memory: %u bytes.\n",  (unsigned int)(gpuprop.totalGlobalMem));
    return 0;
  }
  // Calculate ld_bitsatatime given memory constraints, and possibly nmin-nmax via nstep vs. 2^ld_bitsatatime
  // Things change if nmax-nmin < 1000000 or so, but for now let's go with a constant maximum of ld_bitsatatime<=13.
  i = gpuprop.totalGlobalMem/sizeof(uint64_t); // Total number of 64-bit numbers that can be stored.
#ifndef NDEBUG
  fprintf(stderr, "Available memory = %d bytes\n", (int)(gpuprop.totalGlobalMem));
#endif
  ld_bitsatatime = BITSATATIME;
  ld_bitsmask = BITSMASK+1;

  // Allocate device arrays:
  // TODO: fix this awkward construct.
  while(1) {
    // - d_bitsskip[] (Biggest array first.)
    // Not using cudaMallocPitch because coalescing isn't possible in general.
    if(cudaMalloc((void**)&d_bitsskip, ld_bitsmask*cthread_count*sizeof(uint64_t)) == cudaSuccess) {
      // - P's
      if(cudaMalloc((void**)&d_P, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
        // - K's
        if(cudaMalloc((void**)&d_K, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
          // - d_factor_found[]
          if(cudaMalloc((void**)&d_factor_found, cthread_count*sizeof(unsigned char)) == cudaSuccess) {
#ifndef NDEBUG
            fprintf(stderr, "Allocation successful!\n");
            fprintf(stderr, "ld_bitsatatime = %u\n", ld_bitsatatime);
#endif
            break;  // Allocation successful!
          }
          cudaFree(d_K);
        }
        cudaFree(d_P);
      }
      cudaFree(d_bitsskip);
    }
    fprintf(stderr, "Insufficient available memory on GPU %d.\n", gpuno);
    return 0;
  }

  ld_bitsmask--; // Finalize bitsmask

  // Calculate the values that fit the given ld_bitsatatime.
  // ld_nstep was previously calculated in app_init.
  if (ld_nstep > ld_bitsatatime) {
    ld_bpernstep = ld_nstep/ld_bitsatatime;
    ld_nstep = ld_bpernstep*ld_bitsatatime;
  }
  if (ld_nstep > (nmax-nmin+1))
    ld_nstep = (nmax-nmin+1);

  // Set the constants.
  cudaMemcpyToSymbol(d_bitsatatime, &ld_bitsatatime, sizeof(ld_bitsatatime));

  // The following would be "pitch" if using cudaMallocPitch above.
  i = ld_bitsmask+1;
  cudaMemcpyToSymbol(d_len, &i, sizeof(i));//=(1<<bitsatatime); 
  // But the following would not.
  i >>= 1;
  cudaMemcpyToSymbol(d_halflen, &i, sizeof(i));//=(1<<bitsatatime)/2; 
  cudaMemcpyToSymbol(d_kmax, &kmax, sizeof(kmax));
  cudaMemcpyToSymbol(d_bpernstep, &ld_bpernstep, sizeof(ld_bpernstep));
  cudaMemcpyToSymbol(d_nmin, &nmin, sizeof(nmin));
  cudaMemcpyToSymbol(d_nmax, &nmax, sizeof(nmax));
  cudaMemcpyToSymbol(d_nstep, &ld_nstep, sizeof(ld_nstep));
  cudaMemcpyToSymbol(d_search_proth, &search_proth, sizeof(search_proth));

  return cthread_count;
}



// Set up the lookup tables for all P's.
// TODO: Fix the variable names here to be more descriptive.
__global__ void d_setup_ps(const uint64_t *P, uint64_t *bitsskip) {
  unsigned int n, i;
  uint64_t *bs0 = &bitsskip[blockIdx.x * BLOCKSIZE*d_len + threadIdx.x];
  uint64_t my_P;
  uint64_t kpos;
  SHIFT_CAST mul_shift = 0, off_shift = 0;
  unsigned int my_factor_found = 0;
  // Initialize bitsskip array.
  my_P = P[blockIdx.x * BLOCKSIZE + threadIdx.x];

  // Initialize the first two entries.
  bs0[BLOCKSIZE*d_halflen] = (my_P+1)/2;	// Needed first.
  // bs0[0] will be ignored; it's just 0.

  // Fill in the intervening spaces, two numbers at a time.
  for(i=d_halflen; i > 1; i >>= 1) {
    for(n=i/2; n < d_halflen; n+=i) {
      kpos = bs0[BLOCKSIZE*2*n];
      my_factor_found = ((unsigned int)kpos)&1;
      //printf("Filling n=%d from bs0=%lu\n", n, kpos);
      bs0[BLOCKSIZE*n] = (kpos+((my_factor_found)?my_P:(uint64_t)0))/2;
      //printf("Filling n=%d\n", n+d_halflen);
      bs0[BLOCKSIZE*(n+d_halflen)] = (kpos+1+((my_factor_found)?(uint64_t)0:my_P))/2;
    }
  }
  // Now convert the entries to multiples of 1/8*P, plus some small constant.
  my_P >>= BITSATATIME;
  // Might as well skip 0, which is always 0.
  for(i=1; i < d_len; i++) {
    kpos = bs0[i*BLOCKSIZE];
    n = (unsigned int)(kpos/my_P);
    mul_shift |= ((SHIFT_CAST)n) << (BITSATATIME*i);
    off_shift |= (kpos-(((SHIFT_CAST)n)*my_P)) << (BITSATATIME*i);
  }

#if(BITSATATIME == 3)
  bs0[0] = (((uint64_t)off_shift)<<32) + mul_shift;
#elif(BITSATATIME == 4)
  bs0[0] = mul_shift;
  bs0[BLOCKSIZE] = off_shift;
#else
  #error "Invalid BITSATATIME."
#endif
}

// Check all N's.
__global__ void d_check_ns(const uint64_t *P, const uint64_t *K, unsigned char *factor_found_arr, uint64_t *bitsskip) {
  unsigned int n = d_nmin; // = nmin;
  unsigned int i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  uint64_t k0;
  uint64_t kpos;
  unsigned char my_factor_found = 0;
  uint64_t my_P;
  unsigned int shift;
#ifndef NDEBUG
  uint64_t *bs0 = &bitsskip[blockIdx.x * BLOCKSIZE*d_len + threadIdx.x];
#endif
  SHIFT_CAST mul_shift, off_shift;

  //factor_found_arr[i] = 0;
  i = blockIdx.x * BLOCKSIZE*d_len + threadIdx.x;
#if(BITSATATIME == 3)
  my_P = bitsskip[i];
  mul_shift = (unsigned int)my_P;
  off_shift = (unsigned int)(my_P >> 32);
#elif(BITSATATIME == 4)
  mul_shift = bitsskip[i];
  off_shift = bitsskip[i+BLOCKSIZE];
#else
  #error "Invalid BITSATATIME."
#endif
  i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  k0 = K[i];
  my_P = P[i];
  
  if(d_search_proth) k0 = my_P-k0;
  my_P >>= BITSATATIME;
  my_factor_found = 0;
  do { // Remaining steps are all of equal size nstep
    kpos = k0;
    i = __ffsll(kpos)-1;

    kpos >>= i;
    if (kpos <= d_kmax) {
#ifndef NDEBUG
      fprintf(stderr, "%u | %u*2^%u+1 (P[%d])\n", (unsigned int)my_P, (unsigned int)kpos, n+i, blockIdx.x * BLOCKSIZE + threadIdx.x);
#endif
      // Just flag this if kpos <= d_kmax.
      my_factor_found = 1;
    }

    for(i=0; i < d_bpernstep; i++) {
      shift=BITSATATIME*(((unsigned int)k0)&BITSMASK);
#ifndef NDEBUG
      if(shift > BITSATATIME && (bs0[((unsigned int)k0 & BITSMASK)*BLOCKSIZE] != ((((unsigned int)(mul_shift>>shift))&BITSMASK)*my_P + (((unsigned int)(off_shift>>shift))&BITSMASK)))) {
        fprintf(stderr, "Array lookup[%d], %lu != register lookup %lu\n", (unsigned int)k0 & BITSMASK, bs0[((unsigned int)k0 & BITSMASK)*BLOCKSIZE], ((((unsigned int)(mul_shift>>shift))&BITSMASK)*my_P + (((unsigned int)(off_shift>>shift))&BITSMASK)));
      }
      assert((shift == 0 && ((((unsigned int)(mul_shift>>shift))&BITSMASK)*my_P + (((unsigned int)(off_shift>>shift))&BITSMASK) == 0)) || shift == BITSATATIME || (bs0[((unsigned int)k0 & BITSMASK)*BLOCKSIZE] == ((((unsigned int)(mul_shift>>shift))&BITSMASK)*my_P + (((unsigned int)(off_shift>>shift))&BITSMASK))));
#endif
      k0 = (k0 >> BITSATATIME) + (((unsigned int)(mul_shift>>shift))&BITSMASK)*my_P + (((unsigned int)(off_shift>>shift))&BITSMASK);
    }
    n += d_nstep;
  } while (n < d_nmax);
  factor_found_arr[blockIdx.x * BLOCKSIZE + threadIdx.x] = my_factor_found;
}

// Pass the arguments to the CUDA device, run the code, and get the results.
void check_ns(const uint64_t *P, const uint64_t *K, unsigned char *factor_found, unsigned int cthread_count) {
  cudaError_t res;
#ifndef NDEBUG
  fprintf(stderr, "In check_ns...\n");
#endif
  // Pass P.
  res = cudaMemcpy(d_P, P, cthread_count*sizeof(uint64_t), cudaMemcpyHostToDevice);
  if(res != cudaSuccess) {
    if(res == cudaErrorInvalidValue) fprintf(stderr, "Memcpy error: Invalid value!\n");
    if(res == cudaErrorInvalidDevicePointer) fprintf(stderr, "Memcpy error: Invalid device pointer!\n");
    if(res == cudaErrorInvalidMemcpyDirection) fprintf(stderr, "Memcpy error: Invalid memcpy direction!\n");
    exit(1);
  }
#ifndef NDEBUG
  fprintf(stderr, "Memcpy successful...\n");
#endif
  d_setup_ps<<<cthread_count/BLOCKSIZE,BLOCKSIZE>>>(d_P, d_bitsskip);
#ifndef NDEBUG
  fprintf(stderr, "Setup successful...\n");
#endif
  // Pass K.
  if(cudaMemcpy(d_K, K, cthread_count*sizeof(uint64_t), cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Memcpy2 error!\n");
    exit(1);
  }
#ifndef NDEBUG
  fprintf(stderr, "Memcpy2 successful...\n");
#endif
  d_check_ns<<<cthread_count/128,128>>>(d_P, d_K, d_factor_found, d_bitsskip);
#ifndef NDEBUG
  fprintf(stderr, "Main kernel successful...\n");
#endif
  // Get d_factor_found, into the thread'th factor_found array.
  cudaMemcpy(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost);
#ifndef NDEBUG
  fprintf(stderr, "Retrieve successful...\n");
#endif
}

void cuda_finalize(void) {
  cudaFree(d_bitsskip);
  cudaFree(d_K);
  cudaFree(d_P);
  cudaFree(d_factor_found);
}
