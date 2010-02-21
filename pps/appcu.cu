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

// Extern vars in appcu.h:
unsigned int ld_nstep;
//unsigned int gpuno;

// Device constants
__constant__ unsigned int d_bitsatatime;
__constant__ unsigned int d_len;//=(1<<bitsatatime); 
__constant__ unsigned int d_halflen;//=(1<<bitsatatime)/2; 
__constant__ uint64_t d_kmax;
__constant__ unsigned int d_bitsmask;
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
  cthread_count = (gpuprop.major == 1 && gpuprop.minor < 2)?512:1024;
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
  ld_bitsatatime = 1;
  ld_bitsmask = 2;
  while(ld_bitsmask*cthread_count < i) {
    ld_bitsmask <<= 1;
    ld_bitsatatime++;
    if(ld_bitsatatime >= 13) break;
  }

  if(ld_bitsatatime <= 1 || ld_nstep/ld_bitsatatime >= 12) {
    fprintf(stderr, "Insufficient memory on GPU %d.\n", gpuno);
    return 0;
  }

  // Allocate device arrays:
  while(ld_bitsatatime >= 2) {
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
    // If things didn't fit, try a smaller ld_bitsatatime.
    ld_bitsatatime--;
    ld_bitsmask >>= 1;
  }

  if(ld_bitsatatime <= 1 || ld_nstep/ld_bitsatatime >= 12) {
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
  cudaMemcpyToSymbol(d_bitsmask, &ld_bitsmask, sizeof(ld_bitsmask));
  cudaMemcpyToSymbol(d_bpernstep, &ld_bpernstep, sizeof(ld_bpernstep));
  cudaMemcpyToSymbol(d_nmin, &nmin, sizeof(nmin));
  cudaMemcpyToSymbol(d_nmax, &nmax, sizeof(nmax));
  cudaMemcpyToSymbol(d_nstep, &ld_nstep, sizeof(ld_nstep));
  cudaMemcpyToSymbol(d_search_proth, &search_proth, sizeof(search_proth));

  return cthread_count;
}



// Set up the lookup tables for all P's.
__global__ void d_setup_ps(const uint64_t *P, uint64_t *bitsskip) {
  unsigned int n; // = nmin;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  uint64_t *bs0 = &bitsskip[i*d_len];
  uint64_t k0;
  uint64_t kpos;
  unsigned char my_factor_found = 0;
  // Initialize bitsskip array.
  k0 = P[i];

  // Initialize the first two entries.
  bs0[d_halflen] = (k0+1)/2;	// Needed first.
  bs0[0] = 0;			// Ignored through the end.

  // Fill in the intervening spaces, two numbers at a time.
  for(i=d_halflen; i > 1; i >>= 1) {
    for(n=i/2; n < d_halflen; n+=i) {
      kpos = bs0[2*n];
      my_factor_found = ((unsigned int)kpos)&1;
      //printf("Filling n=%d from bs0=%lu\n", n, kpos);
      bs0[n] = (kpos+((my_factor_found)?k0:(uint64_t)0))/2;
      //printf("Filling n=%d\n", n+d_halflen);
      bs0[n+d_halflen] = (kpos+1+((my_factor_found)?(uint64_t)0:k0))/2;
    }
  }
}

// Check all N's.
__global__ void d_check_ns(const uint64_t *P, const uint64_t *K, unsigned char *factor_found_arr, uint64_t *bitsskip) {
  unsigned int n = d_nmin; // = nmin;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  uint64_t *bs0 = &bitsskip[i*d_len];
  uint64_t k0;
  uint64_t kpos;
  unsigned char my_factor_found = 0;
#ifndef NDEBUG
  uint64_t my_P = P[i];
#endif

  //factor_found_arr[i] = 0;
  k0 = K[i];
  if(d_search_proth) k0 = P[i]-k0;
  my_factor_found = 0;
  do { // Remaining steps are all of equal size nstep
    kpos = k0;
    i = __ffsll(kpos)-1;

    kpos >>= i;
    if (kpos <= d_kmax) {
#ifndef NDEBUG
      fprintf(stderr, "%u | %u*2^%u+1 (P[%d])\n", (unsigned int)my_P, (unsigned int)kpos, n+i, blockIdx.x * blockDim.x + threadIdx.x);
#endif
      // Just flag this if kpos <= d_kmax.
      my_factor_found = 1;
    }

    for(i=0; i < d_bpernstep; i++) {
      k0 = (k0 >> d_bitsatatime) + bs0[(unsigned int)k0 & d_bitsmask];
    }
    n += d_nstep;
  } while (n < d_nmax);
  factor_found_arr[blockIdx.x * blockDim.x + threadIdx.x] = my_factor_found;
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
  d_setup_ps<<<cthread_count/128,128>>>(d_P, d_bitsskip);
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
