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
#include "cuda_sleep_memcpy.h"

#define INLINE static inline
/*
#ifndef BITSATATIME
#define BITSATATIME 4
#endif
#define BITSMASK ((1<<BITSATATIME)-1)*/
// BLOCKSIZE should be a power of two for greatest efficiency.
#define BLOCKSIZE 128
/*
#if(BITSATATIME == 3)
  #define SHIFT_CAST unsigned int
#elif(BITSATATIME == 4)
  #define SHIFT_CAST uint64_t
#else
  #error "Invalid BITSATATIME."
#endif
*/
// Extern vars in appcu.h:
unsigned int ld_nstep;

// Device constants
//__constant__ unsigned int d_bitsatatime;
//__constant__ unsigned int d_len;//=(1<<bitsatatime); 
//__constant__ unsigned int d_halflen;//=(1<<bitsatatime)/2; 
__constant__ uint64_t d_kmax;
//__constant__ unsigned int d_bpernstep;
__constant__ unsigned int d_nmin;
__constant__ unsigned int d_nmax;
__constant__ unsigned int d_nstep;
__constant__ unsigned int d_search_proth;

__constant__ int d_bbits;
__constant__ unsigned int d_mont_nstep;
__constant__ uint64_t d_r0;
// Device arrays
uint64_t *d_P;
//uint64_t *d_K;
//uint64_t *d_bitsskip;
unsigned char *d_factor_found;

// Timing variables:
//const int setup_ps_overlap = 5000;
const int check_ns_overlap = 50000;

// find the log base 2 of a number.  Need not be fast; only done once.
int lg2(uint64_t v) {
	int r = 0; // r will be lg(v)

	while (v >>= 1) // unroll for more speed...
	{
		r++;
	}
	return r;
}

/* This function is called once before any threads are started.
 */
unsigned int cuda_app_init(int gpuno)
{
  unsigned int i;
  struct cudaDeviceProp gpuprop;
  //unsigned int ld_bitsatatime = 0;
  //unsigned int ld_halflen=(1<<bitsatatime)/2; 
  //unsigned int ld_bitsmask;
  //unsigned int ld_bpernstep;
  uint64_t ld_r0;
  int bbits;
  unsigned int cthread_count;

  // Find the GPU's properties.
  if(cudaGetDeviceProperties(&gpuprop, gpuno) != cudaSuccess) {
    fprintf(stderr, "GPU %d not compute-capable.\n", gpuno);
    return 0;
  }
  /* Assume N >= 2^32. */
  if(pmin <= ((uint64_t)1)<<32) {
    fprintf(stderr, "Error: PMin is too small, <= 2^32!\n");
    exit(1);
  }
  cudaSetDevice(gpuno);
  fprintf(stderr, "Detected GPU %d: %s\n", gpuno, gpuprop.name);
  fprintf(stderr, "Detected compute capability: %d.%d\n", gpuprop.major, gpuprop.minor);
  fprintf(stderr, "Detected %d multiprocessors.\n", gpuprop.multiProcessorCount);
  //fprintf(stderr, "Detected %lu bytes of device memory.\n", gpuprop.totalGlobalMem);

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
  //ld_bitsatatime = BITSATATIME;
  //ld_bitsmask = BITSMASK+1;

  // Allocate device arrays:
  // TODO: fix this awkward construct.
  while(1) {
    // - d_bitsskip[] (Biggest array first.)
    // Not using cudaMallocPitch because coalescing isn't possible in general.
    //if(cudaMalloc((void**)&d_bitsskip, ld_bitsmask*cthread_count*sizeof(uint64_t)) == cudaSuccess) {
      // - P's
      if(cudaMalloc((void**)&d_P, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
        // - K's
        //if(cudaMalloc((void**)&d_K, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
          // - d_factor_found[]
          if(cudaMalloc((void**)&d_factor_found, cthread_count*sizeof(unsigned char)) == cudaSuccess) {
#ifndef NDEBUG
            fprintf(stderr, "Allocation successful!\n");
            //fprintf(stderr, "ld_bitsatatime = %u\n", ld_bitsatatime);
#endif
            break;  // Allocation successful!
          }
          //cudaFree(d_K);
        //}
        cudaFree(d_P);
      }
      //cudaFree(d_bitsskip);
    //}
    fprintf(stderr, "Insufficient available memory on GPU %d.\n", gpuno);
    return 0;
  }

  //ld_bitsmask--; // Finalize bitsmask

  if (ld_nstep > (nmax-nmin+1))
    ld_nstep = (nmax-nmin+1);

  //assert((1ul << (64-nstep)) < pmin);
  if((((uint64_t)1) << (64-ld_nstep)) > pmin) {
    fprintf(stderr, "Error: pmin is not large enough (or nmax is close to nmin).\n");
    exit(1);
  }
  // Set the constants.
  //cudaMemcpyToSymbol(d_bitsatatime, &ld_bitsatatime, sizeof(ld_bitsatatime));

  // Prepare constants:
  bbits = lg2(nmin);
  //assert(d_r0 <= 32);
  if(bbits < 6) {
    fprintf(stderr, "Error: nmin too small at %d (must be at least 64).\n", nmin);
    exit(1);
  }
  // r = 2^-i * 2^64 (mod N), something that can be done in a uint64_t!
  // If i is large (and it should be at least >= 32), there's a very good chance no mod is needed!
  ld_r0 = ((uint64_t)1) << (64-(nmin >> (bbits-5)));

  bbits = bbits-6;
  cudaMemcpyToSymbol(d_bbits, &bbits, sizeof(bbits));
  // d_mont_nstep is the montgomerized version of nstep.
  i = 64-ld_nstep;
  cudaMemcpyToSymbol(d_mont_nstep, &i, sizeof(i));
  cudaMemcpyToSymbol(d_r0, &ld_r0, sizeof(ld_r0));

  // The following would be "pitch" if using cudaMallocPitch above.
  //i = ld_bitsmask+1;
  ////cudaMemcpyToSymbol(d_len, &i, sizeof(i));//=(1<<bitsatatime); 
  // But the following would not.
  //i >>= 1;
  //cudaMemcpyToSymbol(d_halflen, &i, sizeof(i));//=(1<<bitsatatime)/2; 
  cudaMemcpyToSymbol(d_kmax, &kmax, sizeof(kmax));
  //cudaMemcpyToSymbol(d_bpernstep, &ld_bpernstep, sizeof(ld_bpernstep));
  cudaMemcpyToSymbol(d_nmin, &nmin, sizeof(nmin));
  cudaMemcpyToSymbol(d_nmax, &nmax, sizeof(nmax));
  cudaMemcpyToSymbol(d_nstep, &ld_nstep, sizeof(ld_nstep));
  cudaMemcpyToSymbol(d_search_proth, &search_proth, sizeof(search_proth));


  return cthread_count;
}


/*** Kernel Helpers ***/
// Special thanks to Alex Kruppa for introducing me to Montgomery REDC math!
/* Compute a^{-1} (mod 2^(32 or 64)), according to machine's word size */

// Inside this #ifdef is code used only to check the faster code below it.
#ifndef NDEBUG
#ifdef __x86_64__
#define DEBUG64
/* Reduce a*2^64+b modulo m. Requires a < m, or the quotient (which we don't care about but the chip does) will overflow. */ 
__device__ uint64_t
longmod (uint64_t a, uint64_t b, const uint64_t m)
{
  //ASSERT (a < m);
  __asm__
  ( "divq %2"
    : "+d" (a), /* Put "a" in %rdx, will also get result of mod */
      "+a" (b)  /* Put "b" in %rax, will also be written to 
                   (quotient, which we don't need) */
    : "rm" (m)  /* Modulus can be in a register or memory location */
    : "cc"      /* Flags are clobbered */
  );
  return a;
}

//Now, for modulus!  From http://www.loria.fr/~kruppaal/factorcyc.20090612.c
/* Multiply a and b, and reduce the product modulo m. Remainder is
   returned */
// Be careful if a and b are >> m, as the quotient overflow from longmod could happen here too. :(
// But if either one is <= m, it's fine.
//#ifdef DEBUG
__device__ uint64_t mulmod (uint64_t a, const uint64_t b, const uint64_t m)
{
  uint64_t q, r, t1, t2;
  __asm__
  ( "mulq %3\n\t"
    : "=a" (t1), "=d" (t2)
    : "0" (a), "rm" (b)
    : "cc");
  __asm__
  ( "divq %4"
    : "=a" (q), "=d" (r)
    : "0" (t1), "1" (t2), "rm" (m)
    : "cc"
  );
  return r;
} 


/* Compute REDC(a*b) for modulus N. We need N*Ns == -1 (mod 2^64) */
__device__ uint64_t
asm_mulmod_REDC (const uint64_t a, const uint64_t b, 
             const uint64_t N, const uint64_t Ns)
{
	uint64_t r;

	// Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
	__asm__
		( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
		  "movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
		  "imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
		  "cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
		  "sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
		  "mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
		  "lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
		  "subq %[N], %%rcx\n\t"	//				Cycle 20/19?
		  "addq %%rdx, %%rcx\n\t"	//				Cycle 21/20?
		  "cmovcq %%rcx, %[r]\n\t"	//				Cycle 22/21?
		  : [r] "=r" (r)
		  : "%a" (a), [b] "rm" (b), [N] "rm" (N), [Ns] "rm" (Ns)
		  : "cc", "%rcx", "%rdx"
		);

#ifdef DEBUG64
	if (longmod (r, 0, N) != mulmod(a, b, N))
	{
		fprintf (stderr, "Error, asm mulredc(%lu,%lu,%lu) = %lu\n", a, b, N, r);
		abort();
	}
#endif

   return r;
}
#endif
#endif
__device__ uint64_t
invmod2pow_ul (const uint64_t n)
{
  uint64_t r;
  //unsigned int ir;
  const unsigned int in = (unsigned int)n;

  //ASSERT (n % 2UL != 0UL);
  
  // Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
  // correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
  // Newton iterations are enough.
  r = (n+n+n) ^ ((uint64_t)2);
  // Newton iteration
  r += r - (unsigned int) r * (unsigned int) r * in;
  r += r - (unsigned int) r * (unsigned int) r * in;
  r += r - (unsigned int) r * (unsigned int) r * in;
  r += r - r * r * n;

  return r;
}
/*
__device__ uint64_t
invmod2pow_ul (const uint64_t n)
{
  uint64_t r;

  //ASSERT (n % 2UL != 0UL);
  
  // Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
  // correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
  // Newton iterations are enough.
  r = (((uint64_t)3) * n) ^ ((uint64_t)2);
  // Newton iteration
  r += r - (unsigned int) r * (unsigned int) r * (unsigned int)n;
  r += r - (unsigned int) r * (unsigned int) r * (unsigned int)n;
  //if (sizeof (uint64_t) == 8)
  r += r - (unsigned int) r * (unsigned int) r * (unsigned int)n;
  r += r - r * r * n;

  return r;
}
*/
__device__ uint64_t mulmod_REDC (const uint64_t a, const uint64_t b, 
             const uint64_t N, const uint64_t Ns)
{
  uint64_t rax, rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax = a*b;
  rcx = __umul64hi(a,b);
  //"movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
  //rcx = rdx;
  //"imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
  rax *= Ns;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx += (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = __umul64hi(rax, N);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rax += rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, b, N))
  {
    fprintf (stderr, "Error, mulredc(%lu,%lu,%lu) = %lu\n", a, b, N, rax);
    exit(1);
  }
#endif

  return rax;
}

// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
__device__ uint64_t onemod_REDC(const uint64_t N, uint64_t rax) {
  uint64_t rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //rcx = 0;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx = (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = __umul64hi(rax, N);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rax += rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  return rax;
}

// Like mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
__device__ uint64_t mod_REDC(const uint64_t a, const uint64_t N, const uint64_t Ns) {
#ifndef DEBUG64
  return onemod_REDC(N, Ns*a);
#else
  const uint64_t r = onemod_REDC(N, Ns*a);

  if (longmod (r, 0, N) != mulmod(a, 1, N)) {
    fprintf (stderr, "Error, redc(%lu,%lu) = %lu\n", a, N, r);
    exit(1);
  }

  return r;
#endif
}

// Compute T=a<<s; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
__device__ uint64_t shiftmod_REDC (const uint64_t a, const unsigned int s, 
             const uint64_t N, const uint64_t Ns)
{
  uint64_t rax, rcx;

  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax = a << s;
  rcx = a >> (64-s);
  //"movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
  //"imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
  rax *= Ns;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx += (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = __umul64hi(rax, N);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rax += rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, ((uint64_t)1)<<s, N))
  {
    fprintf (stderr, "Error, shiftredc(%lu,%u,%lu) = %lu\n", a, s, N, rax);
    exit(1);
  }
#endif

  return rax;
}

// Hybrid powmod, sidestepping several loops and possible mispredicts, and with no more than one longmod!
/* Compute (2^-1)^b (mod m), using Montgomery arithmetic. */
// From NewPGen: NewPGen is much faster if the base is 2. This is because division by 2 modulo a prime is easy (you shift it right if it is even, otherwise you add the prime then shift it). Division by other bases isn't so straightforward, however.
// This function takes the modular inverse of 2 mod P, then exponentiates.  2^-1 mod P = (P+1)/2
// The other thing about this is that when you cross 2^-1 mod N with Montgomery multiplication,
// 2^-32 * 2^64 (mod N) = 2^32 (mod N).  So that kills the last mod if N > 2^32.
// This version runs two N's at once.
// Doing two or more invmod2pow_ul's at once is a little faster.
// Doing two mulmod_REDCs at once is even a little faster.
// This program works with K*2^N, N constant == d_nmin.
// A Left-to-Right version of the powmod.  Calcualtes 2^-(first 6 bits), then just keeps squaring and dividing by 2 when needed.


__device__ uint64_t
invpowmod_REDClr (const uint64_t N, const uint64_t Ns) {
  uint64_t r;
  int bbits = d_bbits;

  r = d_r0;

  // Now work through the other bits of nmin.
  for(; bbits >= 0; --bbits) {
    // Just keep squaring r.
    mulmod_REDC(r, r, N, Ns);
    // If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod N).
    if(d_nmin & (1u << bbits)) {
      r += (r&1)?N:0;
      r >>= 1;
    }
  }

#ifdef DEBUG64
  //assert (mod_REDC (r, N, Ns) == invmod(powmod (d_nmin, N), N));
#endif

  // Convert back to standard.
  //r = mod_REDC (r, N, Ns);

  return r;
}

// Check all N's.
__global__ void d_check_ns(const uint64_t *P, unsigned char *factor_found_arr) {
  unsigned int n = d_nmin; // = nmin;
  unsigned int i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  uint64_t k0;
  uint64_t kpos;
  unsigned char my_factor_found = 0;
  uint64_t my_P, Ps;
  my_P = P[i];
  
  assert(longmod(1, 0, 47) == 25);
  // Better get this done before the first mulmod.
  Ps = -invmod2pow_ul (my_P); /* Ns = -N^{-1} % 2^64 */
  
  // Calculate k0, in Montgomery form.
  k0 = invpowmod_REDClr(my_P, Ps);

  if(d_search_proth) k0 = my_P-k0;

  my_factor_found = 0;
  do { // Remaining steps are all of equal size nstep
    // Get K from the Montgomery form.
    kpos = mod_REDC(k0, my_P, Ps);
    i = __ffsll(kpos)-1;

    kpos >>= i;
    if (kpos <= d_kmax) {
#ifndef NDEBUG
      fprintf(stderr, "%u | %u*2^%u+1 (P[%d])\n", (unsigned int)my_P, (unsigned int)kpos, n+i, blockIdx.x * BLOCKSIZE + threadIdx.x);
#endif
      // Just flag this if kpos <= d_kmax.
      my_factor_found = 1;
    }

    // Proceed to the K for the next N.
    k0 = shiftmod_REDC(k0, d_mont_nstep, my_P, Ps);
    n += d_nstep;
  } while (n < d_nmax);
  factor_found_arr[blockIdx.x * BLOCKSIZE + threadIdx.x] = my_factor_found;
}

// Pass the arguments to the CUDA device, run the code, and get the results.
void check_ns(const uint64_t *P, uint64_t *K, unsigned char *factor_found, unsigned int cthread_count) {
  // timing variables:
  static __thread int check_ns_delay = 0;
  cudaError_t res;
  // Pass P.
  res = cudaMemcpy(d_P, P, cthread_count*sizeof(uint64_t), cudaMemcpyHostToDevice);
  if(res != cudaSuccess) {
    if(res == cudaErrorInvalidValue) fprintf(stderr, "Memcpy error: Invalid value!\n");
    if(res == cudaErrorInvalidDevicePointer) fprintf(stderr, "Memcpy error: Invalid device pointer!\n");
    if(res == cudaErrorInvalidMemcpyDirection) fprintf(stderr, "Memcpy error: Invalid memcpy direction!\n");
    exit(1);
  }
#ifndef NDEBUG
  fprintf(stderr, "Setup successful...\n");
#endif
  // Don't Pass K; it's calculated internally!
  d_check_ns<<<cthread_count/128,128>>>(d_P, d_factor_found);
#ifndef NDEBUG
  fprintf(stderr, "Main kernel successful...\n");
#endif
  // Get d_factor_found, into the thread'th factor_found array.
  cudaSleepMemcpy(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost, &check_ns_delay, check_ns_overlap);
#ifndef NDEBUG
  fprintf(stderr, "Retrieve successful...\n");
#endif
}

void cuda_finalize(void) {
  //cudaFree(d_bitsskip);
  //cudaFree(d_K);
  cudaFree(d_P);
  cudaFree(d_factor_found);
}
