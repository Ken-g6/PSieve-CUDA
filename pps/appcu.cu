/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.cu -- (C) Ken Brazier February - September 2010.

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
//#include <cuda.h>
//#include <assert.h>
#include "cuda_sleep_memcpy.h"
#include "main.h"
#include "putil.h"
#include "app.h"
#include "appcu.h"

#include <unistd.h>

#define INLINE static inline
/*
#ifndef BITSATATIME
#define BITSATATIME 4
#endif
#define BITSMASK ((1<<BITSATATIME)-1)*/
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
int ld_bbits;
uint64_t ld_r0;
// Stores the N to start at for the given bit position in 
//unsigned int n_subsection_start[9];
//static unsigned int *first_n_subsection = &n_subsection_start[7];

// Device constants
//__constant__ unsigned int d_bitsatatime;
//__constant__ unsigned int d_len;//=(1<<bitsatatime); 
//__constant__ unsigned int d_halflen;//=(1<<bitsatatime)/2; 
__constant__ uint64_t d_kmax;
__constant__ uint64_t d_kmin;
//__constant__ unsigned int d_bpernstep;
__constant__ unsigned int d_nmin;
__constant__ unsigned int d_nmax;
__constant__ unsigned int d_nstep;
__constant__ unsigned int d_kernel_nstep;
__constant__ unsigned int d_search_proth;

__constant__ int d_bbits;
__constant__ unsigned int d_mont_nstep;
__constant__ unsigned int d_mont_nstep_small;
__constant__ uint64_t d_r0;
// Device arrays
uint64_t *d_P;
uint64_t *d_Ps, *d_K;
//unsigned int *d_N;
unsigned char *d_factor_found;

// Timing variables:
//const int setup_ps_overlap = 5000;
const int check_ns_overlap = 50000;
bool blocking_sync_ok=true;

static unsigned int *thread_kernel_nstep;
static struct cudaDeviceProp *ccapability;
static int max_ns_delay = 0;
//static int ccapability = 1;

//globals for cuda
//cudaEvent_t stop;
//cudaStream_t stream = NULL;

// A getter for n_subsection_start.  Yes, in C!
//int get_n_subsection_start(int index) { return n_subsection_start[index]; }

// find the log base 2 of a number.  Need not be fast; only done once.
int lg2(uint64_t v) {
  int r = 0; // r will be lg(v)

  while (v >>= 1) // unroll for more speed...
  {
    r++;
  }
  return r;
}

// Clean up and terminate, whether there was an error or no.
void cuda_finalize(void) {
  //cudaFree(d_bitsskip);
  //cudaFree(d_N);
  if(d_K){
    cudaFree(d_K);
  }
  if(d_Ps){
    cudaFree(d_Ps);
  }
  if(d_P){
    cudaFree(d_P);
  }
  if(d_factor_found){
    cudaFree(d_factor_found);
  }
  //if(stop){
    //cudaEventDestroy(stop);
  //}
  //if(stream){
    //cudaStreamDestroy(stream);
  //}
}

void checkCUDAErr(const char* msg) {
  cudaError_t err = cudaGetLastError();
  if(cudaSuccess!=err) {
    fprintf(stderr, "Cuda error: %s: %s\n", msg, cudaGetErrorString(err));
    cuda_finalize();
    bexit(EXIT_FAILURE);
  }
}

bool SetCUDABlockingSync(int device) {
  cudaError_t status = cudaGetLastError();

  if(status != cudaSuccess) return false;

  status = cudaSetDevice(device);
  if(status != cudaSuccess) return false;

  status = cudaSetDeviceFlags(cudaDeviceBlockingSync);
  if(status != cudaSuccess) return false;

  return true;
}

static void sleep_for_server() {
#ifdef USE_BOINC
  printf("Sleeping 10 minutes to give the server a break.\n");
  bmsg("Sleeping 10 minutes to give the server a break.\n");
  sleep(600);
#endif
}

void cuda_init(void) {
  // Set the constants.
  //cudaMemcpyToSymbol(d_bitsatatime, &ld_bitsatatime, sizeof(ld_bitsatatime));
  max_ns_delay = (int)((nmax-nmin+1)*MAX_NS_DELAY_PER_N);
  if (ld_nstep > (nmax-nmin+1))
    ld_nstep = (nmax-nmin+1);

  //assert((1ul << (64-nstep)) < pmin);
  if((((uint64_t)1) << (64-ld_nstep)) > pmin) {
    uint64_t pmin_1 = (((uint64_t)1) << (64-ld_nstep));
    bmsg("Error: pmin is not large enough (or nmax is close to nmin).\n");
    cuda_finalize();
    while((((uint64_t)1) << (64-ld_nstep)) > pmin) {
      pmin *= 2;
      ld_nstep++;
    }
    if(pmin_1 < pmin) pmin = pmin_1;
#ifndef _WIN32
    fprintf(stderr, "This program will work by the time pmin == %lu.\n", pmin);
#endif
    bexit(ERR_INVALID_PARAM);
  }
#ifdef SEARCH_TWIN
  // For TPS, decrease the ld_nstep by one to allow overlap, checking both + and -.
  ld_nstep--;
#endif
  // Use the 32-step algorithm where useful.
  if(ld_nstep >= 32 && ld_nstep < 48 && (((uint64_t)1) << 32) <= pmin) {
    if(ld_nstep != 32) printf("nstep changed to 32\n");
    ld_nstep = 32;
  } else {
#ifdef SEARCH_TWIN
    printf("Changed nstep to %u\n", ld_nstep);
#else
    printf("Didn't change nstep from %u\n", ld_nstep);
#endif
  }

  ld_bbits = lg2(nmin);
  //assert(d_r0 <= 32);
  if(ld_bbits < 6) {
    fprintf(stderr, "%sError: nmin too small at %d (must be at least 64).\n", bmprefix(), nmin);
    cuda_finalize();
    bexit(ERR_INVALID_PARAM);
  }
  // r = 2^-i * 2^64 (mod N), something that can be done in a uint64_t!
  // If i is large (and it should be at least >= 32), there's a very good chance no mod is needed!
  ld_r0 = ((uint64_t)1) << (64-(nmin >> (ld_bbits-5)));

  ld_bbits = ld_bbits-6;

  // Allocate arrays.
  thread_kernel_nstep = (unsigned int *)xmalloc(num_threads*sizeof(unsigned int));
  ccapability = (struct cudaDeviceProp *)xmalloc(num_threads*sizeof(struct cudaDeviceProp));
}

/* This function is called once per thread.
 */
unsigned int cuda_app_init(int gpuno, int th, unsigned int cthread_count)
{
  unsigned int i, ld_kernel_nstep;
  struct cudaDeviceProp *gpuprop;
  //unsigned int ld_bitsatatime = 0;
  //unsigned int ld_halflen=(1<<bitsatatime)/2; 
  //unsigned int ld_bitsmask;
  //unsigned int ld_bpernstep;
  //unsigned int cthread_count;

  // Find the GPU's properties.
  if(cudaGetDeviceProperties(&ccapability[th], gpuno) != cudaSuccess) {
    fprintf(stderr, "%sGPU %d not compute-capable.\n", bmprefix(), gpuno);
#ifdef USE_BOINC
    fprintf(stderr, "Cuda error: getting device properties: %s\n", cudaGetErrorString(cudaGetLastError()));
    //sleep_for_server();
    bexit(ERR_NOT_IMPLEMENTED);
#else
    return 0;
#endif
  }
  gpuprop = &ccapability[th];
  /* Assume N >= 2^32. */
  if(pmin <= ((uint64_t)1)<<32) {
    bmsg("Error: PMin is too small, <= 2^32!\n");
    bexit(ERR_INVALID_PARAM);
  }
  blocking_sync_ok = SetCUDABlockingSync(gpuno);
  if(blocking_sync_ok == false) bmsg("Blocking sync setup failed; try upgrading your drivers.\n");
  //  cudaSetDevice(gpuno);
  fprintf(stderr, "%sDetected GPU %d: %s\n", bmprefix(), gpuno, gpuprop->name);
  fprintf(stderr, "%sDetected compute capability: %d.%d\n", bmprefix(), gpuprop->major, gpuprop->minor);
#ifndef _DEVICEEMU
  if(gpuprop->major == 9999 && gpuprop->minor == 9999) {
    bmsg("Detected emulator!  We can't use that!\n");
    sleep_for_server();
    bexit(ERR_NOT_IMPLEMENTED);  // System call not implemented on this platform: it's a CPU, not a GPU!
  }
#endif
    
  fprintf(stderr, "%sDetected %d multiprocessors.\n", bmprefix(), gpuprop->multiProcessorCount);
  //fprintf(stderr, "%sDetected %lu bytes of device memory.\n", bmprefix(), gpuprop->totalGlobalMem);

  // Use them to set cthread_count.
  // If cthread_count was already set, make sure it's 0 mod BLOCKSIZE.
  if(cthread_count == 0) {
    // Threads per multiprocessor, based on compute capability, if not manually set.
    cthread_count = (gpuprop->major == 1 && gpuprop->minor < 2)?384:768;
    if(gpuprop->major >= 2) cthread_count = 1024*2;
  } else {
    if(cthread_count < BLOCKSIZE) cthread_count *= BLOCKSIZE;
    else cthread_count -= cthread_count % BLOCKSIZE;
  }
  cthread_count *= gpuprop->multiProcessorCount;

  if(gpuprop->totalGlobalMem < cthread_count*(3*sizeof(uint64_t)+sizeof(unsigned char))) {
    fprintf(stderr, "%sInsufficient GPU memory: %u bytes.\n", bmprefix(), (unsigned int)(gpuprop->totalGlobalMem));
#ifdef USE_BOINC
    bexit(ERR_INSUFFICIENT_RESOURCE);
#else
    return 0;
#endif
  }
  // Calculate ld_bitsatatime given memory constraints, and possibly nmin-nmax via nstep vs. 2^ld_bitsatatime
  // Things change if nmax-nmin < 1000000 or so, but for now let's go with a constant maximum of ld_bitsatatime<=13.
  //i = gpuprop->totalGlobalMem/sizeof(uint64_t); // Total number of 64-bit numbers that can be stored.
  //ld_bitsatatime = BITSATATIME;
  //ld_bitsmask = BITSMASK+1;

  // Allocate device arrays:
  i=0;  // Iterate some number of times if not enough memory yet.
  while(1) {
    // - d_bitsskip[] (Biggest array first.)
    //if(cudaMalloc((void**)&d_bitsskip, ld_bitsmask*cthread_count*sizeof(uint64_t)) == cudaSuccess) {
    // - P's
    if(cudaMalloc((void**)&d_P, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
      // - Ps's
      if(cudaMalloc((void**)&d_Ps, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
        // - K's
        if(cudaMalloc((void**)&d_K, cthread_count*sizeof(uint64_t)) == cudaSuccess) {
          // - N's
          //if(cudaMalloc((void**)&d_N, cthread_count*sizeof(unsigned int)) == cudaSuccess) {
          // - d_factor_found[]
          if(cudaMalloc((void**)&d_factor_found, cthread_count*sizeof(unsigned char)) == cudaSuccess) {
#ifndef NDEBUG
            fprintf(stderr, "Allocation successful!\n");
            fprintf(stderr, "ld_bitsatatime = %u\n", ld_bitsatatime);
#endif
            break;  // Allocation successful!
          }
          //cudaFree(d_N);
          //}
          cudaFree(d_K);
        }
        cudaFree(d_Ps);
      }
      cudaFree(d_P);
    }
    //cudaFree(d_bitsskip);
    //}
    fprintf(stderr, "%sInsufficient available memory on GPU %d.\n", bmprefix(), gpuno);
    if(pstart == pmin || (++i) >= 86) {
#ifdef USE_BOINC
      bexit(ERR_INSUFFICIENT_RESOURCE);
#else
      return 0;
#endif
    } else {
      sleep(7);
      bmsg("Trying again...");
    }
  }

  //ld_bitsmask--; // Finalize bitsmask


  cudaMemcpyToSymbol(d_bbits, &ld_bbits, sizeof(ld_bbits));
  // d_mont_nstep is the montgomerized version of nstep.
  i = 64-ld_nstep;
  cudaMemcpyToSymbol(d_mont_nstep, &i, sizeof(i));
  i -= 32;
  cudaMemcpyToSymbol(d_mont_nstep_small, &i, sizeof(i));
  cudaMemcpyToSymbol(d_r0, &ld_r0, sizeof(ld_r0));

  // N's to search each time a kernel is run:
  ld_kernel_nstep = ITERATIONS_PER_KERNEL;
  if(ld_nstep == 32) ld_kernel_nstep /= 2;
  // Adjust for differing block sizes.
  ld_kernel_nstep *= 384;
  ld_kernel_nstep /= (cthread_count/gpuprop->multiProcessorCount);
  // Increase the step size for Fermis, to reduce memory bandwidth,
  // and because they have more threads running per multiprocessor.
  if(gpuprop->major >= 2) ld_kernel_nstep *= 3;
  // But shrink it to give at least four big N sections.
  if(ld_nstep == 32) i = 2;
  else i = 1;
  while((nmax-nmin) < 4*(ld_kernel_nstep*ld_nstep*i) && ld_kernel_nstep >= 100) ld_kernel_nstep /= 2;
  
  // Finally, make sure it's a multiple of ld_nstep!!!
  ld_kernel_nstep *= ld_nstep;
  // When ld_nstep is 32, the special algorithm there effectively needs ld_kernel_nstep divisible by 64.
  if(ld_nstep == 32) ld_kernel_nstep *= 2;

  cudaMemcpyToSymbol(d_kernel_nstep, &ld_kernel_nstep, sizeof(ld_kernel_nstep));
  thread_kernel_nstep[th] = ld_kernel_nstep;
  cudaMemcpyToSymbol(d_kmax, &kmax, sizeof(kmax));
  cudaMemcpyToSymbol(d_kmin, &kmin, sizeof(kmin));
  cudaMemcpyToSymbol(d_nmin, &nmin, sizeof(nmin));
  cudaMemcpyToSymbol(d_nmax, &nmax, sizeof(nmax));
  cudaMemcpyToSymbol(d_nstep, &ld_nstep, sizeof(ld_nstep));
  i = (search_proth == 1)?1:0;	// search_proth is 1 or -1, not 0.
  cudaMemcpyToSymbol(d_search_proth, &i, sizeof(i));

  // Initialize this n_subsection_start
  // In slot 0 we insert nmax, so I'll refer to the bits as 1-8 to avoid confusion.
  // Bit 1 is the highest part.  Usually bit 8 is the lowest, but sometimes it's a lower bit.
  {
    uint64_t test_n = nmin, next_n;
    int j;
    unsigned int *n_subsection_start = thread_subsections[th].n_subsection_start;

    thread_subsections[th].first_n_subsection = &n_subsection_start[7];
    n_subsection_start[8] = nmin;
    for(j=7; j >= 0; j--) {
      // Divide the range into 8 sub-ranges of N's.
      next_n = nmin + ((nmax - nmin + 8)/8)*(8-j);
      // Pretty inefficient, but no more so than one kernel call loop.
      while(test_n < next_n) test_n += ld_kernel_nstep;
      n_subsection_start[j] = test_n;
      // If test_n wasn't changed at all, shrink the range.
      // Horribly inefficient at O(n^2), but n == 9.
      if(test_n == n_subsection_start[j+1]) {
        thread_subsections[th].first_n_subsection--;
        for(i=j; i < 8; i++)
          n_subsection_start[i] = n_subsection_start[i+1];
      }
    }
    // Make sure bit 0 is the highest.
    if(n_subsection_start[0] < nmax) bmsg("Warning: n_subsection_start[0] too small.\n");
    n_subsection_start[0] = nmax;
  }
/*
  printf("Listing N subsections created:\n");
  for(i=0; &n_subsection_start[i] != first_n_subsection; i++)
    printf("Subsection %d: %d-%d.\n", i, get_n_subsection_start(i+1), get_n_subsection_start(i));
  printf("Subsection %d: %d-%d.\n", i, get_n_subsection_start(i+1), get_n_subsection_start(i));
*/
  //cudaStreamCreate(&stream);
  //checkCUDAErr("cudaStreamCreate");

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
    fprintf (stderr, "%sError, asm mulredc(%lu,%lu,%lu) = %lu\n", bmprefix(), a, b, N, r);
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
    fprintf (stderr, "%sError, mulredc(%lu,%lu,%lu) = %lu\n", bmprefix(), a, b, N, rax);
    cuda_finalize();
    bexit(ERR_NEG);
  }
#endif

  return rax;
}

// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
__device__ uint64_t onemod_REDC(const uint64_t N, uint64_t rax) {
  uint64_t rcx;

  // Akruppa's way, Compute T=a; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //if(rax != 0) {
    // Most of the time.
    rax = __umul64hi(rax, N) + 1;
  //} else {
    // I believe this is impossible, here, as it would lead to returning 0, which is impossible.
    //rax = __umul64hi(rax, N);
  //}
    
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  //rax = __umul64hi(rax, N) + rcx;
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
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
    fprintf (stderr, "%sError, redc(%lu,%lu) = %lu\n", bmprefix(), a, N, r);
    cuda_finalize();
    bexit(ERR_NEG);
  }

  return r;
#endif
}

// Compute T=a<<s; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
// rax is passed in as a * Ns.
// rax's original value is destroyed, just to keep the register count down.
__device__ uint64_t shiftmod_REDC (const uint64_t a, 
    const uint64_t N, uint64_t rax)
{
  uint64_t rcx;

  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax <<= d_mont_nstep; // So this is a*Ns*(1<<s) == (a<<s)*Ns (mod 2^64).
  rcx = a >> d_nstep;
  //"movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
  //"imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
  //rax *= Ns;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx += (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = __umul64hi(rax, N) + rcx;
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, ((uint64_t)1)<<d_mont_nstep, N))
  {
    fprintf (stderr, "%sError, shiftredc(%lu,%u,%lu) = %lu\n", bmprefix(), a, d_mont_nstep, N, rax);
    cuda_finalize();
    bexit(ERR_NEG);
  }
#endif

  return rax;
}

// Same function, specifically designed for various-bit shift.
// (And to not clobber rax.)
#define SHIFTMOD_REDCX(NSTEP) \
__device__ uint64_t shiftmod_REDC##NSTEP (const uint64_t a, \
    const uint64_t N, uint64_t rax) \
{ \
  uint64_t rcx; \
  rax <<= (64-NSTEP); \
  rcx = a >> NSTEP; \
  rcx += (rax!=0)?1:0; \
  rax = __umul64hi(rax, N) + rcx; \
  rcx = rax - N; \
  rax = (rax>N)?rcx:rax; \
  return rax; \
}
//SHIFTMOD_REDCX(21)
//SHIFTMOD_REDCX(32)
//SHIFTMOD_REDCX(42)

// Clobbers C to save two registers and to make the example I found work.
static inline __device__ void mad_wide_u32(const unsigned int a, const unsigned int b, uint64_t &c) {
#ifndef _DEVICEEMU
  asm("mad.wide.u32 %0, %1, %2, %0;" : "+l" (c) : "r" (a) , "r" (b));
#else
  // CPU emulation:
  c += ((uint64_t)a) * ((uint64_t)b);
#endif
}

static inline __device__ unsigned int __umul24hi(const unsigned int a, const unsigned int b) {
  unsigned int res;
#ifndef _DEVICEEMU
  asm("mul24.hi.u32 %0, %1, %2;" : "=r" (res) : "r" (a) , "r" (b));
#else
  // CPU emulation:
  res = (unsigned int)((((uint64_t)(a&0xFFFFFF)) * ((uint64_t)(b&0xFFFFFF)))>>32);
#endif
  return res;
}

// Same function, specifically designed for 32-bit shift.
// Using two mads we can save two cycles on Fermi - from the rcx += line of all places!
__device__ void shiftmod_REDC32 (uint64_t &rcx,
    const uint64_t N, unsigned int rax)
{
  //unsigned int temp;
  //uint64_t rcx;

  rcx >>= 32;
  //temp = ((unsigned int)(N>>32));
  // This isn't in an asm, but should be compiled as mad.hi.u32.  One cycle/four on older GPUs.
  // We know this can fit an unsigned int because (N-1)*(N-1) = N^2-2N-1, so adding the equivalent of 1N is OK.
  rcx += __umulhi(rax, (unsigned int)N) + (((rax)!=0)?1:0);
  // A wide multiply should take one cycle/four on older GPUs; but two wouldn't kill Fermi.
  mad_wide_u32((rax),((unsigned int)(N>>32)), rcx);
  // Two cycles for this one add!
  //rcx = ((((uint64_t)__umulhi((unsigned int)rax,(unsigned int)(N>>32))) << 32) | (((unsigned int)rax)*((unsigned int)(N>>32)))) + rcx;
  //rax = ((uint64_t)((unsigned int)rax))*((uint64_t)((unsigned int)(N>>32))) + temp;

  // And the rest is normal, but squashed.
  rcx = (rcx>N)?(rcx-N):rcx;
}

// Same function for nstep < 32.  (SMall.)
// Third argument must be passed in as only the low register, as we're effectively left-shifting 32 plus a small number.
__device__ void shiftmod_REDCsm (uint64_t &rcx, const uint64_t N, unsigned int rax)
{
  //unsigned int temp;
  //uint64_t rcx;

  rax <<= d_mont_nstep_small;
  rcx >>= d_nstep;
  //temp = ((unsigned int)(N>>32));
  // This isn't in an asm, but should be compiled as mad.hi.u32.  One cycle/four on older GPUs.
  // We know this can fit an unsigned int because (N-1)*(N-1) = N^2-2N-1, so adding the equivalent of 1N is OK.
  rcx += __umulhi(rax, (unsigned int)N) + (((rax)!=0)?1:0);
  // A wide multiply should take one cycle/four on older GPUs; but two wouldn't kill Fermi.
  mad_wide_u32((rax),((unsigned int)(N>>32)), rcx);
  // Two cycles for this one add!
  //rcx = ((((uint64_t)__umulhi((unsigned int)rax,(unsigned int)(N>>32))) << 32) | (((unsigned int)rax)*((unsigned int)(N>>32)))) + rcx;
  //rax = ((uint64_t)((unsigned int)rax))*((uint64_t)((unsigned int)(N>>32))) + temp;

  // And the rest is normal, but squashed.
  rcx = (rcx>N)?(rcx-N):rcx;
}

// A Left-to-Right version of the powmod.  Calcualtes 2^-(first 6 bits), then just keeps squaring and dividing by 2 when needed.
__device__ uint64_t
invpowmod_REDClr (const uint64_t N, const uint64_t Ns) {
  uint64_t r;
  int bbits = d_bbits;

  r = d_r0;

  // Now work through the other bits of nmin.
  for(; bbits >= 0; --bbits) {
    // Just keep squaring r.
    r = mulmod_REDC(r, r, N, Ns);
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
  r = mod_REDC (r, N, Ns);

  return r;
}

// Device-local function to iterate over some N's.
// To avoid register pressure, clobbers i, and changes all non-const arguments.

#ifdef SEARCH_TWIN
#ifdef _DEVICEEMU
#define TWIN_PRINT_NEG //fprintf(stderr, "%s%lu | %u*2^%u+/-1 (neg, P[%d])\n", bmprefix(), my_P, (unsigned int)kpos, n, blockIdx.x * BLOCKSIZE + threadIdx.x);
#else
#define TWIN_PRINT_NEG
#endif
// Select the even one here, so as to use the zero count and shift.
// The other side (whether positive or negative) is odd then, with no zeroes on the right.
#define TWIN_CHOOSE_EVEN_K0 kpos = (((unsigned int)k0) & 1)?(my_P - k0):k0;
//#define TWIN_CHOOSE_EVEN if(((unsigned int)kpos) & 1) kpos = my_P - kpos;
#define TWIN_CHOOSE_EVEN kpos = (((unsigned int)kpos) & 1)?(my_P - kpos):kpos;
// No zeroes on the right here.
// Small kmax means testing the low and high bits of kpos separately.
#define TWIN_TEST_NEG_SM
/*
#define TWIN_TEST_NEG_SM \
    kpos = my_P - kpos; \
    if (((unsigned int)(kpos>>32)) == 0 && ((unsigned int)kpos) <= ((unsigned int)d_kmax)) {\
      TWIN_PRINT_NEG \
      if(kpos >= d_kmin && n < d_nmax) my_factor_found |= 1; \
    }
*/
#else
#define TWIN_CHOOSE_EVEN_K0 kpos = k0;
#define TWIN_CHOOSE_EVEN
#define TWIN_TEST_NEG_SM
#endif

#ifdef _DEVICEEMU
#define PRINT_FACTOR_FOUND2(STAGE) //printf(STAGE "%lu | %u*2^%u+/-1 (P[%d])\n", my_P, (unsigned int)(kpos>>i), n+i, blockIdx.x * BLOCKSIZE + threadIdx.x); \
  printf("\tPrints only if %lu >= %lu, %u <(=) %u, and %u < %u\n", (kpos >> i), d_kmin, i, d_nstep, n+i, l_nmax);
#define PRINT_FACTOR_FOUND(STAGE) PRINT_FACTOR_FOUND2(STAGE)
#else
#define PRINT_FACTOR_FOUND(STAGE)
#endif

// Count trailing zeros of kpos and store that count in i.
//i = __ffsll(kpos)-1;
#define CTZLL_KPOS \
    i = (unsigned int)kpos; \
    if(i != 0) { \
      i=(__float_as_int(__uint2float_rz(i & -i))>>23)-0x7f; \
    } else { \
      i = (unsigned int)(kpos>>32); \
      i=63 - __clz (i & -i); \
    }

// Test for a factor by shifting by the result of CTZLL_KPOS
// Small kmax means testing the low and high bits of kpos separately.
// STAGE is just a string printed when printing the status while debugging.
#ifdef SEARCH_TWIN
#define D_NSTEP_COMPARE <=
#else
#define D_NSTEP_COMPARE <
#endif
#define TEST_SHIFT_SMALL_KMAX(STAGE) \
    if ((((unsigned int)(kpos>>32))>>i) == 0) \
      if(((unsigned int)(kpos>>i)) <= ((unsigned int)d_kmax)) { \
        PRINT_FACTOR_FOUND(STAGE) \
        if((kpos>>i) >= d_kmin && i D_NSTEP_COMPARE d_nstep && n+i D_NSTEP_COMPARE l_nmax) my_factor_found |= 1; \
      }

__device__ void d_check_some_ns(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos;
/*#ifdef SEARCH_TWIN
  uint64_t kneg;
#endif*/
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do { // Remaining steps are all of equal size nstep
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    CTZLL_KPOS

    if ((kpos >> i) <= d_kmax) {
#ifdef _DEVICEEMU
      //fprintf(stderr, "%s%lu | %lu*2^%u+/-1 (P[%d])\n", bmprefix(), my_P, (kpos>>i), n+i, blockIdx.x * BLOCKSIZE + threadIdx.x);
#endif
      // Just flag this if kpos <= d_kmax.
      if((kpos >> i) >= d_kmin && i D_NSTEP_COMPARE d_nstep) my_factor_found |= 1;
    }

/*
#ifdef SEARCH_TWIN
    kpos = my_P - kpos;

    if (kpos <= d_kmax) {
#ifdef _DEVICEEMU
      //fprintf(stderr, "%s%lu | %lu*2^%u+/-1 (neg, P[%d])\n", bmprefix(), my_P, kpos, n, blockIdx.x * BLOCKSIZE + threadIdx.x);
#endif
      // Just flag this if kpos <= d_kmax.
      if(kpos >= d_kmin) my_factor_found |= 1;
    }
#endif
*/

    // Proceed to the K for the next N.
    // kpos is destroyed, just to keep the register count down.
    // Despite the Montgomery step, this isn't really in Montgomery form.
    // The step just divides by 2^64 after multiplying by 2^(64-nstep).  (All mod P)
    kpos = k0 * Ps;
    n += d_nstep;   // Calculate n here so the next instruction can run in parallel on GF104.
    k0 = shiftmod_REDC(k0, my_P, kpos);
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)k0);
#endif
}

// Device-local function to iterate over some N's.
// To avoid register pressure, clobbers i, and changes all non-const arguments.
// Small kmax means testing the low and high bits of kpos separately.
__device__ void d_check_some_ns_small_kmax(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos;
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do { // Remaining steps are all of equal size nstep
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    CTZLL_KPOS

    TEST_SHIFT_SMALL_KMAX("")

    TWIN_TEST_NEG_SM
    // Proceed to the K for the next N.
    // kpos is destroyed, just to keep the register count down.
    // Despite the Montgomery step, this isn't really in Montgomery form.
    // The step just divides by 2^64 after multiplying by 2^(64-nstep).  (All mod P)
    //kpos = k0 * Ps;
    n += d_nstep;
    shiftmod_REDCsm(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)k0);
#endif
}
// Device-local function to iterate over some N's.
// To avoid register pressure, clobbers i, and changes all non-const arguments.
// This version works only for nstep (== d_mont_nstep) == 32
__device__ void d_check_some_ns_32(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos; //, kPs;
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do { // Remaining steps are all of equal size nstep
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    CTZLL_KPOS
    TEST_SHIFT_SMALL_KMAX("part 1/2: ")

    TWIN_TEST_NEG_SM
    // Skip 32 N's.
    n += 32;
    shiftmod_REDC32(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
    // Test again, because it's already set up.  Loop unrolling!
    TWIN_CHOOSE_EVEN_K0
    CTZLL_KPOS
    TEST_SHIFT_SMALL_KMAX("part 2/2: ")

    TWIN_TEST_NEG_SM
    // Skip 32 N's.
    n += 32;
    shiftmod_REDC32(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)kpos);
#endif
}

// *** Fermi check_some_ns ***
// Same as functions above, but use a different method to check for a factor.

// FERMI_TEST_SMALL_KMAX uses the "alternate algorithm" from the CPU version to pre-test for factors.
#define TEST_SHIFT_SMALL_KMAX2(STAGE) TEST_SHIFT_SMALL_KMAX(STAGE)
#define FERMI_TEST_SMALL_KMAX(STAGE) \
    i = (unsigned int)kpos; \
    if(i != 0) { \
      i=i & -i; \
      if(__umulhi((unsigned int)d_kmax, i) >= ((unsigned int)(kpos>>32))) { \
        i=(__float_as_int(__uint2float_rz(i))>>23)-0x7f; \
        TEST_SHIFT_SMALL_KMAX2("Fermi " STAGE) \
      } \
    } else { \
      i = (unsigned int)(kpos>>32); \
      i=63 - __clz (i & -i); \
      TEST_SHIFT_SMALL_KMAX2("Fermi rare case " STAGE) \
    }

__device__ void d_check_some_ns_small_kmax_fermi(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos;
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do {
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    FERMI_TEST_SMALL_KMAX("")

    TWIN_TEST_NEG_SM
    // Proceed to the K for the next N.
    // kpos is destroyed, just to keep the register count down.
    // Despite the Montgomery step, this isn't really in Montgomery form.
    // The step just divides by 2^64 after multiplying by 2^(64-nstep).  (All mod P)
    n += d_nstep;
    shiftmod_REDCsm(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)k0);
#endif
}
// Device-local function to iterate over some N's.
// To avoid register pressure, clobbers i, and changes all non-const arguments.
// This version works only for nstep (== d_mont_nstep) == 32
__device__ void d_check_some_ns_32_fermi(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos;
/*
#ifdef _DEVICEEMU
  uint64_t kPs;
  unsigned int rcx_sm;
#endif
*/
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do { // Remaining steps are all of equal size nstep
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    FERMI_TEST_SMALL_KMAX("part 1/2: ")

    TWIN_TEST_NEG_SM
    // Skip 32 N's.
    n += 32;
/*
#ifdef _DEVICEEMU
    kPs = k0*Ps;
    rcx_sm = (unsigned int)(k0>>32);
    kPs = shiftmod_REDC(k0, my_P, kPs);
#endif
*/
    shiftmod_REDC32(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
/*
#ifdef _DEVICEEMU
    if(kPs != k0) printf("Error: Expected k0=%ld, but got %ld! Diff=%ld, rcx=%d\n", kPs, k0, kPs-k0, rcx_sm);
#endif
*/
    // Test again, because it's already set up.  Loop unrolling!
    TWIN_CHOOSE_EVEN_K0
    FERMI_TEST_SMALL_KMAX("part 2/2: ")

    TWIN_TEST_NEG_SM
    // Skip 32 N's.
    n += 32;
/*
#ifdef _DEVICEEMU
    kPs = k0*Ps;
    kPs = shiftmod_REDC(k0, my_P, kPs);
#endif
*/
    shiftmod_REDC32(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
/*
#ifdef _DEVICEEMU
    if(kPs != k0) printf("Error: Expected k0=%ld, but got %ld!\n", kPs, k0);
#endif
*/
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)kpos);
#endif
}

// Same functions again, but with the non-Fermi equivalent of the Fermi test method.  May be faster still!
// These require that kmax < 2^24.
#define UNFERMI_TEST_24(STAGE) \
    if((((unsigned int)kpos) & 0xFFFFFF) != 0) { \
      i = -((unsigned int)kpos); \
      i &= ((unsigned int)kpos); \
      if(__umul24hi((unsigned int)d_kmax, i) >= ((unsigned int)(kpos>>32))) { \
        i=(__float_as_int(__uint2float_rz(i))>>23)-0x7f; \
        TEST_SHIFT_SMALL_KMAX2("UnFermi " STAGE) \
      } \
    } else { \
      CTZLL_KPOS \
      TEST_SHIFT_SMALL_KMAX2("UnFermi rare case " STAGE) \
    }

__device__ void d_check_some_ns_24(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos;
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do { // Remaining steps are all of equal size nstep
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    UNFERMI_TEST_24("")

    TWIN_TEST_NEG_SM
    // Proceed to the K for the next N.
    // kpos is destroyed, just to keep the register count down.
    // Despite the Montgomery step, this isn't really in Montgomery form.
    // The step just divides by 2^64 after multiplying by 2^(64-nstep).  (All mod P)
    //kpos = k0 * Ps;
    n += d_nstep;
    shiftmod_REDCsm(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)k0);
#endif
}
// Device-local function to iterate over some N's.
// To avoid register pressure, clobbers i, and changes all non-const arguments.
// This version works only for nstep (== d_mont_nstep) == 32
__device__ void d_check_some_ns_32_24(const uint64_t my_P, const uint64_t Ps, uint64_t &k0,
    unsigned int &n, unsigned char &my_factor_found, unsigned int &i) {
  uint64_t kpos; //, kPs;
  unsigned int l_nmax = n + d_kernel_nstep;
  if(l_nmax > d_nmax) l_nmax = d_nmax;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (unsigned int)k0, d_kernel_nstep);
#endif
  do { // Remaining steps are all of equal size nstep
    // Montgomery form doesn't matter; it's just k*2^64 mod P.
    // Get a copy of K, which isn't in Montgomery form.
    TWIN_CHOOSE_EVEN_K0
    UNFERMI_TEST_24("part 1/2: ")

    TWIN_TEST_NEG_SM
    // Skip 32 N's.
    n += 32;
    shiftmod_REDC32(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
    // Test again, because it's already set up.  Loop unrolling!
    TWIN_CHOOSE_EVEN_K0
    UNFERMI_TEST_24("part 2/2: ")

    TWIN_TEST_NEG_SM
    // Skip 32 N's.
    n += 32;
    shiftmod_REDC32(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
  } while(n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (unsigned int)kpos);
#endif
}

// *** KERNELS ***

// Start checking N's.
__global__ void d_start_ns(const uint64_t *P, uint64_t *Ps, uint64_t *K, unsigned char *factor_found_arr) {
  unsigned int n = d_nmin; // = nmin;
  unsigned int i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  uint64_t k0;
  //unsigned char my_factor_found = 0;
  uint64_t my_P, my_Ps;
  my_P = P[i];

  // Better get this done before the first mulmod.
  my_Ps = -invmod2pow_ul (my_P); /* Ns = -N^{-1} % 2^64 */

  // Calculate k0, not in Montgomery form.
  k0 = invpowmod_REDClr(my_P, my_Ps);

  //if(my_P == 42070000070587) printf("%lu^-1 = %lu (GPU)\n", my_P, my_Ps);

  if(d_search_proth) k0 = my_P-k0;

  //my_factor_found = 0;
  //d_check_some_ns(my_P, my_Ps, k0, n, my_factor_found, i);

  i = blockIdx.x * BLOCKSIZE + threadIdx.x;
  factor_found_arr[i] = 0;
  if(n < d_nmax) {
    Ps[i] = my_Ps;
    K[i] = k0;
  }
#ifdef SEARCH_TWIN
  // Search the non-even value one time only.
  if((((unsigned int)k0) & 1) == 0) k0 = my_P - k0;
  if (((unsigned int)(k0>>32)) == 0 && ((unsigned int)k0) <= ((unsigned int)d_kmax)) {
    TWIN_PRINT_NEG
      if(k0 >= d_kmin) factor_found_arr[i] = 1;
  }
#endif
}

// Continue checking N's.
#define CHECK_MORE_NS_KERNEL(TYPE) \
__global__ void d_check_more_ns##TYPE (const uint64_t *P, const uint64_t *Ps, uint64_t *K, unsigned int N, unsigned char *factor_found_arr, unsigned int shift) { \
  unsigned int n = N; \
  unsigned int i = blockIdx.x * BLOCKSIZE + threadIdx.x; \
  uint64_t k0 = K[i]; \
  unsigned char my_factor_found = factor_found_arr[i]; \
 \
  if(shift == 1) my_factor_found <<= 1; \
 \
  d_check_some_ns##TYPE (P[i], Ps[i], k0, n, my_factor_found, i); \
 \
  i = blockIdx.x * BLOCKSIZE + threadIdx.x; \
  factor_found_arr[i] = my_factor_found; \
  if(n < d_nmax) { \
    K[i] = k0; \
  } \
}

CHECK_MORE_NS_KERNEL()
// Continue checking N's for small kmax.
CHECK_MORE_NS_KERNEL(_small_kmax)
// Continue checking N's for nstep == 32
CHECK_MORE_NS_KERNEL(_32)
// Fermi versions of the above:
// Continue checking N's for small kmax.
CHECK_MORE_NS_KERNEL(_small_kmax_fermi)
// Continue checking N's for nstep == 32
CHECK_MORE_NS_KERNEL(_32_fermi)
// UnFermi versions of the above, like the Fermi versions, but for 24-bit multiplication:
// Continue checking N's for small kmax.
CHECK_MORE_NS_KERNEL(_24)
// Continue checking N's for nstep == 32
CHECK_MORE_NS_KERNEL(_32_24)

// *** Host Kernel-calling functions ***

#define CALL_LOOP(KERNEL) \
      for(n = nmin; n < nmax; n += ld_kernel_nstep) { \
        if(n >= *this_n_subsection) { \
          if(n > *this_n_subsection) fprintf(stderr, #KERNEL " Warning: N, %u, > expected N, %u\n", n, *this_n_subsection); \
          shift = 1; \
          this_n_subsection--; \
        } else shift = 0; \
        KERNEL <<<cblockcount,BLOCKSIZE,0>>>(d_P, d_Ps, d_K, n, d_factor_found, shift); \
        checkCUDAErr("kernel " #KERNEL "invocation"); \
      }
// Pass the arguments to the CUDA device, run the code, and get the results.
void check_ns(const uint64_t *P, const unsigned int cthread_count, const int th) {
  const unsigned int cblockcount = cthread_count/BLOCKSIZE;
  unsigned int n;
  unsigned int shift = 0;
  unsigned int *this_n_subsection = thread_subsections[th].first_n_subsection;
  unsigned int ld_kernel_nstep = thread_kernel_nstep[th];
  // timing variables:

  // Pass P.
  if(P != NULL) {
    cudaMemcpy(d_P, P, cthread_count*sizeof(uint64_t), cudaMemcpyHostToDevice);
    checkCUDAErr("cudaMemcpy");
  }
#ifndef NDEBUG
  bmsg("Setup successful...\n");
#endif
  //cudaEventCreate(&stop);
  //checkCUDAErr("cudaEventCreate");

  d_start_ns<<<cblockcount,BLOCKSIZE,0>>>(d_P, d_Ps, d_K, d_factor_found);
  checkCUDAErr("kernel invocation");
#ifndef NDEBUG
  bmsg("Main kernel successful...\n");
#endif
  // Continue checking until nmax is reached.
  if(kmax < (((uint64_t)1)<<31) && ld_nstep <= 32) {
    if(ccapability[th].major >= 2
#ifdef _DEVICEEMU
        //&& ccapability[th].major != 9999
#endif
) {
      // Use a Fermi kernel.
      if(ld_nstep == 32) {
        CALL_LOOP(d_check_more_ns_32_fermi)
      } else {
        CALL_LOOP(d_check_more_ns_small_kmax_fermi)
      }
    } else {
      if(kmax <= 0xFFFFFF) {
        if(ld_nstep == 32) {
          CALL_LOOP(d_check_more_ns_32_24)
        } else {
          CALL_LOOP(d_check_more_ns_24)
        }
      } else {
        if(ld_nstep == 32) {
          CALL_LOOP(d_check_more_ns_32)
        } else {
          CALL_LOOP(d_check_more_ns_small_kmax)
        }
      }
    }
    //#endif
  } else {
    CALL_LOOP(d_check_more_ns)
  }
}

void get_factors_found(unsigned char *factor_found, const unsigned int cthread_count, const uint64_t start_t, int *check_ns_delay) {
  // Get d_factor_found, into the thread'th factor_found array.
#ifdef USE_BOINC
  cudaError_t err;
  int count = 0;
#endif

  if(!blocking_sync_ok) {
    // Manually sleep-wait for the result.
    if(*check_ns_delay <= max_ns_delay) {
      cudaSleepMemcpyFromTime(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost, check_ns_delay, check_ns_overlap, start_t);
    } else {
      // Pass in zero seconds to wait, and ignore the result passed out.
      int i=0;
      cudaSleepMemcpyFromTime(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost, &i, check_ns_overlap, start_t);
    }
  } else {
    cudaMemcpy(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost);
  }
  //cudaStreamSynchronize(stream);
#ifdef USE_BOINC
  err = cudaGetLastError();
  while(err != cudaSuccess) {
    //fprintf(stderr, "Warning: A kernel failed with error %s.  Retry %d.\n", cudaGetErrorString(err), count+1);
    // Retry the Memcpy first.
    cudaMemcpy(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if(err == cudaSuccess) break;
    //fprintf(stderr, "Warning: A kernel still failed with error %s.  Retry %d.\n", cudaGetErrorString(err), count+1);
    // Retry the computation.
    check_ns(NULL, cthread_count, 0);
    count++;
    cudaMemcpy(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost);
    // If this is the last try, don't check the result here; that seems to eat it!
    if(count == 10) break;
    err = cudaGetLastError();
  }
#endif
  checkCUDAErr("getting factors found");

#ifndef NDEBUG
  bmsg("Retrieve successful...\n");
#endif
}
