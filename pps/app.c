/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.c -- (C) Geoffrey Reynolds, April-August 2009.
 * With improvements by Ken Brazier August-October 2009.

   Proth Prime Search sieve (for many K).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <ctype.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif
#ifdef __SSE2__
#define EMM
#include <emmintrin.h>
#ifdef __x86_64__
#include <time.h>
#endif
#endif
#include "main.h"
#include "util.h"
#include "app.h"
#include "appcu.h"
#define INLINE static inline

#ifndef __x86_64__
// Macro bsfq (Bit Search Forward Quadword) for 32-bit.
// MPOS = result (32-bit)
// KPOS = input (64-bit; evaluated twice!)
// ID = a unique ID string.
#ifdef __i386__
#define BSFQ(MPOS,KPOS,ID)          asm volatile \
            ( 		"bsfl	%[k0l], %[pos]	\n" \
                "	jnz	bsflok"ID"		\n" \
                "	bsfl	%[k0h], %[pos]	\n" \
                "	addl	$32, %[pos]	\n" \
                "bsflok"ID":" \
                : [pos] "=r" (MPOS) \
                : [k0l] "rm" ((unsigned int)(KPOS)), \
                [k0h] "rm" ((unsigned int)((KPOS) >> 32)) \
                : "cc" )
#else
// If anyone wants to compile on some non-x86 platform somehow...
#define BSFQ(MPOS,KPOS,ID) MPOS = __builtin_ctzll(KPOS);
#endif
#endif

#define FORMAT_NEWPGEN 1
#define FORMAT_ABCD 2

uint64_t kmin = 0, kmax = 0;
#ifdef EMM
static uint64_t xkmax[2] __attribute__((aligned(16)));
static int sse2_in_range = 0;
#endif
static uint64_t b0 = 0, b1 = 0;
static unsigned char **bitmap = NULL;
static const char *input_filename = NULL;
static const char *factors_filename = NULL;
static FILE *factors_file = NULL;
unsigned int nmin = 0, nmax = 0;
static unsigned int nstart = 0, nstep = 0;
static unsigned int factor_count = 0;
static int file_format = FORMAT_ABCD;
static int print_factors = 1;
/*
// use_sse2: 
// 0: Use default algorithm
// 1: Use alternate SSE2 algorithm
// 2: (Default) Benchmark and find the best algorithm.
static int use_sse2 = 2;*/
unsigned int search_proth = 1;  // Search for Proth or Riesel numbers?
static unsigned int bitsatatime = 8; // Bits to process at a time, with v0.4 algorithm.
static unsigned int bitsmask, bpernstep;
static uint64_t** bitsskip;
static unsigned char** factor_found;

#ifdef _WIN32
static CRITICAL_SECTION factors_mutex;
#else
static pthread_mutex_t factors_mutex;
#endif


static void report_factor(uint64_t p, uint64_t k, unsigned int n, int c)
{
#ifdef _WIN32
  EnterCriticalSection(&factors_mutex);
#else
  pthread_mutex_lock(&factors_mutex);
#endif

  if (factors_file != NULL)
  {
    fprintf(factors_file,"%"PRIu64" | %"PRIu64"*2^%u%+d\n",p,k,n,c);
    if(print_factors) printf("%"PRIu64" | %"PRIu64"*2^%u%+d\n",p,k,n,c);
  }
  else printf("UNSAVED: %"PRIu64" | %"PRIu64"*2^%u%+d\n",p,k,n,c);
  factor_count++;

#ifdef _WIN32
  LeaveCriticalSection(&factors_mutex);
#else
  pthread_mutex_unlock(&factors_mutex);
#endif
}

// 1 if a number mod 15 is not divisible by 2 or 3.
//                             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
static const int prime15[] = { 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1 };

static void __attribute__((noinline))
test_factor(uint64_t p, uint64_t k, unsigned int n, int c)
{
  uint64_t b = k/2;

  if((k & 1) && n >= nmin) { // k is odd.
    if (bitmap == NULL) {
      // Check that K*2^N+/-1 is not divisible by 3, 5, or 7, to minimize factors printed.
      // We do 3 and 5 at the same time (15 = 2^4-1), then 7 (=2^3-1).
      // (k*(1<<(n%2))+c)%3 == 0
      if(prime15[(unsigned int)(((k<<(n&3))+(uint64_t)c)%(uint64_t)15)] && 
          (unsigned int)(((k<<(n%3))+(uint64_t)c)%(uint64_t)7) != 0)
        report_factor(p,k,n,c);
    } else {
      if (bitmap[n-nmin][(unsigned int)((b-b0)/8)] & (1<<(b-b0)%8))
        report_factor(p,k,n,c);
    }
  }
}

/* Scan the input file to determine format and values for kmin,kmax,nmin,nmax.
 */
static FILE* scan_input_file(const char *fn)
{
  FILE *file;
  uint64_t k0, k1, k, p0;
  unsigned int n0, n1, d, n;
  char ch;

  if ((file = fopen(fn,"r")) == NULL)
  {
    perror(fn);
    exit(EXIT_FAILURE);
  }


  if (fscanf(file,"ABCD %"SCNu64"*2^$a+1 [%u]",
             &k,&n) == 2)
    file_format = FORMAT_ABCD;
  else if (fscanf(file,"%"SCNu64":P:%*c:2:%c",&p0,&ch) == 2 && ch == '1')
  {
    file_format = FORMAT_NEWPGEN;
    if (fscanf(file," %"SCNu64" %u",&k,&n) != 2)
    {
      fprintf(stderr,"Invalid line 2 in input file `%s'\n",fn);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr,"Invalid header in input file `%s'\n",fn);
    exit(EXIT_FAILURE);
  }

  k0 = k1 = k;
  n0 = n1 = n;

  if (file_format == FORMAT_ABCD)
  {
    //while(getc(file) != '\n');
    printf("Scanning ABCD file...\n");
    while (1)
    {
      while(getc(file) != '\n');
      while(1) {
        char c = getc(file);
        if(!isdigit(c)) {
          ungetc(c, file);
          break;
        }
        d = c-'0';
        while(isdigit(c=getc(file))) {
          d *= 10;
          d += c-'0';
        }
        n += d;
        while(c != '\n') c=getc(file);
      }
      //while (fscanf(file," %u",&d) == 1)

      if (n1 < n)
        n1 = n;

      if (fscanf(file, " ABCD %"SCNu64"*2^$a+1 [%u]",
                 &k,&n) == 2)
      {
        if((((int)k)&15) == 1) printf("\rFound K=%"SCNu64"\r", k);
        fflush(stdout);
        if (k0 > k)
          k0 = k;
        if (k1 < k)
          k1 = k;

        if (n0 > n)
          n0 = n;
        if (n1 < n)
          n1 = n;
      }
      else
        break;
    }
  }
  else /* if (file_format == FORMAT_NEWPGEN) */
  {
    while (fscanf(file," %"SCNu64" %u",&k,&n) == 2)
    {
      if (k0 > k)
        k0 = k;
      if (k1 < k)
        k1 = k;

      if (n0 > n)
        n0 = n;
      if (n1 < n)
        n1 = n;
    }
  }

  if (ferror(file))
  {
    fprintf(stderr,"Error reading input file `%s'\n",fn);
    exit(EXIT_FAILURE);
  }

  rewind(file);
  printf("Found K's from %"SCNu64" to %"SCNu64".\n", k0, k1);
  printf("Found N's from %u to %u.\n", n0, n1);

  //if (file_format == FORMAT_ABCD)
  //{
    //k0 = 6*k0+3;
    //k1 = 6*k1+3;
  //}

  if (kmin < k0)
    kmin = k0;
  if (kmax == 0 || kmax > k1)
    kmax = k1;

  if (nmin < n0)
    nmin = n0;
  if (nmax == 0 || nmax > n1)
    nmax = n1;
  return file;
}

static void read_newpgen_file(const char *fn, FILE* file)
{
  //FILE *file;
  uint64_t k, p0;
  unsigned int n, line, count;
  char ch;

  if(file == NULL) {
    if ((file = fopen(fn,"r")) == NULL)
    {
      perror(fn);
      exit(EXIT_FAILURE);
    }
  }

  if (fscanf(file," %"SCNu64":P:%*c:2:%c",&p0,&ch) != 2 || ch != '1')
  {
    fprintf(stderr,"Invalid header in input file `%s'\n",fn);
    exit(EXIT_FAILURE);
  }

  line = 0;
  count = 0;
  while (fscanf(file," %"SCNu64" %u",&k,&n) == 2)
  {
    line++;
    if ((k&1) != 1)
    {
      fprintf(stderr,"Invalid line %u in input file `%s'\n",line,fn);
      exit(EXIT_FAILURE);
    }
    if (k >= kmin && k <= kmax && n >= nmin && n <= nmax)
    {
      uint64_t bit = k/2-b0;
      bitmap[n-nmin][(unsigned int)(bit/8)] |= (1 << bit%8);
      count++; /* TODO: Don't count duplicates */
    }
  }

  if (ferror(file))
  {
    fprintf(stderr,"Error reading input file `%s'\n",fn);
    exit(EXIT_FAILURE);
  }

  //rewind(file);
  fclose(file);

  printf("Read %u terms from NewPGen format input file `%s'\n",count,fn);
}

static void read_abcd_file(const char *fn, FILE *file)
{
  //FILE *file;
  //char buf[80];
  uint64_t k;
  unsigned int n, count, d;

  if(file == NULL) {
    printf("Opening file %s\n", fn);
    if ((file = fopen(fn,"r")) == NULL)
    {
      perror(fn);
      exit(EXIT_FAILURE);
    }
  }
  if (fscanf(file, "ABCD %"SCNu64"*2^$a+1 [%u]",
        &k,&n) != 2)
  {
    fprintf(stderr,"Invalid header in input file `%s'\n",fn);
    exit(EXIT_FAILURE);
  }

  count = 0;
  while(getc(file) != '\n');
  printf("Reading ABCD file.\n");
  while (1)
  {
    uint64_t bit = (k-kmin)/2;
    unsigned int bo8 = (unsigned int)(bit/8);
    unsigned int bm8 = (unsigned int)(1 << bit%8);
    /*if(k < kmin || k > kmax || bit < 0 || bit > (kmax-kmin)/2) {
      printf("\n\nK error: K = %"SCNu64", which is outside %"SCNu64" - %"SCNu64"\n\n\n", k, kmin, kmax);
      exit(EXIT_FAILURE);
    }*/
    if(n >= nmin) bitmap[n-nmin][bo8] |= bm8;
    count++;
    while(getc(file) != '\n');
    //while (fscanf(file," %u",&d) == 1)
    while(1)
    {
      char c = getc(file);
      if(!isdigit(c)) {
        ungetc(c, file);
        break;
      }
      d = c-'0';
      while(isdigit(c=getc(file))) {
        d *= 10;
        d += c-'0';
      }

      n += d;
      /*if(n > nmax) {
        printf("\n\nN error: N = %u, but nmax = %u\n\n\n", n, nmax);
        if(file == NULL) printf("\n\nError: File was closed!\n");
        exit(EXIT_FAILURE);
      }*/
      if(n >= nmin && n <= nmax) bitmap[n-nmin][bo8] |= bm8;
      count++;
      while(c != '\n') c=getc(file);
    }
    if((((int)k)&15) == 1) printf("\rRead K=%"SCNu64"\r", k);
    fflush(stdout);

    if (fscanf(file, "ABCD %"SCNu64"*2^$a+1 [%u]",
          &k,&n) != 2) {
      break;
    }
  }

  if (ferror(file))
  {
    printf("\nError reading input file `%s'\n",fn);
    exit(EXIT_FAILURE);
  }
  //printf("\n\nDone reading ABCD file!\n");

  fclose(file);

  printf("Read %u terms from ABCD format input file `%s'\n",count,fn);
}

/* This function is called once before anything is done. It should at least
   print the name of the application. It can also assign default option
   values which will be overridden in the configuration file or on the
   command line.
*/
void app_banner(void)
{
  printf("ppsieve version " APP_VERSION " (testing)\n");
#ifdef __GNUC__
  printf("Compiled " __DATE__ " with GCC " __VERSION__ "\n");
#endif
}

/* This function is called for each configuration file or command line
   option matching an entry in APP_SHORT_OPTS or APP_LONG_OPTS.
   opt is the value that would be returned by getopt_long().
   arg the value that getopt_long() would place in optarg.
   source is the name of the configuration file the option was read from, or
   NULL if the option came from the command line.

   If source is not NULL then arg (if not NULL) points to temporary memory
   and must be copied if needed.

   If opt is zero then arg is a non-option argument.

   This function should not depend on the options being read in any
   particular order, and should just do the minimum to parse the option
   and/or argument and leave the rest to app_init().

   Return 0 if the option is OK, -1 if the argument is invalid, -2 if out of
   range.
*/
int app_parse_option(int opt, char *arg, const char *source)
{
  int status = 0;

  switch (opt)
  {
    case 'b':
      status = parse_uint(&bitsatatime,arg,1,(1U<<31)-1);
      break;
      
    case 'k':
      status = parse_uint64(&kmin,arg,1,(UINT64_C(1)<<62)-1);
      break;

    case 'K':
      status = parse_uint64(&kmax,arg,1,(UINT64_C(1)<<62)-1);
      break;

    case 'n':
      status = parse_uint(&nmin,arg,1,(1U<<31)-1);
      break;

    case 'N':
      status = parse_uint(&nmax,arg,1,(1U<<31)-1);
      break;

    case 'i':
      input_filename = (source == NULL)? arg : xstrdup(arg);
      break;

    case 'f':
      factors_filename = (source == NULL)? arg : xstrdup(arg);
      break;

      /*
    case 's':
    case 'a':
      if(arg[0] == 'y' || arg[0] == 'Y') use_sse2 = 1;
      else if(arg[0] == 'n' || arg[0] == 'N') use_sse2 = 0;
      break;*/
      
    //case 'q':
      //print_factors = 0;
    case 'R':
      search_proth = 0;
      break;
  }

  return status;
}

void app_help(void)
{
  printf("-a --alt=yes|no    Force setting of alt. algorithm (64-bit/SSE2)\n");
  printf("-b --bitsatatime=b Bits to use at a time: fiddle with this, 5-9.\n");
  printf("-f --factors=FILE  Write factors to FILE (default `%s')\n",
         FACTORS_FILENAME_DEFAULT);
  printf("-i --input=FILE    Read initial sieve from FILE\n");
  printf("-k --kmin=K0\n");
  printf("-K --kmax=K1       Sieve for primes k*2^n+/-1 with K0 <= k <= K1\n");
  printf("-n --nmin=N0\n");
  printf("-N --nmax=N1       Sieve for primes k*2^n+/-1 with N0 <= n <= N1\n");
  printf("-R --riesel        Sieve for primes k*2^n-1 instead of +1.\n");
}

// find the log base 2 of a number.  Need not be fast; only done twice.
int lg2(uint64_t v) {
	int r = 0; // r will be lg(v)

	while (v >>= 1) r++;
	return r;
}

/* This function is called once before any threads are started.
 */
void app_init(void)
{
  FILE *file = NULL;
  unsigned int i;

  print_factors = (quiet_opt)?0:1;
  if (input_filename == NULL && (kmin == 0 || kmax == 0))
  {
    fprintf(stderr,"Please specify an input file or both of kmin,kmax\n");
    exit(EXIT_FAILURE);
  }

  if (input_filename != NULL
      && (kmin == 0 || kmax == 0 || nmin == 0 || nmax == 0))
    file = scan_input_file(input_filename);

  if (kmin > kmax)
  {
    fprintf(stderr,"kmin <= kmax is required\n");
    exit(EXIT_FAILURE);
  }

  if (kmax >= pmin)
  {
    fprintf(stderr,"kmax < pmin is required\n");
    exit(EXIT_FAILURE);
  }

  if (kmax-kmin >= (UINT64_C(3)<<36))
  {
    fprintf(stderr,"kmax-kmin < 3*2^36 is required\n");
    exit(EXIT_FAILURE);
  }

  if (nmin == 0)
  {
    // We can calculate nmin = at least 2*log2(pmin)-log2(kmax),
    // because any number smaller than this, divisible by this prime,
    // would also have been divisible by a smaller prime.
    nmin = 2*lg2(pmin)-lg2(kmax)-1;

    //fprintf(stderr,"Please specify a value for nmin\n");
    //exit(EXIT_FAILURE);
  }

  if (nmax == 0)
    nmax = nmin;

  if (nmin > nmax)
  {
    fprintf(stderr,"nmin <= nmax is required\n");
    exit(EXIT_FAILURE);
  }

  b0 = kmin/2;
  b1 = kmax/2;
  kmin = b0*2+1;
  kmax = b1*2+1;

  for (nstep = 1; (kmax << nstep) < pmin; nstep++)
    ;
  ld_nstep = nstep;
  // Calculate the values that fit the given bitsatatime.
  if (nstep > bitsatatime) {
    bpernstep = nstep/bitsatatime;
    nstep = bpernstep*bitsatatime;
  }
  if (nstep > (nmax-nmin+1))
    nstep = (nmax-nmin+1);

  nstart = nmin;

  printf("nstart=%u, nstep=%u\n",nstart,nstep);

  // Allocate and fill bitmap.
  if (input_filename != NULL)
  {
    bitmap = xmalloc((nmax-nmin+1)*sizeof(unsigned char *));
    for (i = nmin; i <= nmax; i++)
    {
      bitmap[i-nmin] = xmalloc((unsigned int)((b1-b0+8)/8));
      memset(bitmap[i-nmin],0,(unsigned int)((b1-b0+8)/8));
    }
    if (file_format == FORMAT_ABCD)
      read_abcd_file(input_filename, file);
    else /* if (file_format == FORMAT_NEWPGEN) */
      read_newpgen_file(input_filename, file);
  }

  if (factors_filename == NULL)
    factors_filename = FACTORS_FILENAME_DEFAULT;
  if ((factors_file = fopen(factors_filename,"a")) == NULL)
  {
    fprintf(stderr,"Cannot open factors file `%s'\n",factors_filename);
    exit(EXIT_FAILURE);
  }

  // Allocate the bitsskip arrays.
  bitsskip = xmalloc(num_threads*sizeof(uint64_t*));
  bitsmask = 1<<bitsatatime; // Not finalized here - off by 1.
  for(i=0; i < num_threads; i++) {
    bitsskip[i] = xmalloc(bitsmask*sizeof(uint64_t));
  }
  bitsmask--; // Finalize bitsmask.

  factor_found = xmalloc(num_threads*sizeof(unsigned char*));

#ifdef _WIN32
  InitializeCriticalSection(&factors_mutex);
#else
  pthread_mutex_init(&factors_mutex,NULL);
#endif

  printf("ppsieve initialized: %"PRIu64" <= k <= %"PRIu64", %u <= n <= %u\n",
         kmin,kmax,nmin,nmax);
  fflush(stdout);
}

/* This function is called once in thread th, 0 <= th < num_threads, before
   the first call to app_thread_fun(th, ...).
 */
unsigned int app_thread_init(int th)
{
  uint16_t mode;
  unsigned int cthread_count;

  cthread_count = cuda_app_init(th);

  // Allocate the factor_found arrays.
  if(cthread_count > 0) factor_found[th] = xmalloc(cthread_count*sizeof(unsigned char));

  /* Set FPU to use extended precision and round to zero. This has to be
     done here rather than in app_init() because _beginthreadex() doesn't
     preserve the FPU mode. */

  asm ("fnstcw %0" : "=m" (mode) );
  mode |= 0x0F00;
  asm volatile ("fldcw %0" : : "m" (mode) );
  return cthread_count;
}

// Fill the bitskip array, to multiply by 2^-b at once.
// O(2^bitsatatime) performance, only 2^bitsatatime writes.
void fillbitskip(uint64_t *bitskip, uint64_t p) {
  int len = 1<<bitsatatime;
  int halflen=len/2; 
  int j, k; 

  // Initialize the first two entries.
  bitskip[halflen] = (p+1)/2;	// Needed first.
  bitskip[0] = 0;			// Ignored through the end.

  // Fill in the intervening spaces, two numbers at a time.
  for(j=halflen; j > 1; j >>= 1) {
    for(k=j/2; k < halflen; k+=j) {
      register uint64_t bl = bitskip[2*k];
      //printf("Filling k=%d from bitskip=%lu\n", k, bl);
      bitskip[k] = (bl+((bl&1)?p:(uint64_t)0))/2;
      //printf("Filling k=%d\n", k+halflen);
      bitskip[k+halflen] = (bl+1+((bl&1)?(uint64_t)0:p))/2;
    }
  }
}

// Old CPU version...
/*
inline void h_check_ns(const uint64_t *__attribute__((aligned(16))) P, const uint64_t *__attribute__((aligned(16))) K, int th) {
  uint64_t *bs0 = bitsskip[th];
  uint64_t n; // = nmin;
  int i;
  for (i = 0; i < cthread_count; i++) {
    uint64_t k0 = K[i];
    factor_found[th][i] = 0;
    if(search_proth) k0 = P[i]-k0;
    // Initialize bitsskip array.
    fillbitskip(bs0, P[i]);
    //if(P[i] == 42070000198537ul) fprintf(stderr, "Checking N's for P=42070000198537\n");
    n = nmin;
    do { // Remaining steps are all of equal size nstep
      //int j;
      uint64_t kpos;
      unsigned int mpos;
      kpos = k0;
      mpos = __builtin_ctzll(kpos);

      kpos >>= mpos;
      if (kpos <= kmax)// && kpos >= kmin && mpos < nstep)
        // Just flag this if kpos <= kmax.
        factor_found[th][i] = 1;
        //test_factor(P[i],kpos,n+mpos,+1);

      switch(bpernstep) {
        case 12:k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 11:k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 10:k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 9: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 8: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 7: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 6: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 5: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 4: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 3: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 2: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
        case 1: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
                break;
 */       /*default:
                for(j=0; j < bpernstep; j++) {
                  k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
                }
                break;*//*
      }
      n += nstep;
    } while (n < nmax);
  }
}*/

inline void test_one_p(const uint64_t P, uint64_t k0, int th) {
  uint64_t n; // = nmin;
  //uint64_t k0 = K;
  // Initialize bitsskip array.
  uint64_t *bs0 = bitsskip[th];
  int cands_found = 0;
  fillbitskip(bs0, P);
  if(search_proth) k0 = P-k0;
  n = nmin;
  do { // Remaining steps are all of equal size nstep
    //int j;
    uint64_t kpos;
    unsigned int mpos;
    kpos = k0;
    mpos = __builtin_ctzll(kpos);

    kpos >>= mpos;
    if (kpos <= kmax) {
      cands_found++;
      if (kpos >= kmin && mpos < nstep)
        test_factor(P,kpos,n+mpos,+1);
    }

    switch(bpernstep) {
      case 12:k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 11:k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 10:k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 9: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 8: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 7: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 6: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 5: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 4: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 3: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 2: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
      case 1: k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
              break;
              /*default:
                for(j=0; j < bpernstep; j++) {
                k0 = (k0 >> bitsatatime) + bs0[(unsigned int)k0 & bitsmask];
                }
                break;*/
    }
    n += nstep;
  } while (n < nmax);
  if(cands_found == 0) {
    fprintf(stderr, "Computation Error: no candidates found for p=%"PRIu64".\n", P);
  }
}
// Given APP_BUFLEN P's, calculate the appropriate K's to start with, based on nmin.
void init_ks(const uint64_t *__attribute__((aligned(16))) P, uint64_t *__attribute__((aligned(16))) K)
{
  uint64_t T[APP_BUFLEN] __attribute__((aligned(16)));
#if (APP_BUFLEN >= 7)
  long double INV[APP_BUFLEN-6];
#endif
#ifdef EMM
  static const uint64_t xones[2] __attribute__((aligned(16))) = {1ul, 1ul};
  __m128i mones = _mm_load_si128((__m128i*)xones);
#endif
#ifdef __x86_64__
  uint64_t x;
#ifndef _WIN32
  unsigned int nbits, t;
#endif
#else
  unsigned int x;
#endif
  unsigned int i, n;

  n = nstart;

#if (APP_BUFLEN <= 6)
  for (i = 0; i < APP_BUFLEN; i++)
  {
    asm volatile ("fildll %0\n\t"
                  "fld1\n\t"
                  "fdivp"
                  : : "m" (P[(APP_BUFLEN-1)-i]) );
    //if(P[i] == 42070000198537ul) fprintf(stderr, "Setting K0 for P=42070000198537\n");
#ifndef __x86_64__
#ifndef __SSE2__
    K[i] = (P[i]+1)/2; /* K[i] <-- 2^{-1} mod P[i] */
#endif
#endif
  }
#else
  for (i = 0; i < 6; i++)
    asm volatile ("fildll %0\n\t"
                  "fld1\n\t"
                  "fdivp"
                  : : "m" (P[5-i]) );

  for ( ; i < APP_BUFLEN; i++)
    asm volatile ("fildll %1\n\t"
                  "fld1\n\t"
                  "fdivp\n\t"
                  "fstpt %0"
                  : "=m" (INV[i-6]) : "m" (P[i]) );

#ifndef __x86_64__
#ifndef __SSE2__
  for (i = 0; i < APP_BUFLEN; i++)
    K[i] = (P[i]+1)/2; /* K[i] <-- 2^{-1} mod P[i] */
#endif
#endif
#endif

#ifdef __SSE2__
  x = 1U << (30 - __builtin_clz(n));
  {
    // Load the P's into SSE2 registers.
    __m128i mp0 = _mm_load_si128((__m128i*)(&P[0]));
    __m128i mp1 = _mm_load_si128((__m128i*)(&P[2]));
    __m128i mp2 = _mm_load_si128((__m128i*)(&P[4]));
    __m128i mtemp, mk;

    // The loop is unrolled, since each P pair is in a different register.
    // The first square is just 2^-1*2^-1 = (2^-1)/2.
    // So do that without a mulmod.
    // The first bit should be checked, too.
    if (n & x) {
      //k0 = (p0+1)/2
      mk = _mm_add_epi64(mp0, mones);
      mk = _mm_srli_epi64(mk, 1);
      for(i=0; i < 2; ++i) {
        //k0 += (k0 % 2)?p0:0;
        mtemp = _mm_and_si128(mk, mones);	// (mk % 2)
        mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
        mtemp = _mm_andnot_si128(mtemp, mp0);	// mp if mk%2 == 1; 0 if mk%2 == 0.
        mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
        //k0 /= 2;
        mk = _mm_srli_epi64(mk, 1);
      }
      _mm_store_si128((__m128i*)(&K[0]), mk);

      //k0 = (p0+1)/2
      mk = _mm_add_epi64(mp1, mones);
      mk = _mm_srli_epi64(mk, 1);
      for(i=0; i < 2; ++i) {
        //k0 += (k0 % 2)?p0:0;
        mtemp = _mm_and_si128(mk, mones);	// (mk % 2)
        mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
        mtemp = _mm_andnot_si128(mtemp, mp1);	// mp if mk%2 == 1; 0 if mk%2 == 0.
        mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
        //k0 /= 2;
        mk = _mm_srli_epi64(mk, 1);
      }
      _mm_store_si128((__m128i*)(&K[2]), mk);

      //k0 = (p0+1)/2
      mk = _mm_add_epi64(mp2, mones);
      mk = _mm_srli_epi64(mk, 1);
      for(i=0; i < 2; ++i) {
        //k0 += (k0 % 2)?p0:0;
        mtemp = _mm_and_si128(mk, mones);	// (mk % 2)
        mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
        mtemp = _mm_andnot_si128(mtemp, mp2);	// mp if mk%2 == 1; 0 if mk%2 == 0.
        mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
        //k0 /= 2;
        mk = _mm_srli_epi64(mk, 1);
      }
      _mm_store_si128((__m128i*)(&K[4]), mk);
    } else {
      // Same thing, but only divides by 2 twice.
      //k0 = (p0+1)/2
      mk = _mm_add_epi64(mp0, mones);
      mk = _mm_srli_epi64(mk, 1);
      //k0 += (k0 % 2)?p0:0;
      mtemp = _mm_and_si128(mk, mones);		// (mk % 2)
      mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
      mtemp = _mm_andnot_si128(mtemp, mp0);	// mp if mk%2 == 1; 0 if mk%2 == 0.
      mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
      //k0 /= 2;
      mk = _mm_srli_epi64(mk, 1);
      _mm_store_si128((__m128i*)(&K[0]), mk);

      //k0 = (p0+1)/2
      mk = _mm_add_epi64(mp1, mones);
      mk = _mm_srli_epi64(mk, 1);
      //k0 += (k0 % 2)?p0:0;
      mtemp = _mm_and_si128(mk, mones);		// (mk % 2)
      mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
      mtemp = _mm_andnot_si128(mtemp, mp1);	// mp if mk%2 == 1; 0 if mk%2 == 0.
      mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
      //k0 /= 2;
      mk = _mm_srli_epi64(mk, 1);
      _mm_store_si128((__m128i*)(&K[2]), mk);

      //k0 = (p0+1)/2
      mk = _mm_add_epi64(mp2, mones);
      mk = _mm_srli_epi64(mk, 1);
      //k0 += (k0 % 2)?p0:0;
      mtemp = _mm_and_si128(mk, mones);		// (mk % 2)
      mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
      mtemp = _mm_andnot_si128(mtemp, mp2);	// mp if mk%2 == 1; 0 if mk%2 == 0.
      mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
      //k0 /= 2;
      mk = _mm_srli_epi64(mk, 1);
      _mm_store_si128((__m128i*)(&K[4]), mk);
    }
  }   // Discard the cached P values in the SSE2 registers.
#else
  // Just initialize it straight up.
  x = 1U << (31 - __builtin_clz(n));
#endif

  while ((x >>= 1) > 0)
  {
    asm ("fildll %1\n\t"
         "fmul %%st(0)\n\t"
         "fmul %%st(1)\n\t"
         "fistpll %0"
         : "=m" (T[0]) : "m" (K[0]) );

#if (APP_BUFLEN >= 2)
    asm ("fildll %1\n\t"
         "fmul %%st(0)\n\t"
         "fmul %%st(2)\n\t"
         "fistpll %0"
         : "=m" (T[1]) : "m" (K[1]) );
#endif
#if (APP_BUFLEN >= 3)
    asm ("fildll %1\n\t"
         "fmul %%st(0)\n\t"
         "fmul %%st(3)\n\t"
         "fistpll %0"
         : "=m" (T[2]) : "m" (K[2]) );
#endif
#if (APP_BUFLEN >= 4)
    asm ("fildll %1\n\t"
         "fmul %%st(0)\n\t"
         "fmul %%st(4)\n\t"
         "fistpll %0"
         : "=m" (T[3]) : "m" (K[3]) );
#endif
#if (APP_BUFLEN >= 5)
    asm ("fildll %1\n\t"
         "fmul %%st(0)\n\t"
         "fmul %%st(5)\n\t"
         "fistpll %0"
         : "=m" (T[4]) : "m" (K[4]) );
#endif
#if (APP_BUFLEN >= 6)
    asm ("fildll %1\n\t"
         "fmul %%st(0)\n\t"
         "fmul %%st(6)\n\t"
         "fistpll %0"
         : "=m" (T[5]) : "m" (K[5]) );
#endif
#if (APP_BUFLEN >= 7)
    for (i = 6; i < APP_BUFLEN; i++)
      asm ("fldt %2\n\t"
           "fildll %1\n\t"
           "fmul %%st(0)\n\t"
           "fmulp\n\t"
           "fistpll %0"
           : "=m" (T[i]) : "m" (K[i]), "m" (INV[i-6]) );
#endif

#ifdef __x86_64__
    /* A correction is required more often as P[i] increases, but no more
       than about 1 time in 8 on average, even for the largest P[i]. */
    for (i = 0; i < APP_BUFLEN; i++)
      if (__builtin_expect(((K[i] = K[i]*K[i] - T[i]*P[i]) >= P[i]),0))
        K[i] -= P[i];

    if (n & x)
    {
      for (i = 0; i < APP_BUFLEN; i++)
      {
        K[i] += (K[i] % 2)? P[i] : 0; /* Unpredictable */
        K[i] /= 2;
      }
    }
#else
#ifdef EMM // 32-bit SSE2 code for K^2-TP
    if (n & x)
    {
      for (i = 0; i < APP_BUFLEN; i+=2) {
        register __m128i mt, mtemp, mtemp2;
        register __m128i mk = _mm_load_si128((__m128i*)(&K[i]));
        register __m128i mp = _mm_load_si128((__m128i*)(&P[i])); // Slip this in the latency.
        // K * K
        mtemp = _mm_srli_epi64(mk, 32);		// Get the high doubleword.
        mtemp = _mm_mul_epu32(mtemp, mk);
        mk = _mm_mul_epu32(mk, mk);
        mt = _mm_load_si128((__m128i*)(&T[i]));	// Slip this in the latency.
        mtemp = _mm_slli_epi64(mtemp, 33);	// Move result to other column, multiply by 2.
        mk = _mm_add_epi32(mk, mtemp);		// Add the results; only need high doublewords.
        // T * P
        mtemp = _mm_srli_epi64(mp, 32);		// Get the high doubleword.
        mtemp2 = _mm_srli_epi64(mt, 32);	// Get the high doubleword.
        mtemp = _mm_mul_epu32(mtemp, mt);
        mtemp2 = _mm_mul_epu32(mtemp2, mp);
        mt = _mm_mul_epu32(mt, mp);
        mtemp = _mm_add_epi32(mtemp, mtemp2);	// Just need the low doublewords.
        mtemp = _mm_slli_epi64(mtemp, 32);	// Move result to other column (high doublewords).
        mt = _mm_add_epi32(mt, mtemp);		// Add the results; only need high doublewords.
        // K*K-T*P
        mk = _mm_sub_epi64(mk, mt);
        // In case of (n & x), do the divide by two here.
        //k0 += (k0 % 2)?p0:0;
        mtemp = _mm_and_si128(mk, mones);	// (mk % 2)
        mtemp = _mm_sub_epi64(mtemp, mones);	// 1 goes to 0; 0 goes to FFFFFFFFFFFFFFFF.
        mtemp = _mm_andnot_si128(mtemp, mp);	// mp if mk%2 == 1; 0 if mk%2 == 0.
        mk = _mm_add_epi64(mk, mtemp);		// mk += the result.
        //k0 /= 2;
        mk = _mm_srli_epi64(mk, 1);
        _mm_store_si128((__m128i*)(&K[i]), mk);
        mk = _mm_sub_epi64(mk, mp);     // Negative iff K[i] < P[i]
        unsigned int bits = _mm_movemask_epi8(mk);
        if((bits & 0x80) == 0) K[i] -= P[i];
        if((bits & 0x8000) == 0) K[i+1] -= P[i+1];
      }
    } else {
      for (i = 0; i < APP_BUFLEN; i+=2) {
        register __m128i mt, mtemp, mtemp2;
        register __m128i mk = _mm_load_si128((__m128i*)(&K[i]));
        register __m128i mp = _mm_load_si128((__m128i*)(&P[i])); // Slip this in the latency.
        mtemp = _mm_srli_epi64(mk, 32);
        mtemp = _mm_mul_epu32(mtemp, mk);
        mk = _mm_mul_epu32(mk, mk);
        mt = _mm_load_si128((__m128i*)(&T[i])); // Slip this in the latency.
        mtemp = _mm_slli_epi64(mtemp, 33); // Move result to other column, multiply by 2.
        mk = _mm_add_epi32(mk, mtemp);  // Add the results; only need high doublewords.
        mtemp = _mm_srli_epi64(mp, 32);
        mtemp2 = _mm_srli_epi64(mt, 32);
        mtemp = _mm_mul_epu32(mtemp, mt);
        mtemp2 = _mm_mul_epu32(mtemp2, mp);
        mt = _mm_mul_epu32(mt, mp);
        mtemp = _mm_add_epi32(mtemp, mtemp2); // Just need the low doublewords.
        mtemp = _mm_slli_epi64(mtemp, 32); // Move result to other column (high doublewords).
        mt = _mm_add_epi32(mt, mtemp);  // Add the results; only need high doublewords.
        mk = _mm_sub_epi64(mk, mt);
        // In case of (n & x), do the divide by two here; This is not that case.
        _mm_store_si128((__m128i*)(&K[i]), mk);
        mk = _mm_sub_epi64(mk, mp);     // Negative iff K[i] < P[i]
        unsigned int bits = _mm_movemask_epi8(mk);
        if((bits & 0x80) == 0) K[i] -= P[i];
        if((bits & 0x8000) == 0) K[i+1] -= P[i+1];
      }      
    }
#else // 32-bit only code for K^2-TP
    for (i = 0; i < APP_BUFLEN; i++)
      if (__builtin_expect(((K[i] = K[i]*K[i] - T[i]*P[i]) >= P[i]),0))
        K[i] -= P[i];

    if (n & x)
    {
      for (i = 0; i < APP_BUFLEN; i++)
        // This should force use of the CMOV instruction.
        // It's faster than a compare when K[i]%2 is random.
        K[i] += (((unsigned int)K[i]) & 1)?P[i]:0;
      // When dealing directly with memory, leave the memory
      // latency time to get assigned before re-reading it.
      for (i = 0; i < APP_BUFLEN; i++)
        K[i] /= 2;
    }
#endif
#endif
  }

#if (APP_BUFLEN <= 6)
  for (i = 0; i < APP_BUFLEN; i++)
    asm volatile ("fstp %st(0)");
#else
  for (i = 0; i < 6; i++)
    asm volatile ("fstp %st(0)");
#endif
}

/* This function is called 0 or more times in thread th, 0 <= th < num_threads.
   P is an array of APP_BUFLEN candidate primes.
*/
void app_thread_fun(int th, const uint64_t *__attribute__((aligned(16))) P, uint64_t *__attribute__((aligned(16))) K, unsigned int cthread_count)
{
  unsigned int i;
  // Set up tables on the GPU.
  setup_ps(P, cthread_count);
  // Initialize all K's from P's.
  for(i=0; i < cthread_count-APP_BUFLEN; i+=APP_BUFLEN) {
    init_ks(&P[i], &K[i]);
#ifndef NDEBUG
    //fprintf(stderr, "Inited K's through %d\n", i+APP_BUFLEN-1);
#endif
  }
  // Do the final ones, possibly duplicating work.
  init_ks(&P[cthread_count-APP_BUFLEN], &K[cthread_count-APP_BUFLEN]);

  check_ns(P, K, factor_found[th], cthread_count);
  for(i=0; i < cthread_count; i++) {
    if(factor_found[th][i]) {
      test_one_p(P[i], K[i], th);
    }
  }
}

/* This function is called 0 or more times in thread th, 0 <= th < num_threads.
   P is an array of len candidate primes, 0 <= len < APP_PRIMES_BUFLEN.
   The application should be prepared to checkpoint after returning from
   this function.
*/
void app_thread_fun1(int th, uint64_t *P, uint64_t *K, unsigned int cthread_count, unsigned int len)
{
  unsigned int i;

  if (len > 0)
  {
    i = len;
    // Find a non-prime and stick in in the rest of the locations.
    if(len < cthread_count) {
      P[i] = P[i-1];
      while(P[i] % 3 != 0) P[i] += 2;
      for (i++ ; i < cthread_count; i++)
        P[i] = P[i-1];
    }

    app_thread_fun(th,P,K, cthread_count);
  }
}

/* This function is called once in thread th, 0 <= th < num_threads, after
   the final call to app_thread_fun1(th, ...).
*/
void app_thread_fini(int th)
{
}

/* This function is called at most once, after app_init() but before any
   threads have started, with fin open for reading.
   Return 0 if the checkpoint data is invalid, which will cause this
   checkpoint to be ignored.
*/
int app_read_checkpoint(FILE *fin)
{
  unsigned int n0, n1;

  if (fscanf(fin,"nmin=%u,nmax=%u,factor_count=%u",&n0,&n1,&factor_count) != 3)
    return 0;

  if (n0 != nmin || n1 != nmax)
    return 0;

  return 1;
}

/* This function is called periodically with fout open for writing. The
   application can assume that the following conditions apply:

   Threads are blocked for the duration of this call, or have already
   exited. (And so app_thread_fini(th) may have been called already).

   app_thread_fun(th, ...) has not been called since the last call to
   app_thread_fun1(th, ...).

   All candidates before, and no candidates after the checkpoint have been
   passed to one of the functions app_thread_fun() or app_thread_fun1().
*/
void app_write_checkpoint(FILE *fout)
{
  fflush(factors_file);
  fprintf(fout,"nmin=%u,nmax=%u,factor_count=%u\n",nmin,nmax,factor_count);
}

/* This function is called once after all threads have exited.
 */
void app_fini(void)
{
  unsigned int i;

  fclose(factors_file);
  printf("Found %u factor%s\n",factor_count,(factor_count==1)? "":"s");

  cuda_finalize();

#ifdef _WIN32
  DeleteCriticalSection(&factors_mutex);
#else
  pthread_mutex_destroy(&factors_mutex);
#endif

  if (bitmap != NULL)
  {
    for (i = nmin; i <= nmax; i++)
      free(bitmap[i-nmin]);
    free(bitmap);
    bitmap = NULL;
  }
  for(i=0; i < num_threads; i++) {
    free(bitsskip[i]);
  }
  free(bitsskip);
}
