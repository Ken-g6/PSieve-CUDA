/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.c -- (C) Geoffrey Reynolds, April-August 2009.
 * With improvements by Ken Brazier August-October 2009.

   Proth Prime Search sieve (for many K).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

//#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "main.h"
#ifdef _WIN32
#include <windows.h>
#include "getopt.h"
#include "stdint.h"
#include "inttypes.h"
#else
#include <getopt.h>
#include <stdint.h>
#include <inttypes.h>
#ifndef SINGLE_THREAD
#include <pthread.h>
#endif
#endif
#ifdef __SSE2__
#define EMM
#include <emmintrin.h>
#ifdef __x86_64__
#include <time.h>
#endif
#endif
#include "putil.h"
#include "app.h"
#include "appcu.h"
#include "clock.h"
#include "factor_proth.h"
#ifdef __GCC__
#define INLINE static inline
#else
#define INLINE static __inline
#endif

// Macro bsfq (Bit Search Forward Quadword) for 32-bit.
// MPOS = result (32-bit)
// KPOS = input (64-bit; evaluated twice!)
// ID = a unique ID string.
#if !defined(__x86_64__) && defined(__i386__) && defined(__GNUC__)
#define BSFQ(MPOS,KPOS,ID)          asm volatile \
            ( 		"bsfl	%[k0l], %[pos]	\n" \
                "	jnz	bsflok" #ID "		\n" \
                "	bsfl	%[k0h], %[pos]	\n" \
                "	addl	$32, %[pos]	\n" \
                "bsflok" #ID ":" \
                : [pos] "=r" (MPOS) \
                : [k0l] "rm" ((unsigned int)(KPOS)), \
                [k0h] "rm" ((unsigned int)((KPOS) >> 32)) \
                : "cc" )
#else
// If anyone wants to compile on some non-x86 platform somehow...
#define BSFQ(MPOS,KPOS,ID) MPOS = __builtin_ctzll(KPOS)
#endif

#define FORMAT_NEWPGEN 1
#define FORMAT_ABCD 2

uint64_t kmin = 1, kmax = 0;
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
int search_proth = 1; // Search for Proth or Riesel numbers?
int *check_ns_delay;
//static unsigned int bitsatatime = 8; // Bits to process at a time, with v0.4 algorithm.
//static unsigned int bitsmask, bpernstep;
static uint64_t* gpu_started;
//static uint64_t** bitsskip;
static unsigned char** factor_found;
static int device_opt = -1;
static unsigned int user_cthread_count = 0;
static uint64_t r0arr[9];
static int bbitsarr[9];
static int b_term=0; // 0 or 2.
static int constant_c = 0;  // -1, 1, or 0 for non-constant.

struct n_subsection *thread_subsections;

#ifndef SINGLE_THREAD
#ifdef _WIN32
static CRITICAL_SECTION factors_mutex;
#else
static pthread_mutex_t factors_mutex;
#endif
#endif


static void report_factor(uint64_t p, uint64_t k, unsigned int n, int c)
{
#ifndef SINGLE_THREAD
#ifdef _WIN32
  EnterCriticalSection(&factors_mutex);
#else
  pthread_mutex_lock(&factors_mutex);
#endif
#endif

  if (factors_file != NULL && fprintf(factors_file,"%"PRIu64" | %"PRIu64"*2^%u%+d\n",p,k,n,c) > 0)
  {
    if(print_factors) printf("%"PRIu64" | %"PRIu64"*2^%u%+d\n",p,k,n,c);
  }
  else fprintf(stderr, "%sUNSAVED: %"PRIu64" | %"PRIu64"*2^%u%+d\n",bmprefix(),p,k,n,c);
  factor_count++;

#ifndef SINGLE_THREAD
#ifdef _WIN32
  LeaveCriticalSection(&factors_mutex);
#else
  pthread_mutex_unlock(&factors_mutex);
#endif
#endif
}

// 1 if a number mod 15 is not divisible by 5 or 3.
//                             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
static const int prime15[] = { 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1 };

#ifdef _WIN32
static void
#else
static void __attribute__((noinline))
#endif
test_factor(uint64_t p, const uint64_t k_lg, unsigned int n, int c)
{
  // K's cannot exceed N's here, which cannot be uint64_t's.
  unsigned int k = (unsigned int)k_lg;

  // Normalize the Cullen/Woodall candidate.
  while(k<n) {
    k <<= 1;
    n--;
  }

  if((k == n) && k >= nmin && k <= nmax) {
    if (bitmap == NULL) {
      //uint64_t mod31 = (uint64_t)1;
      // Check that K*2^N+/-1 is not divisible by 3, 5, or 7, to minimize factors printed.
      // We do 3 and 5 at the same time (15 = 2^4-1), then 7 (=2^3-1).
      // Then 17, 11 (and 31), 13, and maybe 19, if there's space. 23 can also go in there, if it's worth it.
      // (k*(1<<(n%2))+c)%3 == 0
      if(prime15[(unsigned int)(((k<<(n&3))+c)%15)] && 
          (unsigned int)(((k<<(n%3))+c)%7) != 0)
        if(try_all_factors(k, n, c) == 0)
          report_factor(p,k,n,c);
    } else {
      //printf("%sUNTESTED: %"PRIu64" | %u*2^%u%+d\n",bmprefix(),p,k,n,c);
      if (bitmap[(c==-1)?0:1][(unsigned int)((k-b0)/8)] & (1<<(k-b0)%8))
        report_factor(p,k,n,c);
    }
  }
}

// find the log base 2 of a number.  Need not be fast; only done twice.
int lg2(uint64_t v) {
  int r = 0; // r will be lg(v)

  while (v >>= 1) r++;
  return r;
}
static void line_error(const char *msg, const char *file_name)
{
  fprintf(stderr, "Error in %s file: %s.",file_name, msg);
}

#define ABCDFORMATSTR "ABC $a*$b^$a$%c"
/* Scan the input file to determine format and values for kmin,kmax,nmin,nmax.
 */
#define MAX_LINE_LENGTH 256
#define COMMENT_CHAR '#'
const char *read_line(FILE *file)
{
  static char buf[MAX_LINE_LENGTH];
  const char *ptr;

  while ((ptr = fgets(buf,sizeof(buf),file)) != NULL)
  {
    while (isspace(*ptr))
      ptr++;
    if (*ptr != COMMENT_CHAR && *ptr != '\0')
      break;
  }

  return ptr;
}

static void check_p(uint64_t p, const char *file_name)
{
  if (pmin == 0)
    pmin = p;
  else if (pmin != p) {
    fprintf(stderr, "--pmin=%"PRIu64" from command line overrides pmin=%"PRIu64
            " from `%s'\n", pmin, p, file_name);
  }
}
static void check_c(int c, const char *file_name)
{
  if (c != 1 && c != -1)
    line_error("|c| != 1 in n*b^n+c",file_name);
}

static unsigned int read_abc_header(FILE *file)
{
  const char *file_name = "input";
  const char *line;
  uint64_t p;
  unsigned int b, n;
  int c;
  char ch1, ch2, ch3, ch4;

  /* Try to read header */
  line = read_line(file);

  switch (sscanf(line,"ABC $%c*$%c^$%c$%c //%*[^0-9]%"SCNu64,
                 &ch1,&ch2,&ch3,&ch4,&p))
  {
    default:
      line_error("Invalid ABC header",file_name);

    case 5:
      check_p(p,file_name);
      /* fall through */

    case 4: /* MultiSieve format */
      if (ch1 != 'a' || ch2 != 'b' || ch3 != 'a' || ch4 != 'c')
        line_error("Wrong ABC format",file_name);
        line = read_line(file);
        if(line == NULL) line_error("No data in file", file_name);
        if (sscanf(line, "%u %u %d",&n,&b,&c) != 3)
          line_error("Malformed line",file_name);
        if (b != 2) line_error("mismatched b in n*b^n+c",file_name);
        check_c(c,file_name);
      break;

    case 1:
      switch (sscanf(line,"ABC $%c*%u^$%c%d //%*[^0-9]%"SCNu64,
                 &ch1,&b,&ch3,&c,&p))
      {
        default:
          line_error("Invalid ABC header",file_name);

        case 5:
          check_p(p,file_name);
          /* fall through */

        case 4: /* Fixed c format */
          if (ch1 != 'a' || ch3 != 'a')
            line_error("Wrong ABC format",file_name);
          if (b != 2)
            line_error("Invalid b (not two) in n*b^n+c",file_name);
          b_term = b;
          check_c(c,file_name);
          constant_c = c;
          line = read_line(file);
          if(line == NULL) line_error("No data in file", file_name);
          if (sscanf(line, "%u",&n) != 1)
              line_error("Malformed line",file_name);
          break;

        case 3:
          /* The "%*[+$]" allows the "+" to be ignored if present. This only
             works because we have already determined that a number was not
             matched in this position. */
          switch (sscanf(line,"ABC $%c*%u^$%c%*[+$]%c //%*[^0-9]%"
                         SCNu64,&ch1,&b,&ch3,&ch4,&p))
          {
            default:
              line_error("Invalid ABC header",file_name);

            case 5:
              check_p(p,file_name);
              /* fall through */

            case 4: /* Variable c format */
              if (ch1 != 'a' || ch3 != 'a' || ch4 != 'b')
                line_error("Wrong ABC format",file_name);
              if (b != 2)
                line_error("Invalid b (not two) in n*b^n+c",file_name);
              b_term = b;
              line = read_line(file);
              if(line == NULL) line_error("No data in file", file_name);
              if (sscanf(line, "%u %d",&n,&c) != 2)
                line_error("Malformed line",file_name);
              check_c(c,file_name);
              break;
          }
      }
      break;
  }
  return n;
}

static FILE* scan_input_file(const char *fn)
{
  FILE *file;
  unsigned int n0, n1, n;
  int fileaddsign = 1;
  char ch;

  if ((file = bfopen(fn,"r")) == NULL)
  {
    perror(fn);
    bexit(ERR_FOPEN);
  }
  n = read_abc_header(file);

  ch = 'c';

  n0 = n1 = n;
  while (1)
  {
    while(1) {
      // Read N.
      char c = getc(file);
      if(!isdigit(c)) {
        ungetc(c, file);
        break;
      }
      n = c-'0';
      while(isdigit(c=getc(file))) {
        n *= 10;
        n += c-'0';
      }
      // Skip anything else.
      while(c != '\n' && c != '\r') c=getc(file);
      while(c == '\n' || c == '\r') c=getc(file);
      ungetc(c, file);
    }
    //while (fscanf(file," %u",&d) == 1)

    if (n1 < n) {
      printf("Read maximum n %u\n", n);
      n1 = n;
    }
    else
      break;
  }


  if (ferror(file))
  {
    fprintf(stderr,"%sError reading input file `%s'\n",bmprefix(),fn);
    bexit(ERR_SCANF);
  }

  rewind(file);
  printf("Found N's from %u to %u.\n", n0, n1);

  if (kmax == 0 || kmax > n1)
    kmax = n1;

  if (nmin < n0)
    nmin = n0;
  if (nmax == 0 || nmax > n1)
    nmax = n1 + lg2(n1)+1;
  return file;
}

static void read_abcd_file(const char *fn, FILE *file)
{
  unsigned int n, count;
  char sign;

  if(file == NULL) {
    printf("Opening file %s\n", fn);
    if ((file = bfopen(fn,"r")) == NULL)
    {
      perror(fn);
      bexit(ERR_FOPEN);
    }
  }
  // Skip the first line and any comments before it.
  read_line(file);

  count = 0;
  // Pre-set the sign; it's overwritten if not constant.
  sign = (constant_c == -1)?'-':'+';
  //while(getc(file) != '\n');
  printf("Reading ABC file.\n");
  while(1)
  {
    unsigned int bit;
    char sign;
    char c = getc(file);
    while(c == 10) c = getc(file);
    if(!isdigit(c)) {
      //ungetc(c, file);
      //printf("Read char %d, which is not a digit.\n", c);
      break;
    }
    n = c-'0';
    while(isdigit(c=getc(file))) {
      n *= 10;
      n += c-'0';
    }
    // Read a space before anything else.
    if(b_term == 0 || constant_c == 0) {
      if(c != ' ') { fprintf(stderr, "Invalid file format, line %u\n", count+1); bexit(ERR_SCANF); }
      c = getc(file);
    }

    if(b_term == 0) {
      // Read the base if expected.
      if(c != '2') { fprintf(stderr, "Invalid base (not 2), line %u\n", count+1); bexit(ERR_SCANF); }
      c = getc(file);
      if(c != ' ' && (constant_c == 0 || c != '\r' && c != '\n')) { fprintf(stderr, "Invalid base (not 2), line %u\n", count+1); bexit(ERR_SCANF); }
    }
    if(constant_c == 0) {
      // Read the sign if expected.
      if(c == ' ') c = getc(file);
      sign = c;
    }

    bit = n-nmin;
    if(n >= nmin && n <= nmax) bitmap[(sign=='-')?0:1][bit/8] |= (1 << bit%8);
    //printf("Read %lu*2^%d+/-1\n", k, n);
    count++;
    while(c != '\n' && c != '\r') c=getc(file);
    while(c == '\n' || c == '\r') c=getc(file);
    ungetc(c, file);
  }

  if (ferror(file))
  {
    printf("\nError reading input file `%s'\n",fn);
    bexit(ERR_SCANF);
  }
  //printf("\n\nDone reading ABC file!\n");

  fclose(file);

  printf("Read %u terms from ABC format input file `%s'\n",count,fn);
}

/* This function is called once before anything is done. It should at least
   print the name of the application. It can also assign default option
   values which will be overridden in the configuration file or on the
   command line.
*/
void app_banner(void)
{
  printf(APP_NAME " version " APP_VERSION " (testing)\n");
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
    /*
    case 'b':
      status = parse_uint(&bitsatatime,arg,1,(1U<<31)-1);
      break;
      */
      
    case 'k':
      status = parse_uint64(&kmin,arg,1,(UINT64_C(1)<<62)-1);
      break;

    case 'K':
      status = parse_uint64(&kmax,arg,1,(UINT64_C(1)<<62)-1);
      break;

    case 'm':
      status = parse_uint(&user_cthread_count,arg,1,(1U<<31));
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

    case 'd':
      status = parse_uint((unsigned int *)(&device_opt),arg,0,INT32_MAX);
      break;
      
#ifndef SEARCH_TWIN
    case 'R':
      search_proth = -1;
      break;
#endif
    //case 'q':
      //print_factors = 0;
      //break;
  }

  return status;
}

void app_help(void)
{
  //printf("-b --bitsatatime=b Bits to use at a time: fiddle with this, 5-9.\n");
  printf("-f --factors=FILE  Write factors to FILE (default `%s')\n",
         FACTORS_FILENAME_DEFAULT);
  printf("-i --input=FILE    Read initial sieve from FILE\n");
  printf("-k --kmin=K0\n");
  printf("-K --kmax=K1       Sieve for primes k*2^n+/-1 with K0 <= k <= K1\n");
  printf("-m --mthreads=M    Force M threads or blocks/multiprocessor.\n");
  printf("-n --nmin=N0\n");
  printf("-N --nmax=N1       Sieve for primes k*2^n+/-1 with N0 <= n <= N1\n");
#ifndef SEARCH_TWIN
  printf("-R --riesel        Sieve for primes k*2^n-1 instead of +1.\n");
#endif
  printf("-d --device=N      Use GPU N instead of 0-threads\n");
}

/* This function is called once before any threads are started.
 */
void app_init(void)
{
  FILE *file = NULL;
  unsigned int i;

  print_factors = (quiet_opt)?0:1;
  if (input_filename == NULL && (nmax == 0))
  {
    bmsg("Please specify an input file or nmax\n");
    bexit(ERR_INVALID_PARAM);
  }
  kmax = nmax;

  if (input_filename != NULL
      && (nmax == 0))
    file = scan_input_file(input_filename);

  if (kmin > kmax)
  {
    bmsg("kmin <= kmax is required\n");
    bexit(ERR_INVALID_PARAM);
  }

  if (kmax >= pmin)
  {
    bmsg("kmax < pmin is required\n");
    bexit(ERR_INVALID_PARAM);
  }

  if (nmin == 0)
  {
    // k*2^n is prime if k*2^n < p^2
    // 2^n < p^2/k
    // We can calculate nmin = at least 2*log2(pmin)-log2(kmax),
    // because any number smaller than this, divisible by this prime,
    // would also have been divisible by a smaller prime.
    nmin = 2*lg2(pmin)-lg2(kmax)-1;

    //bmsg("Please specify a value for nmin\n");
    //bexit(ERR_INVALID_PARAM);
  }

  if (nmax == 0)
    nmax = nmin;

  if (nmin > nmax)
  {
    bmsg("nmin <= nmax is required\n");
    bexit(ERR_INVALID_PARAM);
  }

  b0 = nmin;
  b1 = nmax;

  for (nstep = 1; (kmax << nstep) < pmin; nstep++)
    ;
  ld_nstep = nstep;
  // Calculate the values that fit the given bitsatatime.
  /*
  if (nstep > bitsatatime) {
    bpernstep = nstep/bitsatatime;
    nstep = bpernstep*bitsatatime;
  }
  if (nstep > (nmax-nmin+1))
    nstep = (nmax-nmin+1);
  */
  nstart = nmin;

  printf("nstart=%u, nstep=%u\n",nstart,ld_nstep);

  // Allocate and fill bitmap.
  if (input_filename != NULL)
  {
    bitmap = xmalloc(2*sizeof(unsigned char *));
    for (i = 0; i <= 1; i++)
    {
      bitmap[i] = xmalloc((unsigned int)((b1-b0+8)/8));
      memset(bitmap[i],0,(unsigned int)((b1-b0+8)/8));
    }
    if (file_format == FORMAT_ABCD)
      read_abcd_file(input_filename, file);
    else /* if (file_format == FORMAT_NEWPGEN) */
    {
      fprintf(stderr, "NewPGen file format not available.\n");
      bexit(EXIT_FAILURE);
      //read_newpgen_file(input_filename, file);
    }
  } else sieve_small_primes(11);

  if (factors_filename == NULL)
    factors_filename = FACTORS_FILENAME_DEFAULT;
  if ((factors_file = bfopen(factors_filename,"a")) == NULL)
  {
    fprintf(stderr,"%sCannot open factors file `%s'\n",bmprefix(),factors_filename);
    bexit(ERR_FOPEN);
  }

  // Allocate the bitsskip arrays.
  /*
  bitsskip = xmalloc(num_threads*sizeof(uint64_t*));
  bitsmask = 1<<bitsatatime; // Not finalized here - off by 1.
  for(i=0; i < num_threads; i++) {
    bitsskip[i] = xmalloc(bitsmask*sizeof(uint64_t));
  }
  bitsmask--; // Finalize bitsmask.
  */
  // Allocate check_ns_delay;
  check_ns_delay = xmalloc(num_threads*sizeof(int));
  for(i=0; i < (unsigned int)num_threads; i++) check_ns_delay[i] = 0;

  factor_found = xmalloc(num_threads*sizeof(unsigned char*));
  gpu_started = xmalloc(num_threads*sizeof(uint64_t));
  for(i=0; i < (unsigned int)num_threads; i++) {
    gpu_started[i] = (uint64_t)0;
  }

  // Allocate the N subsections for each thread.
  thread_subsections = xmalloc(num_threads*sizeof(struct n_subsection));

  cuda_init();

#ifndef SINGLE_THREAD
#ifdef _WIN32
  InitializeCriticalSection(&factors_mutex);
#else
  pthread_mutex_init(&factors_mutex,NULL);
#endif
#endif

  printf(APP_NAME " initialized: %u <= n <= %u\n", nmin,nmax);
  fflush(stdout);
}

// Replacing the getter in appcu.cu with a macro.
#define get_n_subsection_start(X) thread_subsections[th].n_subsection_start[X]
/* This function is called once in thread th, 0 <= th < num_threads, before
   the first call to app_thread_fun(th, ...).
 */
unsigned int app_thread_init(int th)
{
  unsigned int i, cthread_count;

  if(device_opt >= 0) {
    cthread_count = cuda_app_init(device_opt, th, user_cthread_count);
  } else {
    cthread_count = cuda_app_init(th, th, user_cthread_count);
  }
  // Create r0arr (which gives starting values for the REDC code.)
  //printf("ld_r0[%d] = %lu\n", nmin, ld_r0);
  for(i=0; i < 8; i++) {
    int l_bbits;
    unsigned int n = get_n_subsection_start(i+1);
    l_bbits = lg2(n);

    bbitsarr[i] = l_bbits - 6;
    r0arr[i] = ((uint64_t)1) << (64-(n >> (l_bbits-5)));
    //printf("r0arr[%d] = %lu\n", n, r0arr[i]);
  }

  // Allocate the factor_found arrays.
  if(cthread_count > 0) {
    factor_found[th] = xmalloc(cthread_count*sizeof(unsigned char));
    for(i=0; i < cthread_count; i++) {
      factor_found[th][i] = 0;
    }
  } else factor_found[th] = NULL;

  return cthread_count;
}

/*  Multiplies for REDC code  */

#if defined(GCC) && defined(__x86_64__)
static uint64_t __umul64hi(const uint64_t a, const uint64_t b)
{
  uint64_t t1, t2;
  __asm__
  ( "mulq %3\n\t"
    : "=a" (t1), "=d" (t2)
    : "0" (a), "rm" (b)
    : "cc");
  return t2;
}
#else
#if defined(GCC) && !defined(__x86_64__)
static unsigned int __umulhi(const unsigned int a, const unsigned int b)
{
  unsigned int t1, t2;
  __asm__
  ( "mull %3\n\t"
    : "=a" (t1), "=d" (t2)
    : "0" (a), "rm" (b)
    : "cc");
  return t2;
}
static uint64_t __umul64(const unsigned int a, const unsigned int b)
{
  unsigned int t1, t2;
  __asm__
  ( "mull %3\n\t"
    : "=a" (t1), "=d" (t2)
    : "0" (a), "rm" (b)
    : "cc");
  return (((uint64_t)t2)<<32)+t1;
}
#else
static unsigned int __umulhi(const unsigned int a, const unsigned int b)
{
  uint64_t c = (uint64_t)a * (uint64_t)b;

  return (unsigned int)(c >> 32);
}
static uint64_t __umul64(const unsigned int a, const unsigned int b)
{
  return (uint64_t)a * (uint64_t)b;
}
#endif

static uint64_t __umul64hi(const uint64_t a, const uint64_t b)
{
  const unsigned int a_lo = (unsigned int)a;
  const unsigned int a_hi = (unsigned int)(a >> 32);
  const unsigned int b_lo = (unsigned int)b;
  const unsigned int b_hi = (unsigned int)(b >> 32);
  uint64_t m1 = __umul64(a_lo, b_hi);
  uint64_t m2 = __umul64(a_hi, b_lo);
  unsigned int           carry;

  carry = (((uint64_t)__umulhi(a_lo, b_lo)) + (unsigned int)m1 + (unsigned int)m2) >> 32;

  return __umul64(a_hi, b_hi) + (m1 >> 32) + (m2 >> 32) + carry;
}
#endif

/*  BEGIN REDC CODE  */

static uint64_t invmod2pow_ul (const uint64_t n)
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

static uint64_t mulmod_REDC (const uint64_t a, const uint64_t b, 
             const uint64_t N, const uint64_t Ns)
{
  uint64_t rax, rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  rax = a*b;
  rcx = __umul64hi(a,b);
  rax *= Ns;
  rcx += (rax!=0)?1:0;
  rax = __umul64hi(rax, N);
  rax += rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, b, N))
  {
    fprintf (stderr, "%sError, mulredc(%lu,%lu,%lu) = %lu\n", bmprefix(), a, b, N, rax);
    bexit(ERR_NEG);
  }
#endif

  return rax;
}

// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
static uint64_t onemod_REDC(const uint64_t N, uint64_t rax) {
  uint64_t rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  rcx = (rax!=0)?1:0;
  rax = __umul64hi(rax, N) + rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  return rax;
}

// Like mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
static uint64_t mod_REDC(const uint64_t a, const uint64_t N, const uint64_t Ns) {
  return onemod_REDC(N, Ns*a);
}

// Compute T=a<<s; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
// rax is passed in as a * Ns.
static uint64_t shiftmod_REDC (const uint64_t a, 
             const uint64_t N, uint64_t rax)
{
  uint64_t rcx;
  unsigned int d_mont_nstep = 64-ld_nstep;

  rax <<= d_mont_nstep; // So this is a*Ns*(1<<s) == (a<<s)*Ns.
  rcx = a >> ld_nstep;
  rcx += (rax!=0)?1:0;
  rax = __umul64hi(rax, N) + rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, ((uint64_t)1)<<d_mont_nstep, N))
  {
    fprintf (stderr, "%sError, shiftredc(%lu,%u,%lu) = %lu\n", bmprefix(), a, d_mont_nstep, N, rax);
    bexit(ERR_NEG);
  }
#endif

  return rax;
}
// Same function, but for 32 bits or less nstep.
static uint64_t shiftmod_REDC_sm (uint64_t rcx, 
             const uint64_t N, unsigned int rax)
{
  unsigned int d_mont_nstep = 32-ld_nstep;

  rax <<= d_mont_nstep; // So this is (a*Ns*(1<<s) == (a<<s)*Ns) >> 32.
  rcx >>= ld_nstep;
  rcx += __umulhi(rax, (unsigned int)N) + ((rax!=0)?1:0);
  rcx += __umul64(rax, (unsigned int)(N>>32));
  //rcx = rax - N;
  rcx = (rcx>N)?(rcx-N):rcx;

  /*
#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, ((uint64_t)1)<<d_mont_nstep, N))
  {
    fprintf (stderr, "%sError, shiftredc(%lu,%u,%lu) = %lu\n", bmprefix(), a, d_mont_nstep, N, rax);
    bexit(ERR_NEG);
  }
#endif
*/

  return rcx;
}

// Hybrid powmod, sidestepping several loops and possible mispredicts, and with no longmod!
/* Compute (2^-1)^b (mod m), using Montgomery arithmetic. */
static uint64_t invpowmod_REDClr (const uint64_t N, const uint64_t Ns, const unsigned int l_nmin, uint64_t r, int bbits) {
  // Now work through the other bits of nmin.
  for(; bbits >= 0; --bbits) {
    //if(N == 42070000070587) printf("r = %lu (CPU)\n", r);
    // Just keep squaring r.
    r = mulmod_REDC(r, r, N, Ns);
    // If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod N).
    if(l_nmin & (1u << bbits)) {
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

void test_one_p(const uint64_t my_P, const unsigned int l_nmin, const unsigned int l_nmax, const uint64_t r0, const int l_bbits) {
  unsigned int n = l_nmin; // = nmin;
  unsigned int i;
  uint64_t k0, kPs;
  uint64_t kpos;
//#ifdef SEARCH_TWIN
  //uint64_t kneg;
//#endif
  uint64_t Ps;
  int cands_found = 0;

  //printf("I think it found %lu divides N in {%u, %u}\n", my_P, l_nmin, l_nmax);
  
  // Better get this done before the first mulmod.
  Ps = -invmod2pow_ul (my_P); /* Ns = -N^{-1} % 2^64 */
  
  // Calculate k0, in Montgomery form.
  k0 = invpowmod_REDClr(my_P, Ps, l_nmin, r0, l_bbits);
  //printf("k0[%u] = %lu\n", l_nmin, k0);

  //if(my_P == 42070000070587) printf("%lu^-1 = %lu (CPU)\n", my_P, Ps);
  /*
  // Verify the first result.
  kpos = 1;
  for(i=0; i < l_nmin; i++) {
      kpos += (kpos&1)?my_P:0;
      kpos >>= 1;
  }
  kPs = k0 * Ps;
  assert(kpos == onemod_REDC(my_P, kPs));
  if(kpos != onemod_REDC(my_P, kPs)) {
    fprintf(stderr, "Error: %lu != %lu!\n", kpos, onemod_REDC(my_P, kPs));
    //bexit(ERR_NEG);
  } */
#ifdef SEARCH_TWIN
    // Select the first odd one.  All others are tested by overlap.
    kpos = (k0&1)?k0:(my_P - k0);
    if (kpos <= l_nmin) {
      cands_found++;
      test_factor(my_P,kpos,n,(kpos==k0)?-1:1);
    }
#endif

#ifndef SEARCH_TWIN
  if(search_proth == 1) k0 = my_P-k0;
#endif

  do { // Remaining steps are all of equal size nstep
    // Get K from the Montgomery form.
    // This is equivalent to mod_REDC(k, my_P, Ps), but the intermediate kPs value is kept for later.
    kPs = __umul64(k0, Ps);
    kpos = k0;
#ifdef SEARCH_TWIN
    // Select the even one.
    kpos = (kpos&1)?(my_P - kpos):kpos;
#endif
    //i = __ffsll(kpos)-1;
    //i = __builtin_ctzll(kpos);
    BSFQ(i, kpos, 1);

#ifdef SEARCH_TWIN
    if ((kpos>>i) <= n+i && i <= ld_nstep) {
#else
    if ((kpos>>i) <= n+i && i < ld_nstep) {
#endif
      unsigned int new_k = (kpos>>i);
      unsigned int new_n = n+i;
      while(new_k < new_n) {
        new_k <<= 1;
        new_n--;
      }
      if(new_k == new_n) {
        cands_found++;
        //if (i < ld_nstep)
#ifdef SEARCH_TWIN
        test_factor(my_P,(kpos>>i),n+i,(kpos==k0)?-1:1);
#else
        test_factor(my_P,(kpos>>i),n+i,search_proth);
#endif
      }
    }

    // Proceed to the K for the next N.
    n += ld_nstep;
    k0 = shiftmod_REDC(k0, my_P, kPs);
  } while (n < l_nmax);
  if(cands_found == 0) {
    fprintf(stderr, "%sComputation Error: no candidates found for p=%"PRIu64" in %u<=n<%u.\n", bmprefix(), my_P, l_nmin, l_nmax);
#ifdef USE_BOINC
    bexit(ERR_NEG);
#endif
  }
}

void test_one_p_sm(const uint64_t my_P, const unsigned int l_nmin, const unsigned int l_nmax, const uint64_t r0, const int l_bbits) {
  unsigned int n = l_nmin; // = nmin;
  unsigned int i;
  uint64_t k0; //, kPs;
  uint64_t kpos;
//#ifdef SEARCH_TWIN
  //uint64_t kneg;
//#endif
  uint64_t Ps;
  int cands_found = 0;

  //printf("I think it found %lu divides N in {%u, %u}\n", my_P, l_nmin, l_nmax);
  
  // Better get this done before the first mulmod.
  Ps = -invmod2pow_ul (my_P); /* Ns = -N^{-1} % 2^64 */
  
  // Calculate k0, in Montgomery form.
  k0 = invpowmod_REDClr(my_P, Ps, l_nmin, r0, l_bbits);
  //printf("k0[%u] = %lu\n", l_nmin, k0);

  //if(my_P == 42070000070587) printf("%lu^-1 = %lu (CPU)\n", my_P, Ps);
#ifdef SEARCH_TWIN
    // Select the first odd one.  All others are tested by overlap.
    kpos = (k0&1)?k0:(my_P - k0);
    if (kpos <= l_nmin) {
      cands_found++;
      test_factor(my_P,kpos,n,(kpos==k0)?-1:1);
    }
#endif

#ifndef SEARCH_TWIN
  if(search_proth == 1) k0 = my_P-k0;
#endif

  do { // Remaining steps are all of equal size nstep
    // Get K from the Montgomery form.
    // This is equivalent to mod_REDC(k, my_P, Ps), but the intermediate kPs value is kept for later.
    kpos = k0;
#ifdef SEARCH_TWIN
    // Select the even one.
    kpos = (kpos&1)?(my_P - kpos):kpos;
#endif
    //i = __ffsll(kpos)-1;
    //i = __builtin_ctzll(kpos);
    BSFQ(i, kpos, 2);

#ifdef SEARCH_TWIN
    if ((kpos>>i) <= kmax && (kpos>>i) >= kmin && i <= ld_nstep) {
#else
    if ((kpos>>i) <= kmax && (kpos>>i) >= kmin && i < ld_nstep) {
#endif
      unsigned int new_k = (kpos>>i);
      unsigned int new_n = n+i;
      while(new_k < new_n) {
        new_k <<= 1;
        new_n--;
      }
      if(new_k == new_n) {
        cands_found++;
      //if (i < ld_nstep)
#ifdef SEARCH_TWIN
        test_factor(my_P,(kpos>>i),n+i,(kpos==k0)?-1:1);
#else
        test_factor(my_P,(kpos>>i),n+i,search_proth);
#endif
      }
    }

    // Proceed to the K for the next N.
    n += ld_nstep;
    //kpos = shiftmod_REDC(k0, my_P, k0*Ps);
    k0 = shiftmod_REDC_sm(k0, my_P, ((unsigned int)k0)*((unsigned int)Ps));
    //if(kpos != k0) printf("Expected k0=%lu, but got %lu!\n", kpos, k0);
  } while (n < l_nmax);
  if(cands_found == 0) {
    fprintf(stderr, "%sComputation Error: no candidates found for p=%"PRIu64" in %u<=n<%u.\n", bmprefix(), my_P, l_nmin, l_nmax);
#ifdef USE_BOINC
    bexit(ERR_NEG);
#endif
  }
}

INLINE void check_factors_found(const int th, const uint64_t *P, const unsigned int cthread_count) {
  unsigned int i, j;
  char factorlist;
  //fprintf(stderr, "Checking factors starting with P=%llu\n", P[0]);
  // Check the previous results.
  for(i=0; i < cthread_count; i++) {
    if((factorlist=factor_found[th][i]) != 0) {
      // Test that P, at the location(s) specified.
      for(j=0; j < 8; j++) {
        //if(factorlist & 1) test_one_p(P[i], nmin, get_n_subsection_start(j), ld_r0, ld_bbits);
        if(factorlist & 1) {
#ifndef __x86_64__
          if(nstep <= 32) 
            test_one_p_sm(P[i], get_n_subsection_start(j+1), get_n_subsection_start(j), r0arr[j], bbitsarr[j]);
          else
#endif
            test_one_p(P[i], get_n_subsection_start(j+1), get_n_subsection_start(j), r0arr[j], bbitsarr[j]);
        }
        factorlist >>= 1;
      }
    }
  }
}

/* This function is called 0 or more times in thread th, 0 <= th < num_threads.
   P is an array of APP_BUFLEN candidate primes.
*/
void app_thread_fun(int th, const uint64_t *P, uint64_t *lastP, const unsigned int cthread_count)
{
  unsigned int i;
  uint64_t new_start_time;

  // If there was a kernel running, get its results first.
  if(gpu_started[th] != (uint64_t)0) {
    //printf("Getting factors from iteration at %d\n", gpu_started[th]);
    get_factors_found(factor_found[th], cthread_count, gpu_started[th], &check_ns_delay[th]);
  }

  // Start the next kernel.
  check_ns(P, cthread_count, th);
  new_start_time = elapsed_usec();
  //printf("Checking N's for iteration starting at %d with P=%lu\n", new_start_time, P[0]);

  if(gpu_started[th] != (uint64_t)0) {
    check_factors_found(th, lastP, cthread_count);
    //printf("Checking factors for iteration starting at %d with P=%lu\n", gpu_started[th], lastP[0]);
  }

  // Copy the new P's over the old.
  for(i=0; i < cthread_count; i++) {
    lastP[i] = P[i];
  }
  gpu_started[th] = new_start_time;
}

/* This function is called 0 or more times in thread th, 0 <= th < num_threads.
   P is an array of len candidate primes, 0 <= len < APP_PRIMES_BUFLEN.
   The application should be prepared to checkpoint after returning from
   this function.
*/
void app_thread_fun1(int th, uint64_t *P, uint64_t *lastP, const unsigned int cthread_count, unsigned int len)
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

    app_thread_fun(th,P,lastP, cthread_count);
  }
  // Finish the last kernel.
  if(gpu_started[th] != (uint64_t)0) {
#ifdef TRACE
    printf("Getting factors from iteration at %d\n", (unsigned int)gpu_started[th]);
#endif
    get_factors_found(factor_found[th], cthread_count, gpu_started[th], &check_ns_delay[th]);
#ifdef TRACE
    printf("Checking factors for iteration starting at %d with P=%lu\n", (unsigned int)gpu_started[th], lastP[0]);
#endif
    check_factors_found(th, lastP, cthread_count);
    gpu_started[th] = (uint64_t)0;
  }
}

/* This function is called once in thread th, 0 <= th < num_threads, after
   the final call to app_thread_fun1(th, ...).
*/
void app_thread_fini(int th)
{
  cuda_finalize();
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

#ifndef SINGLE_THREAD
#ifdef _WIN32
  DeleteCriticalSection(&factors_mutex);
#else
  pthread_mutex_destroy(&factors_mutex);
#endif
#endif

  if (bitmap != NULL)
  {
    for (i = 0; i <= 1; i++)
      free(bitmap[i]);
    free(bitmap);
    bitmap = NULL;
  }
  //for(i=0; i < num_threads; i++) {
    //free(bitsskip[i]);
  //}
  //free(bitsskip);
  for(i=0; i < (unsigned int)num_threads; i++) {
    if(factor_found[i] != NULL) free(factor_found[i]);
  }
  free(factor_found);
}
