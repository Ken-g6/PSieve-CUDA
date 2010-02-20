/* gfn_main.c -- (C) Geoffrey Reynolds, March 2009.

   Multithreaded sieve application for algorithms of the form:

   For each odd k in 1 <= k0 <= k < k1 < 2^62
     If k*2^n+1 has no factors < q (n,q fixed)
       Do something with k


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <signal.h>
#include <math.h>
#include <getopt.h>

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif

#include "gfn_main.h"
#include "clock.h"
#include "sieve.h"
#include "util.h"

#include "gfn_app.h"


/* Global variables
 */
unsigned int num_threads = 1;
uint64_t kmin = 0, kmax = 0;
unsigned int nmin = 0, nmax = 0;

/* Local variables */
static uint64_t kstart;
static unsigned int nstart;
static unsigned int qmax = QMAX_MAX;
static unsigned int blocksize_opt = BLOCKSIZE_OPT_DEFAULT;
static unsigned int chunksize_opt = CHUNKSIZE_OPT_DEFAULT;
static unsigned int blocks_opt = BLOCKS_OPT_DEFAULT;
static sieve_t *sv;
static unsigned int report_opt = REPORT_OPT_DEFAULT;
static unsigned int checkpoint_opt = CHECKPOINT_OPT_DEFAULT;
static unsigned int priority_opt = 0;

static uint64_t cand_count = 0;
static uint64_t cand_sum = 0;
static uint64_t report_period;
static uint64_t checkpoint_period;

static uint64_t program_start_time;
static uint64_t sieve_start_time;
static uint64_t sieve_start_processor_time;
static uint64_t last_checkpoint_time;
static uint64_t last_checkpoint_progress;
static uint64_t last_report_time;
static uint64_t last_report_processor_time;
static uint64_t last_report_processor_cycles;
static uint64_t last_report_progress;

static volatile char no_more_checkpoints = 0;
static volatile char checkpointing = 0;
static volatile char stopping = 0;


/* Per-thread local variables
 */

typedef struct {
  uint64_t count;
  uint64_t sum;
  char exiting;
} thread_data_t;

static thread_data_t thread_data[MAX_THREADS];


/* Thread shared variables
 */
#ifdef _WIN32
static HANDLE checkpoint_semaphoreA;
static HANDLE checkpoint_semaphoreB;
#else
static pthread_mutex_t exiting_mutex;
static pthread_cond_t exiting_cond;
static sem_t checkpoint_semaphoreA;
static sem_t checkpoint_semaphoreB;
#endif


static void handle_signal(int signum)
{
  switch (signum)
  {
    case SIGINT:
    case SIGTERM:
#ifdef SIGHUP
    case SIGHUP:
#endif
      stopping = 1;
      break;
  }
}

static void (*old_sigint_handler)(int);
static void (*old_sigterm_handler)(int);
#ifdef SIGHUP
static void (*old_sighup_handler)(int);
#endif

static void init_signals(void)
{
  if ((old_sigint_handler = signal(SIGINT,handle_signal)) == SIG_IGN)
    signal(SIGINT,SIG_IGN);
  if ((old_sigterm_handler = signal(SIGTERM,handle_signal)) == SIG_IGN)
    signal(SIGTERM,SIG_IGN);
#ifdef SIGHUP
  if ((old_sighup_handler = signal(SIGHUP,handle_signal)) == SIG_IGN)
    signal(SIGHUP,SIG_IGN);
#endif
}

static void fini_signals(void)
{
  if (old_sigint_handler != SIG_IGN)
    signal(SIGINT,old_sigint_handler);
  if (old_sigterm_handler != SIG_IGN)
    signal(SIGTERM,old_sigterm_handler);
#ifdef SIGHUP
  if (old_sighup_handler != SIG_IGN)
    signal(SIGHUP,old_sighup_handler);
#endif
}


static void report_status(uint64_t now, uint64_t processor_time,
                          uint64_t cycles, uint64_t progress)
{
  double rate, cpus, freq, done;
  const char *unit;
  int prec;

  rate = (double)(progress-last_report_progress)/(now-last_report_time);
  unit = "M";

  if (rate < 1.0)
    rate *= 1000.0, unit = "K";

  if (rate < 1.0)
    rate *= 1000.0, unit = "";

  if (rate < 10.0)
    prec = 3;
  else if (rate < 100.0)
    prec = 2;
  else if (rate < 1000.0)
    prec = 1;
  else
    prec = 0;

  cpus = (double)(processor_time-last_report_processor_time)/(now-last_report_time);

  freq = (double)(cycles-last_report_processor_cycles)/(now-last_report_time);

  done = (double)(progress-kmin)/(kmax-kmin);

  printf("k=%"PRIu64", %.*f%s k/sec, %.2fx%.0fMHz CPU, %.1f%% done\n",
         progress,prec,rate,unit,cpus,freq,100.0*done);
}

static void write_checkpoint(uint64_t k, unsigned int n)
{
  FILE *fout;
  unsigned int i;
  uint64_t count, sum, checksum;

  if ((fout = fopen(CHECKPOINT_FILENAME,"w")) != NULL)
  {
    for (i = 0, count = cand_count, sum = cand_sum; i < num_threads; i++)
    {
      count += thread_data[i].count;
      sum += thread_data[i].sum;
    }

    checksum = kmin + k + nmin + n + count + sum;
    fprintf(fout,"kmin=%"PRIu64",k=%"PRIu64",nmin=%u,n=%u"
            ",count=%"PRIu64",sum=0x%016"PRIx64",checksum=0x%016"PRIx64"\n",
            kmin,k,nmin,n,count,sum,checksum);

    app_write_checkpoint(fout);
    fclose(fout);
  }
}

/* Try to read the checkpoint file and return the starting point. If the
   file cannot be read or the starting point (k,n) is not in kmin <= k < kmax,
   nmin <= n <= nmax, then return kmin,nmin.
*/ 
static void read_checkpoint(uint64_t *kret, unsigned int *nret)
{
  uint64_t k0, k, count, sum, checksum;
  unsigned int n0, n;
  FILE *fin;
  int valid;

  if ((fin = fopen(CHECKPOINT_FILENAME,"r")) == NULL)
  {
    *kret = kmin;
    *nret = nmin;
    return;
  }

  valid = 0;
  if (fscanf(fin,"kmin=%"SCNu64",k=%"SCNu64",nmin=%u,n=%u"
             ",count=%"SCNu64",sum=0x%"SCNx64",checksum=0x%"SCNx64"\n",
             &k0,&k,&n0,&n,&count,&sum,&checksum) == 7)
    if (k0 == kmin && k > kmin && k < kmax &&
        n0 == nmin && n >= nmin && n <= nmax)
      valid = app_read_checkpoint(fin);

  fclose(fin);

  if (valid && k0 + k + n0 + n + count + sum == checksum)
  {
    fprintf(stderr,"Resuming from checkpoint k=%"PRIu64",n=%u in %s\n",
            k, n, CHECKPOINT_FILENAME);
    cand_count = count;
    cand_sum = sum;
    *kret = k;
    *nret = n;
    return;
  }
  else
  {
    fprintf(stderr,"Ignoring invalid checkpoint in %s\n", CHECKPOINT_FILENAME);
    *kret = kmin;
    *nret = nmin;
    return;
  }
}


static const char *short_opts = "k:K:n:N:Q:B:C:c:r:t:z:h" APP_SHORT_OPTS;

static const struct option long_opts[] = {
  {"kmin",        required_argument, 0, 'k'},
  {"kmax",        required_argument, 0, 'K'},
  {"nmin",        required_argument, 0, 'n'},
  {"nmax",        required_argument, 0, 'N'},
  {"qmax",        required_argument, 0, 'Q'},
  {"blocksize",   required_argument, 0, 'B'},
  {"chunksize",   required_argument, 0, 'C'},
  {"blocks",      required_argument, 0, 256},
  {"checkpoint",  required_argument, 0, 'c'},
  {"report",      required_argument, 0, 'r'},
  {"threads",     required_argument, 0, 't'},
  {"priority",    required_argument, 0, 'z'},
  {"help",        no_argument,       0, 'h'},
  APP_LONG_OPTS
  {0,0,0,0}
};

static void help(void)
{
  printf("-k --kmin=K0       Sieve start: 1 <= K0 <= k (default K0=1)\n");
  printf("-K --kmax=K1       Sieve end: k < K1 <= 2^62 (default K1=K0+10^9)\n");
  printf("-n --nmin=N0       Sieve start: 1 <= N0 <= n\n");
  printf("-N --nmax=N1       Sieve end: n <= N1 < 2^31 (default N1=N0)\n");
  printf("-Q --qmax=Q1       Sieve only with odd primes q < Q1 <= 2^31\n");
  printf("-B --blocksize=N   Sieve in blocks of N bytes (default N=%d)\n",
         BLOCKSIZE_OPT_DEFAULT);
  printf("-C --chunksize=N   Process blocks in chunks of N bytes (default N=%d)\n", CHUNKSIZE_OPT_DEFAULT);
  printf("   --blocks=N      Sieve up to N blocks ahead (default N=%d)\n",
         BLOCKS_OPT_DEFAULT);
  printf("-c --checkpoint=N  Checkpoint every N seconds (default N=%d)\n",
         CHECKPOINT_OPT_DEFAULT);
  printf("-r --report=N      Report status every N seconds (default N=%d)\n",
         REPORT_OPT_DEFAULT);
  printf("-t --threads=N     Start N child threads (default N=1)\n");
  printf("-z --priority=N    Set process priority to nice N or {idle,low,normal}\n");
  printf("-h --help          Print this help\n");
}

static int help_opt = 0;

static int parse_option(int opt, char *arg, const char *source)
{
  int status = 0;

  switch (opt)
  {
    case 'k':
      status = parse_uint64(&kmin,arg,1,KMAX_MAX-1);
      break;

    case 'K':
      status = parse_uint64(&kmax,arg,2,KMAX_MAX);
      break;

    case 'n':
      status = parse_uint(&nmin,arg,1,NMAX_MAX);
      break;

    case 'N':
      status = parse_uint(&nmax,arg,1,NMAX_MAX);
      break;

    case 'Q':
      status = parse_uint(&qmax,arg,3,QMAX_MAX);
      break;

    case 'B':
      status = parse_uint(&blocksize_opt,arg,BLOCKSIZE_OPT_MIN,BLOCKSIZE_OPT_MAX);
      break;

    case 'C':
      status = parse_uint(&chunksize_opt,arg,CHUNKSIZE_OPT_MIN,CHUNKSIZE_OPT_MAX);
      break;

    case 256:
      status = parse_uint(&blocks_opt,arg,BLOCKS_OPT_MIN,BLOCKS_OPT_MAX);
      break;

    case 'c':
      status = parse_uint(&checkpoint_opt,arg,0,UINT32_MAX);
      break;

    case 'r':
      status = parse_uint(&report_opt,arg,0,UINT32_MAX);
      break;

    case 't':
      status = parse_uint(&num_threads,arg,1,MAX_THREADS);
      break;

    case 'z':
      if (strcmp(arg,"idle") == 0)
        priority_opt = 19;
      else if (strcmp(arg,"low") == 0)
        priority_opt = 10;
      else if (strcmp(arg,"normal") == 0)
        priority_opt = 0;
      else if (strcmp(arg,"none") == 0)
        priority_opt = -1;
      else
        status = parse_uint(&priority_opt,arg,0,19);
      priority_opt++;
      break;

    case 'h':
      help_opt = 1;
      break;

    case '?':
      status = -3;
      break;

    default:
      status = app_parse_option(opt,arg,source);
      break;
  }

  return status;
}

/* Process command-line options using getopt_long().
   Non-option arguments are treated as if they belong to option zero.
   Returns the number of options processed.
 */
static int process_args(int argc, char *argv[])
{
  int count = 0, ind = -1, opt;

  while ((opt = getopt_long(argc,argv,short_opts,long_opts,&ind)) != -1)
    switch (parse_option(opt,optarg,NULL))
    {
      case 0:
        ind = -1;
        count++;
        break;

      case -1:
        /* If ind is unchanged then this is a short option, otherwise long. */
        if (ind == -1)
          fprintf(stderr,"%s: invalid argument -%c %s\n",
                  argv[0],opt,optarg);
        else
          fprintf(stderr,"%s: invalid argument --%s %s\n",
                  argv[0],long_opts[ind].name,optarg);
        exit(EXIT_FAILURE);

      case -2:
        /* If ind is unchanged then this is a short option, otherwise long. */
        if (ind == -1)
          fprintf(stderr,"%s: out of range argument -%c %s\n",
                  argv[0],opt,optarg);
        else
          fprintf(stderr,"%s: out of range argument --%s %s\n",
                  argv[0],long_opts[ind].name,optarg);
        exit(EXIT_FAILURE);

      default:
        exit(EXIT_FAILURE);
    }

  while (optind < argc)
    switch (parse_option(0,argv[optind],NULL))
    {
      case 0:
        optind++;
        count++;
        break;

      case -1:
        fprintf(stderr,"%s: invalid non-option argument %s\n",
                  argv[0],argv[optind]);
        exit(EXIT_FAILURE);

      case -2:
        fprintf(stderr,"%s: out of range non-option argument %s\n",
                  argv[0],argv[optind]);
        exit(EXIT_FAILURE);

      default:
        exit(EXIT_FAILURE);
    }

  if (help_opt)
  {
    help();
    app_help();
    exit(EXIT_SUCCESS);
  }

  return count;
}

/* Read and parse options from configuration file fn.
   Returns the number of options read, or zero if the file cannot be opened.
*/
static int read_config_file(const char *fn)
{
  const char comment_character = '#';
  const char *delimiters = "= \n\r\t\v";
  char line[128];
  char *str, *arg;
  int ind, count;
  FILE *file;

  assert(fn != NULL);

  if ((file = fopen(fn,"r")) == NULL)
    return 0;

  for (count = 0; fgets(line,sizeof(line),file) != NULL; )
  {
    str = strtok(line,delimiters);

    if (str == NULL || str[0] == comment_character)
      continue;

    arg = strtok(NULL,delimiters);

    for (ind = 0; long_opts[ind].name; ind++)
      if (strcmp(str,long_opts[ind].name) == 0)
        break;

    if (long_opts[ind].name == NULL)
    {
      fprintf(stderr,"%s: unrecognised option `%s'\n",fn,str);
      exit(EXIT_FAILURE);
    }

    if (long_opts[ind].has_arg == no_argument && arg != NULL)
    {
      fprintf(stderr,"%s: option `%s' doesn't allow an argument\n",fn,str);
      exit(EXIT_FAILURE);
    }
    else if (long_opts[ind].has_arg == required_argument && arg == NULL)
    {
      fprintf(stderr,"%s: option `%s' requires an argument\n",fn,str);
      exit(EXIT_FAILURE);
    }

    if (long_opts[ind].flag != NULL)
    {
      *long_opts[ind].flag = long_opts[ind].val;
    }
    else switch (parse_option(long_opts[ind].val,arg,fn))
    {
      case 0:
        break;

      case -1:
        fprintf(stderr,"%s: invalid argument %s %s\n",fn,str,arg);
        exit(EXIT_FAILURE);

      case -2:
        fprintf(stderr,"%s: out of range argument %s %s\n",fn,str,arg);
        exit(EXIT_FAILURE);

      default:
        exit(EXIT_FAILURE);
    }

    count++;
  }

  fclose(file);

  return count;
}


#ifndef _WIN32
/* Child thread cleanup handler signals parent before child thread exits.
   This is needed because the pthreads API lacks the equivalent of a select
   function to wait (with timeout) for one of a number of threads to exit.
*/
static void thread_cleanup(void *arg)
{
  int th = (int)arg;

#ifndef NDEBUG
  fprintf(stderr,"thread_cleanup: %d\n",th);
#endif

  pthread_mutex_lock(&exiting_mutex);
  thread_data[th].exiting = 1;
  pthread_cond_signal(&exiting_cond);
  pthread_mutex_unlock(&exiting_mutex);
}
#endif

static void *thread_fun(void *arg)
{
  int th = (int)arg;
  uint64_t k0 = 0, k, count = 0, sum = 0;
  uint64_t K[APP_BUFLEN] __attribute__ ((aligned(16)));
  unsigned int klen, len;
  unsigned long *buf;

  assert(((unsigned int)&K & 15) == 0); /* check stack alignment */

#ifdef _WIN32
  if (priority_opt)
  {
    if (priority_opt > 14)
      SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_IDLE);
    else if (priority_opt > 7)
      SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_BELOW_NORMAL);
    else
      SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_NORMAL);
  }
#endif

#ifndef _WIN32
  pthread_cleanup_push(thread_cleanup,arg);
#endif

  fprintf(stderr,"Thread %d starting\n",th);

  klen = 0;
  len = sv->chunk_size;
  buf = NULL;

  app_thread_init(th);

  while (!stopping)
  {
    unsigned int i;

#if TRACE
    if ((k0 = get_chunk(th,sv,&buf)) >= kmax)
#else
    if ((k0 = get_chunk(sv,&buf)) >= kmax)
#endif
      break;

    for (i = 0; i < len; i++)
    {
      unsigned int j;
      unsigned long u;

      for (j = 0, u = buf[i]; u != 0; u >>= 1, j++)
      {
#ifdef __GNUC__
        int z = __builtin_ctzl(u);
        u >>= z;
        j += z;
#else
        if (!(u & 1UL))
          continue;
#endif
        count++;
        k = k0+(uint64_t)(ULONG_BIT*i+j)*2;
        sum += k;
        K[klen++] = k;
        if (klen == APP_BUFLEN)
          app_thread_fun(th,K), klen = 0;
      }
    }

#if TRACE
    free_chunk(th,sv,k0);
#else
    free_chunk(sv,k0);
#endif

    if (checkpointing)
    {
      app_thread_fun1(th,K,klen), klen = 0;
      thread_data[th].count = count;
      thread_data[th].sum = sum;
#ifdef TRACE
      printf("Thread %d: Synchronising for checkpoint\n",th);
#endif
#ifdef _WIN32
      ReleaseSemaphore(checkpoint_semaphoreA,1,NULL);
      WaitForSingleObject(checkpoint_semaphoreB,INFINITE);
#else
      sem_post(&checkpoint_semaphoreA);
      sem_wait(&checkpoint_semaphoreB);
#endif
#ifdef TRACE
      printf("Thread %d: Continuing after checkpoint\n",th);
#endif
    }
  }

  app_thread_fun1(th,K,klen);
  thread_data[th].count = count;
  thread_data[th].sum = sum;
  app_thread_fini(th);

  if (k0 >= kmax)
    fprintf(stderr,"Thread %d completed\n",th);
  else
    fprintf(stderr,"Thread %d interrupted\n",th);

  /* Just in case a checkpoint is signalled before other threads exit */
  no_more_checkpoints = 1;
#ifdef _WIN32
  ReleaseSemaphore(checkpoint_semaphoreA,1,NULL);
#else
  sem_post(&checkpoint_semaphoreA);
#endif

#ifndef _WIN32
  pthread_cleanup_pop(1);
#endif

  return 0;
}

#ifdef _WIN32
static unsigned int __stdcall thread_fun_wrapper(void *arg)
{
#ifdef __i386__
  /* _beginthreadex doesn't align the stack */
  asm ("push %%ebp\n\t"
       "mov %%esp, %%ebp\n\t"
       "and $-16, %%esp\n\t"
       "sub $12, %%esp\n\t"
       "push %0\n\t"
       "call _thread_fun\n\t"
       "mov %%ebp, %%esp\n\t"
       "pop %%ebp"
       : "+a" (arg) : "i" (thread_fun) : "%edx", "%ecx", "cc");

  return (unsigned int)arg;
#else
  return (unsigned int)thread_fun(arg);
#endif
}
#endif


int main(int argc, char *argv[])
{
#ifdef _WIN32
  HANDLE tid[MAX_THREADS];
  DWORD thread_ret;
#else
  pthread_t tid[MAX_THREADS];
  void *thread_ret;
  int joined;
#endif
  uint64_t kstop;
  unsigned int ncurr;
  int th, process_ret = EXIT_SUCCESS;

  program_start_time = elapsed_usec();
  app_banner();

  read_config_file(CONFIG_FILENAME);
  process_args(argc,argv);

  if (kmin < KMIN_MIN)
    kmin = KMIN_MIN;
  if (kmax > KMAX_MAX)
    kmax = KMAX_MAX;
  if (kmin >= kmax)
  {
    if (kmax == 0 && kmin < KMAX_MAX-1000000000)
    {
      fprintf(stderr,"kmax not specified, using default kmax = kmin + 1e9\n");
      kmax = kmin + 1000000000; /* Default range */
    }
    else
    {
      fprintf(stderr, "Option out of range: kmax must be greater than kmin\n");
      exit(EXIT_FAILURE);
    }
  }

  if (nmin < NMIN_MIN)
    nmin = NMIN_MIN;
  if (nmax > NMAX_MAX)
    nmax = NMAX_MAX;
  if (nmin > nmax)
  {
    if (nmax == 0 && nmin <= NMAX_MAX)
    {
      fprintf(stderr,"nmax not specified, using default nmax = nmin\n");
      nmax = nmin; /* Default range */
    }
    else
    {
      fprintf(stderr, "Option out of range: nmax must be no less than nmin\n");
      exit(EXIT_FAILURE);
    }
  }

  checkpoint_period = (uint64_t)checkpoint_opt * 1000000; /* usec */
  report_period = (uint64_t)report_opt * 1000000; /* usec */

  if (blocksize_opt/chunksize_opt < num_threads)
  {
    chunksize_opt = blocksize_opt/num_threads;
    if (chunksize_opt < CHUNKSIZE_OPT_MIN)
    {
      chunksize_opt = CHUNKSIZE_OPT_MIN;
      blocksize_opt = CHUNKSIZE_OPT_MIN*num_threads;
    }
  }

  if (priority_opt)
  {
#ifdef _WIN32
    if (priority_opt > 14)
      SetPriorityClass(GetCurrentProcess(),IDLE_PRIORITY_CLASS);
    else if (priority_opt > 7)
      SetPriorityClass(GetCurrentProcess(),BELOW_NORMAL_PRIORITY_CLASS);
    else
      SetPriorityClass(GetCurrentProcess(),NORMAL_PRIORITY_CLASS);
#else
    setpriority(PRIO_PROCESS,0,priority_opt-1);
#endif
  }

  if (nmax >= 62 || ((kmax<<nmax)>>nmax) != kmax || (kmax<<nmax) >= UINT64_C(1)<<62)
  {
    fprintf(stderr,"kmax,nmax out of range (k*2^n+1 must be less than 2^62)\n");
    exit(EXIT_FAILURE);
  }

  uint64_t pmax = (kmax<<nmax)+1;
  if ((uint64_t)qmax*qmax > pmax)
    qmax = sqrt(pmax);
  init_sieve_primes(qmax);

  app_init();

  fprintf(stderr,"Sieve started: %"PRIu64" <= k < %"PRIu64", %u <= n <= %u\n",
          kmin,kmax,nmin,nmax);

  read_checkpoint(&kstart,&nstart);
  kstart |= 1; /* Must be odd! */
  ncurr = nstart;

  if (nstart != nmin || nmax != nmin)
  {
    fprintf(stderr, "Multiple n not implemented\n");
    exit(EXIT_FAILURE);
  }

  /* For each n */

  sv = create_gfn_sieve(kstart,kmax,nstart,qmax,chunksize_opt,blocksize_opt,blocks_opt);

  init_signals();

#ifdef _WIN32
  checkpoint_semaphoreA = CreateSemaphore(NULL,0,2147483647,NULL);
  checkpoint_semaphoreB = CreateSemaphore(NULL,0,2147483647,NULL);
#else
  pthread_mutex_init(&exiting_mutex,NULL);
  pthread_cond_init(&exiting_cond,NULL);
  sem_init(&checkpoint_semaphoreA,0,0);
  sem_init(&checkpoint_semaphoreB,0,0);
#endif

  sieve_start_time = elapsed_usec();
  sieve_start_processor_time = processor_usec();

  /* Start child threads */
#ifdef _WIN32
  for (th = 0; th < num_threads; th++)
    if ((tid[th] = (HANDLE)
         _beginthreadex(NULL,0,thread_fun_wrapper,(void *)th,0,NULL)) == 0)
    {
      perror("_beginthreadex");
      exit(EXIT_FAILURE);
    }
#else
  pthread_mutex_lock(&exiting_mutex);
  for (th = 0; th < num_threads; th++)
    if (pthread_create(&tid[th],NULL,thread_fun,(void *)th) != 0)
    {
      perror("pthread_create");
      exit(EXIT_FAILURE);
    }
#endif

  last_checkpoint_time = sieve_start_time;
  last_checkpoint_progress = kstart;
  last_report_time = sieve_start_time;
  last_report_processor_time = sieve_start_processor_time;
  last_report_processor_cycles = processor_cycles();
  last_report_progress = kstart;

  while (!stopping)
  {
    uint64_t current_time, processor_time, cycles, progress;
    uint64_t next_checkpoint_time, next_report_time, timeout;

    current_time = elapsed_usec();
    processor_time = processor_usec();
    cycles = processor_cycles();

    next_report_time = last_report_time + report_period;
    if (current_time >= next_report_time)
    {
      progress = next_chunk(sv);
      report_status(current_time,processor_time,cycles,progress);
      last_report_time = current_time;
      next_report_time = current_time + report_period;
      last_report_processor_time = processor_time;
      last_report_processor_cycles = cycles;
      last_report_progress = progress;
    }

    if (checkpoint_opt)
    {
      next_checkpoint_time = last_checkpoint_time + checkpoint_period;
      if (current_time >= next_checkpoint_time && !no_more_checkpoints)
      {
#if TRACE
        printf("Main: checkpointing\n");
#endif
        checkpointing = 1;
#ifdef _WIN32
        for (th = 0; th < num_threads; th++)
          WaitForSingleObject(checkpoint_semaphoreA,INFINITE);
#else
        for (th = 0; th < num_threads; th++)
          sem_wait(&checkpoint_semaphoreA);
#endif
        progress = next_chunk(sv);
        write_checkpoint(progress,ncurr);
        last_checkpoint_progress = progress;
        last_checkpoint_time = current_time;
        next_checkpoint_time = current_time + checkpoint_period;
        checkpointing = 0;
#if TRACE
        printf("Main: finished checkpointing\n");
#endif
#ifdef _WIN32
        ReleaseSemaphore(checkpoint_semaphoreB,num_threads,NULL);
#else
        for (th = 0; th < num_threads; th++)
          sem_post(&checkpoint_semaphoreB);
#endif
      }
    }

    if (checkpoint_opt && next_checkpoint_time < next_report_time)
      timeout = next_checkpoint_time;
    else
      timeout = next_report_time;

    /* Wait until timeout, or some child thread exits. */
#ifdef _WIN32
    DWORD timeout_interval = (DWORD)(timeout - current_time+999)/1000;
    /* Wait for any thread */
    if (WaitForMultipleObjects(num_threads,tid,0,timeout_interval)
        == WAIT_OBJECT_0)
    {
      /* Find which thread exited */
      for (th = 0; th < num_threads; th++)
        if (WaitForSingleObject(tid[th],0) == WAIT_OBJECT_0)
        {
          /* If this thread failed, stop the others too. */
          if (GetExitCodeThread(tid[th],&thread_ret) && thread_ret != 0)
            stopping = 1;
          break;
        }
      break;
    }
#else
    struct timespec wait_timespec;
    wait_timespec.tv_sec = timeout/1000000;
    wait_timespec.tv_nsec = (timeout%1000000)*1000;
    if (pthread_cond_timedwait(&exiting_cond,&exiting_mutex,&wait_timespec)==0)
      break;
#endif
  }

  /* Restore signal handlers in case some thread fails to join below. */
  fini_signals();

  fprintf(stderr,"Waiting for threads to exit\n");

#ifdef _WIN32
  /* Wait for all threads, then examine return values and close. */
  WaitForMultipleObjects(num_threads,tid,1,INFINITE);
  for (th = 0; th < num_threads; th++)
  {
    if (GetExitCodeThread(tid[th],&thread_ret) && thread_ret != 0)
    {
      fprintf(stderr,"Thread %d failed: %lX\n",th,thread_ret);
      process_ret = EXIT_FAILURE;
    }
    CloseHandle(tid[th]);
  }
#else
  pthread_mutex_unlock(&exiting_mutex);
  /* Find an exiting thread, if there is one (there might not be if the loop
     above was exited because of a signal that set the stopping flag) */
  joined = num_threads;
  for (th = 0; th < num_threads; th++)
  {
    if (thread_data[th].exiting)
    {
      pthread_join(tid[th],&thread_ret);
      if (thread_ret != 0)
      {
        /* This thread exited with an error, so stop the others too */
        fprintf(stderr,"Thread %d failed: %p\n",th,thread_ret);
        process_ret = EXIT_FAILURE;
        stopping = 1;
      }
      joined = th; /* Note which thread was joined. */
      break;
    }
  }

  /* Join any remaining threads. If joined < num_threads then skip the
     thread that was already joined above. */
  for (th = 0; th < num_threads; th++)
  {
    if (th != joined)
    {
      pthread_join(tid[th],&thread_ret);
      if (thread_ret != 0)
      {
        fprintf(stderr,"Thread %d failed: %p\n",th,thread_ret);
        process_ret = EXIT_FAILURE;
      }
    }
  }
#endif

#ifdef _WIN32
  CloseHandle(checkpoint_semaphoreA);
  CloseHandle(checkpoint_semaphoreB);
#else
  pthread_cond_destroy(&exiting_cond);
  pthread_mutex_destroy(&exiting_mutex);
  sem_destroy(&checkpoint_semaphoreA);
  sem_destroy(&checkpoint_semaphoreB);
#endif

  if (process_ret == EXIT_SUCCESS)
  {
    kstop = next_chunk(sv);
    write_checkpoint(kstop,ncurr);
  }
  else
  {
    kstop = last_checkpoint_progress;
  }

  if (kstop >= kmax)
  {
    fprintf(stderr,"Sieve complete: %"PRIu64" <= k < %"PRIu64"\n",kmin,kmax);
    remove(CHECKPOINT_FILENAME);
  }
  else
  {
   fprintf(stderr,"Sieve incomplete: %"PRIu64" <= k < %"PRIu64"\n",kmin,kstop);
  }

  app_fini();

  destroy_sieve(sv);
  free_sieve_primes();

  if (process_ret == EXIT_SUCCESS)
  {
    /* Print final candidate count/sum */
    for (th = 0; th < num_threads; th++)
    {
      cand_count += thread_data[th].count;
      cand_sum += thread_data[th].sum;
    }
    fprintf(stderr,"count=%"PRIu64",sum=0x%016"PRIx64"\n",cand_count,cand_sum);
  }

  /* Print statistics for this run */
  uint64_t stop_time = elapsed_usec();
  uint64_t stop_processor_time = processor_usec();
  fprintf(stderr,"Elapsed time: %.2f sec. (%.2f init + %.2f sieve)"
          " at %.0f k/sec.\n",
          (stop_time-program_start_time)/1000000.0,
          (sieve_start_time-program_start_time)/1000000.0,
          (stop_time-sieve_start_time)/1000000.0,
          (double)(kstop-kstart)/(stop_time-sieve_start_time)*1000000);
  fprintf(stderr,"Processor time: %.2f sec. (%.2f init + %.2f sieve)"
          " at %.0f k/sec.\n",
          (stop_processor_time)/1000000.0,
          (sieve_start_processor_time)/1000000.0,
          (stop_processor_time-sieve_start_processor_time)/1000000.0,
          (double)(kstop-kstart)/(stop_processor_time-sieve_start_processor_time)*1000000);
  fprintf(stderr,"Average processor utilization: %.2f (init), %.2f (sieve)\n",
          (double)(sieve_start_processor_time)
          /(sieve_start_time-program_start_time),
          (double)(stop_processor_time-sieve_start_processor_time)
          /(stop_time-sieve_start_time));

  return process_ret;
}
