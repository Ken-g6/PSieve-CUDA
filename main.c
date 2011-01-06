/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* main.c -- (C) Geoffrey Reynolds, March 2009.
 * and some  (C) Ken Brazier, October 2009-July 2010.

   Multithreaded sieve application for algorithms of the form:

   For each prime p in 3 <= p0 <= p < p1 < 2^62
     Do something with p


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stdint.h"
#include "inttypes.h"
#include <signal.h>
#include <math.h>
#include "getopt.h"
#include <time.h>
#include "main.h"

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#ifndef SINGLE_THREAD
#include <pthread.h>
#include <semaphore.h>
#endif
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif

#include "clock.h"
#include "sieve.h"
#include "putil.h"

#include "app.h"
#include "appcu.h"
#ifdef USE_BOINC
#ifdef _WIN32                //  Stuff we only need on Windows: 
#include "boinc_win.h"
#endif

/* BOINC API */
#if defined(__APPLE__) || defined(_WIN32)
#include "boinc_api.h"
//#include "diagnostics.h"
#include "filesys.h"
#else
#include "boinc_api.h"
//#include "diagnostics.h"     // boinc_init_diagnostics()
// Note: filesys.h has some improperly #ifdef-ed C++ #includes, when building on Win32 with MinGW/MSys.
#include "filesys.h"         // boinc_fopen(), etc...
#endif
#include "do_boinc_init.h"
#endif

/* Global variables
 */
int num_threads = 1;
uint64_t pmin = 0, pmax = 0;
unsigned int quiet_opt = 0;
unsigned int priority_opt = 0;
uint64_t pstart;

/* Local variables */
static char *pmin_str = NULL;
static char *checkpoint_filename = NULL;
static /*const*/ char *empty_string = "";
static unsigned int qmax = QMAX_MAX;
static unsigned int blocksize_opt = BLOCKSIZE_OPT_DEFAULT;
static unsigned int chunksize_opt = CHUNKSIZE_OPT_DEFAULT;
static unsigned int blocks_opt = BLOCKS_OPT_DEFAULT;
static sieve_t *sv;
static unsigned int report_opt = REPORT_OPT_DEFAULT;
static unsigned int checkpoint_opt = CHECKPOINT_OPT_DEFAULT;

static uint64_t cand_count = 0;
static uint64_t cand_sum = 0;
static uint64_t report_period;
static uint64_t checkpoint_period;

static uint64_t program_start_time;
static time_t sieve_start_date;
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
#ifndef SINGLE_THREAD
#ifdef _WIN32
static HANDLE checkpoint_semaphoreA;
static HANDLE checkpoint_semaphoreB;
#else
static pthread_mutex_t exiting_mutex;
static pthread_cond_t exiting_cond;
static sem_t checkpoint_semaphoreA;
static sem_t checkpoint_semaphoreB;
#endif
#endif

static int got_sigint = 0, got_sigterm=0, got_sighup=0;
static void handle_signal(int signum)
{
  switch (signum)
  {
    case SIGINT: got_sigint = 1; stopping = 1; break;
    case SIGTERM: got_sigterm = 1; stopping = 1; break;
#ifdef SIGHUP
    case SIGHUP: got_sighup = 1; stopping = 1; break;
#endif
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

void boincordie(int res, char *message) {
  if(res) {
    fprintf(stderr, "%s", message);
    bexit(EXIT_FAILURE);
  }
}

#ifdef USE_BOINC
FILE* bfopen(const char *filename, const char *mode) {
  char resolved_name[512];
  if(boinc_resolve_filename(filename, resolved_name, 512)) return NULL;
  return boinc_fopen(resolved_name, mode);
}
void bmsg(const char *msg) {
  //fprintf(stderr, "%s %s", boinc_msg_prefix(), msg);
  fprintf(stderr, "%s", msg);
}
char* bmprefix() {
  //return boinc_msg_prefix();
  return "";
}
/*
 *  Dummy graphics API entry points.  This app does not do graphics, 
 *  but it still must provide these empty callbacks, for BOINC 5.  
 */
#if BOINC_MAJOR_VERSION < 6
void app_graphics_init() {}
void app_graphics_resize(int width, int height){}
void app_graphics_render(int xs, int ys, double time_of_day) {}
void app_graphics_reread_prefs() {}
void boinc_app_mouse_move(int x, int y, int left, int middle, int right ){}
void boinc_app_mouse_button(int x, int y, int which, int is_down){}
void boinc_app_key_press(int x, int y){}
void boinc_app_key_release(int x, int y){}
#endif

#else
FILE* bfopen(const char *filename, const char *mode) {
  return fopen(filename, mode);
}
void bmsg(const char *msg) {
  fprintf(stderr, "%s", msg);
}
char* bmprefix() {
  return "";
}
int boinc_init() {
  fprintf(stderr, "%sStarting non-BOINC version of PPSieve", bmprefix());
  return 0;
}
void boinc_finish(int status) {
  exit(status);
}
int boinc_time_to_checkpoint() {
  return 1;
}
void boinc_checkpoint_completed() {
}
#endif
void bexit(int status) {
  boinc_finish(status);
}

#define STRFTIME_FORMAT "ETA %d %b %H:%M"
static void report_status(uint64_t now, uint64_t processor_time,
                          uint64_t cycles, uint64_t progress)
{
  double done;
#ifndef USE_BOINC
  double rate, cpus, freq, done_eta;
  const char *unit;
  int prec;
  char buf[32];
#endif

  done = (double)(progress-pmin)/(pmax-pmin);
#ifdef USE_BOINC
  boinc_fraction_done(done);
//#endif
#else
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

  done_eta = (double)(progress-pstart)/(pmax-pstart);

  // Calculate ETA.
  buf[0] = '\0';
  if (done_eta > 0.0) /* Avoid division by zero */
  {
    time_t finish_date = sieve_start_date+(time(NULL)-sieve_start_date)/done_eta;
    struct tm *finish_tm = localtime(&finish_date);
    if (!finish_tm || !strftime(buf,sizeof(buf),STRFTIME_FORMAT,finish_tm))
      buf[0] = '\0';
  }
  printf("p=%"PRIu64", %.*f%s p/sec, %.2f CPU cores, %.1f%% done. %s  ",
         progress,prec,rate,unit,cpus,100.0*done, buf);
  if(quiet_opt) { 
    putchar('\r');
    fflush(stdout);
  } else putchar('\n');
#endif
}

static void write_checkpoint(uint64_t p)
{
  FILE *fout;
  unsigned int i;
  uint64_t count, sum, checksum;

  if ((fout = fopen(checkpoint_filename,"w")) != NULL)
  {
    for (i = 0, count = cand_count, sum = cand_sum; i < (unsigned int)num_threads; i++)
    {
      count += thread_data[i].count;
      sum += thread_data[i].sum;
    }

    checksum = pmin + p + count + sum;
    fprintf(fout,"pmin=%"PRIu64",p=%"PRIu64
            ",count=%"PRIu64",sum=0x%016"PRIx64",checksum=0x%016"PRIx64"\n",
            pmin,p,count,sum,checksum);

    app_write_checkpoint(fout);
    fclose(fout);
  }
}

/* Try to read the checkpoint file and return the starting point. If the
   file cannot be read or the starting point p is not in pmin <= p < pmax,
   then return pmin.
*/ 
static uint64_t read_checkpoint(void)
{
  uint64_t p0, p, count, sum, checksum;
  char *cpf = checkpoint_filename;
  FILE *fin;
  int valid;

  if ((fin = fopen(cpf,"r")) == NULL) {
#ifdef OLD_CHECKPOINT_FILENAME
    if ((fin = fopen(OLD_CHECKPOINT_FILENAME,"r")) == NULL)
#endif
      return pmin;
#ifdef OLD_CHECKPOINT_FILENAME
    else
      cpf = OLD_CHECKPOINT_FILENAME;
#endif
  }

  valid = 0;
  if (fscanf(fin,"pmin=%"SCNu64",p=%"SCNu64
             ",count=%"SCNu64",sum=0x%"SCNx64",checksum=0x%"SCNx64"\n",
             &p0,&p,&count,&sum,&checksum) == 5)
    if (p0 == pmin && p > pmin && p < pmax)
      valid = app_read_checkpoint(fin);

  fclose(fin);

  if (valid && p0 + p + count + sum == checksum)
  {
    fprintf(stderr,"%sResuming from checkpoint p=%"PRIu64" in %s\n",
            bmprefix(), p, cpf);
    cand_count = count;
    cand_sum = sum;
    return p;
  }
  else
  {
    fprintf(stderr,"%sIgnoring invalid checkpoint in %s\n", bmprefix(), cpf);
    return pmin;
  }
}


#ifdef SINGLE_THREAD
static const char *short_opts = "p:P:Q:B:C:c:s:r:z:qh" APP_SHORT_OPTS;
#else
static const char *short_opts = "p:P:Q:B:C:c:s:r:t:z:qh" APP_SHORT_OPTS;
#endif
static const struct option long_opts[] = {
  {"pmin",        required_argument, 0, 'p'},
  {"pmax",        required_argument, 0, 'P'},
  {"qmax",        required_argument, 0, 'Q'},
  {"blocksize",   required_argument, 0, 'B'},
  {"chunksize",   required_argument, 0, 'C'},
  {"blocks",      required_argument, 0, 256},
  {"checkpoint",  required_argument, 0, 'c'},
  {"savepoint",   required_argument, 0, 's'},
  {"report",      required_argument, 0, 'r'},
#ifndef SINGLE_THREAD
  {"nthreads",    required_argument, 0, 't'},
#endif
  {"priority",    required_argument, 0, 'z'},
  //{"device",      required_argument, 0, 'd'},
  {"help",        no_argument,       0, 'h'},
  {"quiet",       no_argument,       0, 'q'},
  APP_LONG_OPTS
  {0,0,0,0}
};

static void help(void)
{
  printf("-p --pmin=P0       Sieve start: 3 <= P0 <= p (default P0=3)\n");
  printf("-P --pmax=P1       Sieve end: p < P1 <= 2^62 (default P1=P0+10^9)\n");
  printf("-Q --qmax=Q1       Sieve only with odd primes q < Q1 <= 2^31\n");
  printf("-B --blocksize=N   Sieve in blocks of N bytes (default N=%d)\n",
         BLOCKSIZE_OPT_DEFAULT);
  printf("-C --chunksize=N   Process blocks in chunks of N bytes (default N=%d)\n", CHUNKSIZE_OPT_DEFAULT);
  printf("   --blocks=N      Sieve up to N blocks ahead (default N=%d)\n",
         BLOCKS_OPT_DEFAULT);
  printf("-c --checkpoint=N  Checkpoint every N seconds (default N=%d)\n",
         CHECKPOINT_OPT_DEFAULT);
  printf("-s --savepoint=N   Checkpoint every N minutes\n");
  printf("-q --quiet         Don't print factors to screen\n");
  printf("-r --report=N      Report status every N seconds (default N=%d)\n",
         REPORT_OPT_DEFAULT);
#ifndef SINGLE_THREAD
  printf("-t --nthreads=N    Start N threads on N GPUs (default N=1)\n");
#endif
  printf("-z --priority=N    Set process priority to nice N or {idle,low,normal}\n");
  printf("-h --help          Print this help\n");
}

static int help_opt = 0;

static int parse_option(int opt, char *arg, const char *source)
{
  int status = 0;

  switch (opt)
  {
    case 'p':
      astrcpy(&pmin_str, arg);
      status = parse_uint64(&pmin,arg,3,PMAX_MAX-1);
      break;

    case 'P':
      status = parse_uint64(&pmax,arg,4,PMAX_MAX);
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

    case 's':
      status = parse_uint(&checkpoint_opt,arg,0,UINT32_MAX);
      checkpoint_opt *= 60;
      break;

    case 'r':
      status = parse_uint(&report_opt,arg,0,UINT32_MAX);
      break;

#ifndef SINGLE_THREAD
    case 't':
      status = parse_uint(&num_threads,arg,1,MAX_THREADS);
      break;
#endif

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

    case 'q':
      quiet_opt = 1;
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
          fprintf(stderr,"%s%s: invalid argument -%c %s\n", bmprefix(),
                  argv[0],opt,optarg);
        else
          fprintf(stderr,"%s%s: invalid argument --%s %s\n", bmprefix(),
                  argv[0],long_opts[ind].name,optarg);
        boinc_finish(EXIT_FAILURE);

      case -2:
        /* If ind is unchanged then this is a short option, otherwise long. */
        if (ind == -1)
          fprintf(stderr,"%s%s: out of range argument -%c %s\n", bmprefix(),
                  argv[0],opt,optarg);
        else
          fprintf(stderr,"%s%s: out of range argument --%s %s\n", bmprefix(),
                  argv[0],long_opts[ind].name,optarg);
        boinc_finish(EXIT_FAILURE);

      default:
        boinc_finish(EXIT_FAILURE);
    }

  while (optind < argc)
    switch (parse_option(0,argv[optind],NULL))
    {
      case 0:
        optind++;
        count++;
        break;

      case -1:
        fprintf(stderr,"%s%s: invalid non-option argument %s\n", bmprefix(),
                  argv[0],argv[optind]);
        boinc_finish(EXIT_FAILURE);

      case -2:
        fprintf(stderr,"%s%s: out of range non-option argument %s\n", bmprefix(),
                  argv[0],argv[optind]);
        boinc_finish(EXIT_FAILURE);

      default:
        boinc_finish(EXIT_FAILURE);
    }

  if (help_opt)
  {
    help();
    app_help();
    boinc_finish(EXIT_SUCCESS);
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

  if ((file = bfopen(fn,"r")) == NULL)
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
      fprintf(stderr,"%s%s: unrecognised option `%s'\n", bmprefix(),fn,str);
      boinc_finish(EXIT_FAILURE);
    }

    if (long_opts[ind].has_arg == no_argument && arg != NULL)
    {
      fprintf(stderr,"%s%s: option `%s' doesn't allow an argument\n", bmprefix(),fn,str);
      boinc_finish(EXIT_FAILURE);
    }
    else if (long_opts[ind].has_arg == required_argument && arg == NULL)
    {
      fprintf(stderr,"%s%s: option `%s' requires an argument\n", bmprefix(),fn,str);
      boinc_finish(EXIT_FAILURE);
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
        fprintf(stderr,"%s%s: invalid argument %s %s\n", bmprefix(),fn,str,arg);
        boinc_finish(EXIT_FAILURE);

      case -2:
        fprintf(stderr,"%s%s: out of range argument %s %s\n", bmprefix(),fn,str,arg);
        boinc_finish(EXIT_FAILURE);

      default:
        fprintf(stderr,"%s%s: weird argument %s %s\n", bmprefix(),fn,str,arg);
        boinc_finish(EXIT_FAILURE);
    }

    count++;
  }

  fclose(file);

  return count;
}


#ifndef SINGLE_THREAD
#ifndef _WIN32
/* Child thread cleanup handler signals parent before child thread exits.
   This is needed because the pthreads API lacks the equivalent of a select
   function to wait (with timeout) for one of a number of threads to exit.
*/
static void thread_cleanup(void *arg)
{
  int th = (int)((long)arg);

#ifndef NDEBUG
  fprintf(stderr,"%sthread_cleanup: %d\n",bmprefix(),th);
#endif

  pthread_mutex_lock(&exiting_mutex);
  thread_data[th].exiting = 1;
  pthread_cond_signal(&exiting_cond);
  pthread_mutex_unlock(&exiting_mutex);
}
#endif
#endif

static void *thread_fun(void *arg)
{
  int th = (int)((long)arg);
  uint64_t p0 = 0, p, count = 0, sum = 0;
  uint64_t *P, *K, *P0, *K0;
  unsigned char *P1;
  unsigned int plen, len;
  unsigned long *buf;
  unsigned int cthread_count;
#ifdef SINGLE_THREAD
  uint64_t progress;
  // 1 if a CUDA block just finished, making this a good point to report/save status.
  int good_break_point = 0;
#endif

  //fprintf(stderr,"Thread %d starting\n",th);

#ifndef SINGLE_THREAD
#ifdef _WIN32
#ifdef USE_BOINC
  //SetPriorityClass(GetCurrentProcess(), NORMAL_PRIORITY_CLASS);
  //SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_NORMAL);
#else
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
#endif

#ifndef _WIN32
  pthread_cleanup_push(thread_cleanup,arg);
#endif
#endif

  fprintf(stderr,"%sThread %d starting\n", bmprefix(),th);

  plen = 0;
  len = sv->chunk_size;
  buf = NULL;

  cthread_count = app_thread_init(th);
  if(cthread_count > 0) {
    // Create the P array, to match the number of CUDA threads.
    P0 = (uint64_t*)xmalloc((cthread_count+2)*sizeof(uint64_t));
    P1 = (unsigned char*)P0;
    while((((unsigned long)P1) & 15) != 0) P1++; /* set stack alignment */
    P = (uint64_t*) P1;
    //fprintf(stderr,"Malloc 1 done.\n",th);

    // Create the K array, to match the P array.
    K0 = (uint64_t*)xmalloc((cthread_count+2)*sizeof(uint64_t));
    P1 = (unsigned char*)K0;
    while((((unsigned long)P1) & 15) != 0) P1++; /* set stack alignment */
    K = (uint64_t*) P1;
    //fprintf(stderr,"Malloc 2 done.\n",th);
    while (!stopping)
    {
      unsigned int i;

#if TRACE
      if ((p0 = get_chunk(th,sv,&buf)) >= pmax)
#else
        if ((p0 = get_chunk(sv,&buf)) >= pmax)
#endif
          break;

#ifdef SINGLE_THREAD
      good_break_point = 0;
#endif
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
          p = p0+(uint64_t)(ULONG_BIT*i+j)*2;
          sum += p;
          P[plen++] = p;
          if (plen == cthread_count) {
            app_thread_fun(th,P,K,cthread_count);
            plen = 0;
#ifdef SINGLE_THREAD
            good_break_point = 1;
#endif
          }
        }
      }

#if TRACE
      free_chunk(th,sv,p0);
#else
      free_chunk(sv,p0);
#endif

#ifdef SINGLE_THREAD
      if(good_break_point) {
        progress = next_chunk(sv);
        report_status(0,0,0,progress);
#ifdef USE_BOINC
        if(boinc_time_to_checkpoint())
#else
      if (checkpointing)
#endif
#else
      if (checkpointing)
#endif
      {
        app_thread_fun1(th,P,K,cthread_count,plen), plen = 0;
        thread_data[th].count = count;
        thread_data[th].sum = sum;
#ifdef TRACE
        printf("Thread %d: Synchronising for checkpoint\n",th);
#endif
#ifdef USE_BOINC
        /* Without threads, write a checkpoint whenever done with a chunk and BOINC is ready. */
      write_checkpoint(progress);
      last_checkpoint_progress = progress;
      boinc_checkpoint_completed();
#else
#ifdef _WIN32
        ReleaseSemaphore(checkpoint_semaphoreA,1,NULL);
        WaitForSingleObject(checkpoint_semaphoreB,INFINITE);
#else
        sem_post(&checkpoint_semaphoreA);
        sem_wait(&checkpoint_semaphoreB);
#endif
#endif
#ifdef TRACE
        printf("Thread %d: Continuing after checkpoint\n",th);
#endif
      }
#ifdef TRACE
      else {
        printf("Thread %d: Not checkpointing.\n",th);
	  }
#endif
#ifdef SINGLE_THREAD
    }
#endif
    }

    app_thread_fun1(th,P,K,cthread_count,plen);
    thread_data[th].count = count;
    thread_data[th].sum = sum;
    app_thread_fini(th);

    if (p0 >= pmax)
      fprintf(stderr,"\n%sThread %d completed", bmprefix(),th);
    else
    fprintf(stderr,"\n%sThread %d interrupted", bmprefix(),th);

    /* Just in case a checkpoint is signalled before other threads exit */
    no_more_checkpoints = 1;
  }
#ifndef SINGLE_THREAD
#ifdef _WIN32
  ReleaseSemaphore(checkpoint_semaphoreA,1,NULL);
#else
  sem_post(&checkpoint_semaphoreA);
#endif

#ifndef _WIN32
  pthread_cleanup_pop(1);
#endif
#endif

  if(cthread_count > 0) free(P0);
  if(cthread_count > 0) free(K0);
  return 0;
}

#ifndef SINGLE_THREAD
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
#endif


int main(int argc, char *argv[])
{
#ifndef SINGLE_THREAD
#ifdef _WIN32
  HANDLE tid[MAX_THREADS];
  DWORD thread_ret;
#else
  pthread_t tid[MAX_THREADS];
  void *thread_ret;
  int joined;
#endif
#endif
  uint64_t pstop;
  int th, process_ret = EXIT_SUCCESS;

  program_start_time = elapsed_usec();
  app_banner();

  read_config_file(CONFIG_FILENAME);
  process_args(argc,argv);
#ifdef USE_BOINC
  boincordie(do_boinc_init(), "BOINC initialization failed!\n");
#endif
#ifdef USE_BOINC
  // Get the priority.
  //priority_opt = getThreadPriority();
  //printf("Got priority %d.\n", priority_opt);
#endif

  if (pmin < PMIN_MIN)
    pmin = PMIN_MIN;
  if (pmax > PMAX_MAX)
    pmax = PMAX_MAX;
  if (pmin >= pmax)
  {
    if (pmax == 0 && pmin < PMAX_MAX-1000000000)
    {
      fprintf(stderr,"%spmax not specified, using default pmax = pmin + 1e9\n", bmprefix());
      pmax = pmin + 1000000000; /* Default range */
    }
    else
    {
      fprintf(stderr, "%sOption out of range: pmax must be greater than pmin\n", bmprefix());
      boinc_finish(EXIT_FAILURE);
    }
  }
  if(pmin_str == NULL) pmin_str = empty_string;
  asprintf(&checkpoint_filename, CHECKPOINT_FILENAME, pmin_str);

  checkpoint_period = (uint64_t)checkpoint_opt * 1000000; /* usec */
  report_period = (uint64_t)report_opt * 1000000; /* usec */

  if (blocksize_opt/chunksize_opt < (unsigned int)num_threads)
  {
    chunksize_opt = blocksize_opt/num_threads;
    if (chunksize_opt < CHUNKSIZE_OPT_MIN)
    {
      chunksize_opt = CHUNKSIZE_OPT_MIN;
      blocksize_opt = CHUNKSIZE_OPT_MIN*num_threads;
    }
  }

#ifdef USE_BOINC
#ifdef _WIN32
  SetPriorityClass(GetCurrentProcess(),NORMAL_PRIORITY_CLASS);
#else
  setpriority(PRIO_PROCESS,0,0);
#endif
#else
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
#endif

  if ((uint64_t)qmax*qmax > pmax)
    qmax = (unsigned int)sqrt((double)pmax);
  init_sieve_primes(qmax);

  app_init();

  fprintf(stderr,"%sSieve started: %"PRIu64" <= p < %"PRIu64"\n", bmprefix(),pmin,pmax);

  pstart = read_checkpoint();
  pstart |= 1; /* Must be odd! */

  sv = create_sieve(pstart,pmax,qmax,chunksize_opt,blocksize_opt,blocks_opt);

  init_signals();

#ifdef SINGLE_THREAD
  last_checkpoint_time = sieve_start_time;
  last_checkpoint_progress = pstart;
  last_report_time = sieve_start_time;
  last_report_processor_time = sieve_start_processor_time;
  last_report_processor_cycles = processor_cycles();
  last_report_progress = pstart;
  sieve_start_date = time(NULL);
  sieve_start_time = elapsed_usec();
  sieve_start_processor_time = processor_usec();
  thread_fun((void *)0);
  fini_signals();
#else
#ifdef _WIN32
  checkpoint_semaphoreA = CreateSemaphore(NULL,0,2147483647,NULL);
  checkpoint_semaphoreB = CreateSemaphore(NULL,0,2147483647,NULL);
#else
  pthread_mutex_init(&exiting_mutex,NULL);
  pthread_cond_init(&exiting_cond,NULL);
  sem_init(&checkpoint_semaphoreA,0,0);
  sem_init(&checkpoint_semaphoreB,0,0);
#endif

  sieve_start_date = time(NULL);
  sieve_start_time = elapsed_usec();
  sieve_start_processor_time = processor_usec();

  /* Start child threads */
#ifdef _WIN32
  for (th = 0; th < num_threads; th++)
    if ((tid[th] = (HANDLE)
         _beginthreadex(NULL,0,thread_fun_wrapper,(void *)th,0,NULL)) == 0)
    {
      perror("_beginthreadex");
      boinc_finish(EXIT_FAILURE);
    }
#else
  pthread_mutex_lock(&exiting_mutex);
  {
  long thl;
  for (thl = 0; thl < num_threads; thl++)
    if (pthread_create(&tid[thl],NULL,thread_fun,(void *)thl) != 0)
    {
      perror("pthread_create");
      boinc_finish(EXIT_FAILURE);
    }
  }
#endif

  last_checkpoint_time = sieve_start_time;
  last_checkpoint_progress = pstart;
  last_report_time = sieve_start_time;
  last_report_processor_time = sieve_start_processor_time;
  last_report_processor_cycles = processor_cycles();
  last_report_progress = pstart;

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
        if(boinc_time_to_checkpoint()) {
          write_checkpoint(progress);
          last_checkpoint_progress = progress;
          boinc_checkpoint_completed();
        }
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
    {
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

  fprintf(stderr,"\n%sWaiting for threads to exit", bmprefix());

#ifdef _WIN32
  /* Wait for all threads, then examine return values and close. */
  WaitForMultipleObjects(num_threads,tid,1,INFINITE);
  for (th = 0; th < num_threads; th++)
  {
    if (GetExitCodeThread(tid[th],&thread_ret) && thread_ret != 0)
    {
      fprintf(stderr,"%sThread %d failed: %lX\n", bmprefix(),th,thread_ret);
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
        fprintf(stderr,"%sThread %d failed: %p\n", bmprefix(),th,thread_ret);
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
        fprintf(stderr,"%sThread %d failed: %p\n", bmprefix(),th,thread_ret);
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
#endif
  if (process_ret == EXIT_SUCCESS)
  {
    pstop = next_chunk(sv);
    write_checkpoint(pstop);
  }
  else
  {
    pstop = last_checkpoint_progress;
  }

  if (pstop >= pmax)
  {
    fprintf(stderr,"\n%sSieve complete: %"PRIu64" <= p < %"PRIu64"\n", bmprefix(),pmin,pmax);
    remove(checkpoint_filename);
  }
  else
  {
   fprintf(stderr,"\n%sSieve incomplete: %"PRIu64" <= p < %"PRIu64"\n", bmprefix(),pmin,pstop);
  }

  app_fini(pstop);

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
    fprintf(stderr,"%scount=%"PRIu64",sum=0x%016"PRIx64"\n", bmprefix(),cand_count,cand_sum);
  }

  /* Print statistics for this run */
  {
    uint64_t stop_time = elapsed_usec();
    uint64_t stop_processor_time = processor_usec();
    fprintf(stderr,"%sElapsed time: %.2f sec. (%.2f init + %.2f sieve)"
        " at %.0f p/sec.\n", bmprefix(),
        (stop_time-program_start_time)/1000000.0,
        (sieve_start_time-program_start_time)/1000000.0,
        (stop_time-sieve_start_time)/1000000.0,
        (double)(pstop-pstart)/(stop_time-sieve_start_time)*1000000);
    fprintf(stderr,"%sProcessor time: %.2f sec. (%.2f init + %.2f sieve)"
        " at %.0f p/sec.\n", bmprefix(),
        (stop_processor_time)/1000000.0,
        (sieve_start_processor_time)/1000000.0,
        (stop_processor_time-sieve_start_processor_time)/1000000.0,
        (double)(pstop-pstart)/(stop_processor_time-sieve_start_processor_time)*1000000);
    fprintf(stderr,"%sAverage processor utilization: %.2f (init), %.2f (sieve)\n",
        bmprefix(),
        (double)(sieve_start_processor_time)
        /(sieve_start_time-program_start_time),
        (double)(stop_processor_time-sieve_start_processor_time)
        /(stop_time-sieve_start_time));
  }

  if(got_sigint) raise(SIGINT);
  if(got_sigterm) raise(SIGTERM);
#ifdef SIGHUP
  if(got_sighup) raise(SIGHUP);
#endif
  boinc_finish(process_ret);
  return process_ret;
}
