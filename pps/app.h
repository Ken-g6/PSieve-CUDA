/* app.h -- (C) Geoffrey Reynolds, April 2009.
 * and Ken Brazier October 2009 - April 2010.

   Proth Prime Search sieve.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _APP_H
#define _APP_H 1

#include <stdio.h>
#include "stdint.h"

#ifndef USE_OPENCL
#define APP_VERSION "cuda-0.1.5a"
#else
#define APP_VERSION "cl-0.1.0-beta"
#endif

/* Number of primes to buffer between calls to app_thread_fun()
 */
#define APP_BUFLEN 6

#define CHECKPOINT_FILENAME "ppcheck%s.txt"
#define OLD_CHECKPOINT_FILENAME "ppcheckpoint.txt"

#define CONFIG_FILENAME "ppconfig.txt"

#define FACTORS_FILENAME_DEFAULT "ppfactors.txt"

#ifdef USE_OPENCL
#define APP_SHORT_OPTS "k:K:n:N:i:f:b:d:m:Rv:"
#else
#define APP_SHORT_OPTS "k:K:n:N:i:f:b:d:m:R"
#endif
#define APP_LONG_OPTS \
  {"kmin",          required_argument, 0, 'k'}, \
  {"kmax",          required_argument, 0, 'K'}, \
  {"nmin",          required_argument, 0, 'n'}, \
  {"nmax",          required_argument, 0, 'N'}, \
  {"input",         required_argument, 0, 'i'}, \
  {"factors",       required_argument, 0, 'f'}, \
  {"device",        required_argument, 0, 'd'}, \
  {"bitsatatime",   required_argument, 0, 'b'}, \
  {"mthreads",      required_argument, 0, 'm'}, \
  {"vecsize",       required_argument, 0, 'v'}, \
  {"riesel",        no_argument,       0, 'R'},

#ifdef USE_BOINC
#include <error_numbers.h>
#else
// Error codes, matching BOINC's.
#define ERR_FOPEN -108
#define ERR_SCANF -140
#define ERR_INVALID_PARAM -178
#define ERR_NEG -117
#define ERR_INSUFFICIENT_RESOURCE -198
#define ERR_NOT_IMPLEMENTED -150
#endif

void app_banner(void);
int app_parse_option(int opt, char *arg, const char *source);
void app_help(void);
void app_init(void);
unsigned int app_thread_init(int th);
void app_thread_fun(int th, const uint64_t *P, uint64_t *lastP, const unsigned int cthread_count);
//void app_thread_fun(int th, uint64_t *__attribute__((aligned(16))) P, uint64_t *__attribute__((aligned(16))) K);
void app_thread_fun1(int th, uint64_t *P, uint64_t *lastP, const unsigned int cthread_count, unsigned int len);
void app_thread_fini(int th);
int app_read_checkpoint(FILE *f);
void app_write_checkpoint(FILE *f);
void app_fini(void);
#ifdef __cplusplus
extern "C" {
#endif
extern unsigned int nmin, nmax;
extern int search_proth;
extern uint64_t kmax, kmin;
#ifdef __cplusplus
}
#endif
#endif /* _APP_H */
