/* sieve.h -- (C) Geoffrey Reynolds, April 2009.

   Simple multi-threaded sieve for prime (or probable prime) candidates:
     p,        where 3 <= p < 2^63, p odd
     k*2^n+1,  where 1 <= k < 2^63, k odd, 0 < n < 2^31.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SIEVE_H
#define _SIEVE_H

#include <limits.h>
#include <stdint.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#ifndef ULONG_BIT
# if (ULONG_MAX == UINT_MAX)
#  define ULONG_BIT 32
# else
#  define ULONG_BIT 64
# endif
#endif

#ifndef SMALL_PRIMES
# if (ULONG_BIT==32)
#  define SMALL_PRIMES 10
# else
#  define SMALL_PRIMES 17
# endif
#endif

#define MAX_SIEVE_BLOCKS 4

typedef struct
{
  uint64_t cand_next;         /* start of next block */
  uint64_t cand_max;          /* pmax/kmax */
  unsigned int res_max;
  unsigned int res_num;
#ifdef __GNUC__
  unsigned int res[0];
#else
  unsigned int res[1];
#endif
} sieve_data_t;

typedef struct
{
#ifdef _WIN32
  CRITICAL_SECTION mutexA;
  CRITICAL_SECTION mutexB;
  HANDLE condC;
  unsigned int condC_waiting;
#else
  pthread_mutex_t mutexA;
  pthread_mutex_t mutexB;
  pthread_cond_t condC;
#endif

  uint64_t cand_next;         /* next chunk to be returned by get_chunk() */
  unsigned char free_blocks;  /* Number of free blocks */
  unsigned char curr_block;   /* Block to try to get next chunk from */
  unsigned char sieve_done;   /* Set if there are no more blocks to sieve */
  unsigned char all_done;     /* Set if there are no more chunks to get */

  unsigned int chunk_size;    /* number of unsigned longs per chunk */
  unsigned int block_size;    /* number of unsigned longs per block */
  unsigned int chunk_bits;    /* number of bits per chunk */
  unsigned int block_bits;    /* number of bits per block */
  unsigned int num_chunks;
  unsigned int num_blocks;

  struct
  {
    unsigned int free_chunks; /* free chunks (= num_chunks if block free) */
    unsigned int next_chunk;  /* next ready chunk (= num_chunks if none) */
    uint64_t base;            /* First candidate */
    unsigned long *bits;      /* bitmap of odd candidates p or k */
  } block[MAX_SIEVE_BLOCKS];

  sieve_data_t sieve_data;

} sieve_t;

/* These variables should be hidden, but they are used in other places. */
extern unsigned int *sieve_primes; /* Master list of odd primes */
extern unsigned int num_sieve_primes; /* Number of primes in list */

/* Initialize/free tables of sieve primes */
void init_sieve_primes(unsigned int qmax);
void free_sieve_primes(void);

/* Create/destroy sieve objects */
sieve_t *create_sieve(uint64_t pmin,
                      uint64_t pmax,
                      unsigned int qmax,
                      unsigned int chunk_bytes,
                      unsigned int block_bytes,
                      unsigned int blocks);

sieve_t *create_gfn_sieve(uint64_t kmin,
                          uint64_t kmax,
                          unsigned int n,
                          unsigned int qmax,
                          unsigned int chunk_bytes,
                          unsigned int block_bytes,
                          unsigned int blocks);

void destroy_sieve(sieve_t *sv);


/* Get/free a chunk of candidates */
#if TRACE
uint64_t get_chunk(int th, sieve_t *sv, unsigned long **chunk);
void free_chunk(int th,sieve_t *sv, uint64_t chunk);
#else
uint64_t get_chunk(sieve_t *sv, unsigned long **chunk);
void free_chunk(sieve_t *sv, uint64_t chunk);
#endif


/* Get the next candidate to be handed out by get_block */
uint64_t next_chunk(sieve_t *sv);


#endif /* _SIEVE_H */
