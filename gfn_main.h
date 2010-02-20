/* gfn_main.h -- (C) Geoffrey Reynolds, May 2009.

   Multithreaded sieve application for algorithms of the form:

   For each odd k in 1 <= k0 <= k < k1 < 2^62
     If k*2^n+1 has no factors < q (n,q fixed)
       Do something with k


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GFN_MAIN_H
#define _GFN_MAIN_H 1

#include <stdint.h>


#define MAX_THREADS 32  /* Excluding parent thread */

#define BLOCKSIZE_OPT_DEFAULT 32768  /* Bytes per sieve block */
#define BLOCKSIZE_OPT_MIN 1024
#define BLOCKSIZE_OPT_MAX (1<<27)  /* 128Mb = 2^30 bits */

#define CHUNKSIZE_OPT_DEFAULT 128  /* Bytes per sieve chunk */
#define CHUNKSIZE_OPT_MIN 16
#define CHUNKSIZE_OPT_MAX (1<<17)  /* 128Kb = 2^20 bits */

#define BLOCKS_OPT_DEFAULT 2       /* Number of sieve blocks */
#define BLOCKS_OPT_MIN 2
#define BLOCKS_OPT_MAX 4


#define KMIN_MIN 1
#define KMAX_MAX (UINT64_C(1)<<62)
#define NMIN_MIN 1
#define NMAX_MAX ((1U<<31)-1)
#define QMAX_MAX (1U<<31)

#define REPORT_OPT_DEFAULT 60  /* seconds between status reports */
#define CHECKPOINT_OPT_DEFAULT 300  /* seconds between checkpoints */


extern unsigned int num_threads; /* Excluding parent thread */
extern uint64_t kmin, kmax;
extern unsigned int nmin, nmax;


#endif /* _GFN_MAIN_H */
