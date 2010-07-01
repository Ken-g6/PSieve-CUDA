/* appcu.h -- (C) Ken Brazier February 2010.

   Proth Prime Search sieve CUDA portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifndef _APPCU_H
#define _APPCU_H 1

// BLOCKSIZE should be a power of two for greatest efficiency.
#define BLOCKSIZE 128
// Maximum iterations to be run by any one kernel.  Breaking up the kernel this way should improve display effects.
#define ITERATIONS_PER_KERNEL 3000

#ifdef __cplusplus
extern "C" {
#endif

unsigned int cuda_app_init(int gpuno);
//void setup_ps(const uint64_t *P, unsigned int cthread_count);
void check_ns(const uint64_t *P, const unsigned int cthread_count);
void get_factors_found(unsigned char *factor_found, const unsigned int cthread_count);
void cuda_finalize(void);

extern unsigned int ld_nstep;
extern int ld_bbits;
extern uint64_t ld_r0;
#ifdef __cplusplus
}
#endif

#endif
