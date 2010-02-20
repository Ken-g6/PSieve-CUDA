/* appcu.h -- (C) Ken Brazier February 2010.

   Proth Prime Search sieve CUDA portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifndef _APPCU_H
#define _APPCU_H 1

#ifdef __cplusplus
extern "C" {
#endif

unsigned int cuda_app_init(int gpuno);
void check_ns(const uint64_t *P, const uint64_t *K, unsigned char *factor_found, unsigned int cthread_count);
void cuda_finalize(void);
#ifdef __cplusplus
}
#endif

extern unsigned int ld_nstep;
extern unsigned int gpuno;

#endif
