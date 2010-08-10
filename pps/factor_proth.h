/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* factor_proth.h -- (C) Ken Brazier August 2010.
   Factor a Proth number with small primes, and see if it breaks.
   To be used to test whether potential larger factors are useful.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _FACTOR_PROTH
#define _FACTOR_PROTH

void sieve_small_primes(int min);

int try_all_factors(uint64_t K, unsigned int N, int sign);

#endif
