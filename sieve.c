/* sieve.c -- (C) Geoffrey Reynolds, April 2009.

   Simple multi-threaded sieve for prime (or probable prime) candidates:
     p,        where 3 <= p < 2^63, p odd
     k*2^n+1,  where 1 <= k < 2^63, k odd, 0 < n < 2^31.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include "stdint.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inttypes.h"
#include <math.h>
#include "putil.h"
#include "sieve.h"


unsigned int *sieve_primes = NULL; /* Master list of odd primes */
unsigned int num_sieve_primes = 0; /* Number of primes in list */

static unsigned int lim_sieve_primes = 0; /* List contains all primes <= lim */


#if SMALL_PRIMES
#ifdef _WIN32
static void init_residues(sieve_data_t *sv)
#else
static void __attribute__((noinline)) init_residues(sieve_data_t *sv)
#endif
{
  uint64_t p0;
  unsigned int i, n, q, r;

  p0 = sv->cand_next;
  n = sv->res_num;

  for (i = 0; i < n; i++)
  {
    q = sieve_primes[i];
    r = p0 % q;

    if (r != 0)
    {
      r = q-r;
      if (r % 2)
        r += q;
      r /= 2;
    }

    if (p0 <= q)
      r += q;

    sv->res[i] = r;
  }
}
#endif /* SMALL_PRIMES */


/* Create a new sieve object for odd candidates p in pmin <= p < pmax using
   sieve primes no larger than qmax. If qmax = 0 then use the greatest
   available sieve prime.
 */
static sieve_data_t *
create_sieve_data(uint64_t pmin, uint64_t pmax, unsigned int qmax)
{
  sieve_data_t *sv;
  unsigned int num;

  assert(pmin < pmax);
  assert(pmax > 3);
  assert(sieve_primes != NULL);
  assert(num_sieve_primes > 0);

  if (pmin < 2)
    pmin = 2;

  pmin |= 1;
  pmax |= 1;

  if (qmax == 0)
  {
    /* Restrict qmax to 2^31-1, the greatest possible sieve prime */
    if (pmax > (uint64_t)((1U<<31)-1)*((1U<<31)-1))
      qmax = (1U<<31)-1;
    else
      qmax = sqrt(pmax);
  }

#if SMALL_PRIMES
  if (qmax < sieve_primes[SMALL_PRIMES-1])
    qmax = sieve_primes[SMALL_PRIMES-1];
#endif

  /* Find the number of initialized sieve primes <= qmax */
  for (num = num_sieve_primes; num > 0; num--)
    if (sieve_primes[num-1] <= qmax)
      break;

  /* Allocate a new sieve object with room for num residues */
  sv = xmalloc(sizeof(sieve_data_t) + num*sizeof(unsigned int));

  sv->cand_next = pmin;
  sv->cand_max = pmax;
  sv->res_max = num;
#if SMALL_PRIMES
  sv->res_num = SMALL_PRIMES;
  init_residues(sv);
#else
  sv->res_num = 0;
  /* Residues will be initialized during the first call to sieve(). */
#endif

#ifndef NDEBUG
  printf("create_sieve_data: Using %u odd primes, max=%u, lim=%u\n",
         num,sieve_primes[num-1],qmax);
#endif

  return sv;
}

#if SMALL_PRIMES
#if (ULONG_BIT==32)
#define MASK3 (1UL<<0|1UL<<3|1UL<<6|1UL<<9|1UL<<12|1UL<<15|1UL<<18|1UL<<21|1UL<<24|1UL<<27|1UL<<30)
#define MASK5 (1UL<<0|1UL<<5|1UL<<10|1UL<<15|1UL<<20|1UL<<25|1UL<<30)
#define MASK7 (1UL<<0|1UL<<7|1UL<<14|1UL<<21|1UL<<28)
#define MASK11 (1UL<<0|1UL<<11|1UL<<22)
#define MASK13 (1UL<<0|1UL<<13|1UL<<26)
#define MASK17 (1UL<<0|1UL<<17)
#define MASK19 (1UL<<0|1UL<<19)
#define MASK23 (1UL<<0|1UL<<23)
#define MASK29 (1UL<<0|1UL<<29)
#define MASK31 (1UL<<0|1UL<<31)
#else
#define MASK3 (1UL<<0|1UL<<3|1UL<<6|1UL<<9|1UL<<12|1UL<<15|1UL<<18|1UL<<21|1UL<<24|1UL<<27|1UL<<30|1UL<<33|1UL<<36|1UL<<39|1UL<<42|1UL<<45|1UL<<48|1UL<<51|1UL<<54|1UL<<57|1UL<<60|1UL<<63)
#define MASK5 (1UL<<0|1UL<<5|1UL<<10|1UL<<15|1UL<<20|1UL<<25|1UL<<30|1UL<<35|1UL<<40|1UL<<45|1UL<<50|1UL<<55|1UL<<60)
#define MASK7 (1UL<<0|1UL<<7|1UL<<14|1UL<<21|1UL<<28|1UL<<35|1UL<<42|1UL<<49|1UL<<56|1UL<<63)
#define MASK11 (1UL<<0|1UL<<11|1UL<<22|1UL<<33|1UL<<44|1UL<<55)
#define MASK13 (1UL<<0|1UL<<13|1UL<<26|1UL<<39|1UL<<52)
#define MASK17 (1UL<<0|1UL<<17|1UL<<34|1UL<<51)
#define MASK19 (1UL<<0|1UL<<19|1UL<<38|1UL<<57)
#define MASK23 (1UL<<0|1UL<<23|1UL<<46)
#define MASK29 (1UL<<0|1UL<<29|1UL<<58)
#define MASK31 (1UL<<0|1UL<<31|1UL<<62)
#define MASK37 (1UL<<0|1UL<<37)
#define MASK41 (1UL<<0|1UL<<41)
#define MASK43 (1UL<<0|1UL<<43)
#define MASK47 (1UL<<0|1UL<<47)
#define MASK53 (1UL<<0|1UL<<53)
#define MASK59 (1UL<<0|1UL<<59)
#define MASK61 (1UL<<0|1UL<<61)
#endif
#endif

/* Sieve the next block of len unsigned longs in sv.
   Stores the start of the block in *pp, bitmap of candidates in buf.
   Returns the number of elements of buf that were used (normally len).
   If there are no more blocks to sieve then the sieve end is stored in *pp
   and zero returned.
*/
static unsigned int
sieve(sieve_data_t *sv, uint64_t *pp, unsigned long buf[], unsigned int len)
{
  uint64_t p0, p1;
  unsigned int i, bits, q, r, rnum, rmax;

  assert(sv != NULL);
  assert(pp != NULL);
  assert(buf != NULL);
  assert(len > 0);
  assert(len <= (1U<<31)/ULONG_BIT);

  bits = len*ULONG_BIT;
  p0 = sv->cand_next;
  if (sv->cand_max - p0 >= (uint64_t)bits*2)
    p1 = p0 + (uint64_t)bits*2;
  else
  {
    p1 = sv->cand_max;
    bits = (p1-p0)/2;
    len = (bits+ULONG_BIT-1)/ULONG_BIT;
  }
  sv->cand_next = p1;

#if SMALL_PRIMES
  assert(SMALL_PRIMES <= 17 && (ULONG_BIT==64 || SMALL_PRIMES <= 10));
  {
    unsigned long x;
    int r3 = sv->res[0];
#if (SMALL_PRIMES >= 2)
    int r5 = sv->res[1];
#endif
#if (SMALL_PRIMES >= 3)
    int r7 = sv->res[2];
#endif
#if (SMALL_PRIMES >= 4)
    int r11 = sv->res[3];
#endif
#if (SMALL_PRIMES >= 5)
    int r13 = sv->res[4];
#endif
#if (SMALL_PRIMES >= 6)
    int r17 = sv->res[5];
#endif
#if (SMALL_PRIMES >= 7)
    int r19 = sv->res[6];
#endif
#if (SMALL_PRIMES >= 8)
    int r23 = sv->res[7];
#endif
#if (SMALL_PRIMES >= 9)
    int r29 = sv->res[8];
#endif
#if (SMALL_PRIMES >= 10)
    int r31 = sv->res[9];
#endif
#if (SMALL_PRIMES >= 11)
    int r37 = sv->res[10];
#endif
#if (SMALL_PRIMES >= 12)
    int r41 = sv->res[11];
#endif
#if (SMALL_PRIMES >= 13)
    int r43 = sv->res[12];
#endif
#if (SMALL_PRIMES >= 14)
    int r47 = sv->res[13];
#endif
#if (SMALL_PRIMES >= 15)
    int r53 = sv->res[14];
#endif
#if (SMALL_PRIMES >= 16)
    int r59 = sv->res[15];
#endif
#if (SMALL_PRIMES >= 17)
    int r61 = sv->res[16];
#endif

    /* Allow for the situation where residue r => prime q */
    x = MASK3<<r3;
    if (r3 >= 3)
      r3 -= 3;
    if ((r3 -= ULONG_BIT%3) < 0)
      r3 += 3;
#if (SMALL_PRIMES >= 2)
    x |= MASK5<<r5;
    if (r5 >= 5)
      r5 -= 5;
    if ((r5 -= ULONG_BIT%5) < 0)
      r5 += 5;
#endif
#if (SMALL_PRIMES >= 3)
    x |= MASK7<<r7;
    if (r7 >= 7)
      r7 -= 7;
    if ((r7 -= ULONG_BIT%7) < 0)
      r7 += 7;
#endif
#if (SMALL_PRIMES >= 4)
    x |= MASK11<<r11;
    if (r11 >= 11)
      r11 -= 11;
    if ((r11 -= ULONG_BIT%11) < 0)
      r11 += 11;
#endif
#if (SMALL_PRIMES >= 5)
    x |= MASK13<<r13;
    if (r13 >= 13)
      r13 -= 13;
    if ((r13 -= ULONG_BIT%13) < 0)
      r13 += 13;
#endif
#if (SMALL_PRIMES >= 6)
#if (ULONG_BIT==32)
    if (r17 >= 17)
      r17 -= 17, x |= (MASK17 & ~1UL)<<r17;
    else
      x |= MASK17<<r17;
#else
    x |= MASK17<<r17;
    if (r17 >= 17)
      r17 -= 17;
#endif
    if ((r17 -= ULONG_BIT%17) < 0)
      r17 += 17;
#endif
#if (SMALL_PRIMES >= 7)
#if (ULONG_BIT==32)
    if (r19 >= 19)
      r19 -= 19, x |= (MASK19 & ~1UL)<<r19;
    else
      x |= MASK19<<r19;
#else
    x |= MASK19<<r19;
    if (r19 >= 19)
      r19 -= 19;
#endif
    if ((r19 -= ULONG_BIT%19) < 0)
      r19 += 19;
#endif
#if (SMALL_PRIMES >= 8)
#if (ULONG_BIT==32)
    if (r23 >= 23)
      r23 -= 23, x |= (MASK23 & ~1UL)<<r23;
    else
      x |= MASK23<<r23;
#else
    x |= MASK23<<r23;
    if (r23 >= 23)
      r23 -= 23;
#endif
    if ((r23 -= ULONG_BIT%23) < 0)
      r23 += 23;
#endif
#if (SMALL_PRIMES >= 9)
#if (ULONG_BIT==32)
    if (r29 >= 29)
      r29 -= 29, x |= (MASK29 & ~1UL)<<r29;
    else
      x |= MASK29<<r29;
#else
    x |= MASK29<<r29;
    if (r29 >= 29)
      r29 -= 29;
#endif
    if ((r29 -= ULONG_BIT%29) < 0)
      r29 += 29;
#endif
#if (SMALL_PRIMES >= 10)
#if (ULONG_BIT==32)
    if (r31 >= 31)
      r31 -= 31, x |= (MASK31 & ~1UL)<<r31;
    else
      x |= MASK31<<r31;
#else
    x |= MASK31<<r31;
    if (r31 >= 31)
      r31 -= 31;
#endif
    if ((r31 -= ULONG_BIT%31) < 0)
      r31 += 31;
#endif
#if (SMALL_PRIMES >= 11)
    if (r37 >= 37)
      r37 -= 37, x |= (MASK37 & ~1UL)<<r37;
    else
      x |= MASK37<<r37;
    if ((r37 -= ULONG_BIT%37) < 0)
      r37 += 37;
#endif
#if (SMALL_PRIMES >= 12)
    if (r41 >= 41)
      r41 -= 41, x |= (MASK41 & ~1UL)<<r41;
    else
      x |= MASK41<<r41;
    if ((r41 -= ULONG_BIT%41) < 0)
      r41 += 41;
#endif
#if (SMALL_PRIMES >= 13)
    if (r43 >= 43)
      r43 -= 43, x |= (MASK43 & ~1UL)<<r43;
    else
      x |= MASK43<<r43;
    if ((r43 -= ULONG_BIT%43) < 0)
      r43 += 43;
#endif
#if (SMALL_PRIMES >= 14)
    if (r47 >= 47)
      r47 -= 47, x |= (MASK47 & ~1UL)<<r47;
    else
      x |= MASK47<<r47;
    if ((r47 -= ULONG_BIT%47) < 0)
      r47 += 47;
#endif
#if (SMALL_PRIMES >= 15)
    if (r53 >= 53)
      r53 -= 53, x |= (MASK53 & ~1UL)<<r53;
    else
      x |= MASK53<<r53;
    if ((r53 -= ULONG_BIT%53) < 0)
      r53 += 53;
#endif
#if (SMALL_PRIMES >= 16)
    if (r59 >= 59)
      r59 -= 59, x |= (MASK59 & ~1UL)<<r59;
    else
      x |= MASK59<<r59;
    if ((r59 -= ULONG_BIT%59) < 0)
      r59 += 59;
#endif
#if (SMALL_PRIMES >= 17)
    if (r61 >= 61)
      r61 -= 61, x |= (MASK61 & ~1UL)<<r61;
    else
      x |= MASK61<<r61;
    if ((r61 -= ULONG_BIT%61) < 0)
      r61 += 61;
#endif
    buf[0] = ~x;

    for (i = 1; i < len; i++)
    {
      x = MASK3<<r3;
      if ((r3 -= ULONG_BIT%3) < 0)
        r3 += 3;
#if (SMALL_PRIMES >= 2)
      x |= MASK5<<r5;
      if ((r5 -= ULONG_BIT%5) < 0)
        r5 += 5;
#endif
#if (SMALL_PRIMES >= 3)
      x |= MASK7<<r7;
      if ((r7 -= ULONG_BIT%7) < 0)
        r7 += 7;
#endif
#if (SMALL_PRIMES >= 4)
      x |= MASK11<<r11;
      if ((r11 -= ULONG_BIT%11) < 0)
        r11 += 11;
#endif
#if (SMALL_PRIMES >= 5)
      x |= MASK13<<r13;
      if ((r13 -= ULONG_BIT%13) < 0)
        r13 += 13;
#endif
#if (SMALL_PRIMES >= 6)
      x |= MASK17<<r17;
      if ((r17 -= ULONG_BIT%17) < 0)
        r17 += 17;
#endif
#if (SMALL_PRIMES >= 7)
      x |= MASK19<<r19;
      if ((r19 -= ULONG_BIT%19) < 0)
        r19 += 19;
#endif
#if (SMALL_PRIMES >= 8)
      x |= MASK23<<r23;
      if ((r23 -= ULONG_BIT%23) < 0)
        r23 += 23;
#endif
#if (SMALL_PRIMES >= 9)
      x |= MASK29<<r29;
      if ((r29 -= ULONG_BIT%29) < 0)
        r29 += 29;
#endif
#if (SMALL_PRIMES >= 10)
      x |= MASK31<<r31;
      if ((r31 -= ULONG_BIT%31) < 0)
        r31 += 31;
#endif
      buf[i] = ~x;
    }

#if (SMALL_PRIMES >= 11)
    for (i = 1; i < len; i++)
    {
      x = MASK37<<r37;
      if ((r37 -= ULONG_BIT%37) < 0)
        r37 += 37;
#if (SMALL_PRIMES >= 12)
      x |= MASK41<<r41;
      if ((r41 -= ULONG_BIT%41) < 0)
        r41 += 41;
#endif
#if (SMALL_PRIMES >= 13)
      x |= MASK43<<r43;
      if ((r43 -= ULONG_BIT%43) < 0)
        r43 += 43;
#endif
#if (SMALL_PRIMES >= 14)
      x |= MASK47<<r47;
      if ((r47 -= ULONG_BIT%47) < 0)
        r47 += 47;
#endif
#if (SMALL_PRIMES >= 15)
      x |= MASK53<<r53;
      if ((r53 -= ULONG_BIT%53) < 0)
        r53 += 53;
#endif
#if (SMALL_PRIMES >= 16)
      x |= MASK59<<r59;
      if ((r59 -= ULONG_BIT%59) < 0)
        r59 += 59;
#endif
#if (SMALL_PRIMES >= 17)
      x |= MASK61<<r61;
      if ((r61 -= ULONG_BIT%61) < 0)
        r61 += 61;
#endif

      buf[i] &= ~x;
    }
#endif

    sv->res[0] = r3;
#if (SMALL_PRIMES >= 2)
    sv->res[1] = r5;
#endif
#if (SMALL_PRIMES >= 3)
    sv->res[2] = r7;
#endif
#if (SMALL_PRIMES >= 4)
    sv->res[3] = r11;
#endif
#if (SMALL_PRIMES >= 5)
    sv->res[4] = r13;
#endif
#if (SMALL_PRIMES >= 6)
    sv->res[5] = r17;
#endif
#if (SMALL_PRIMES >= 7)
    sv->res[6] = r19;
#endif
#if (SMALL_PRIMES >= 8)
    sv->res[7] = r23;
#endif
#if (SMALL_PRIMES >= 9)
    sv->res[8] = r29;
#endif
#if (SMALL_PRIMES >= 10)
    sv->res[9] = r31;
#endif
#if (SMALL_PRIMES >= 11)
    sv->res[10] = r37;
#endif
#if (SMALL_PRIMES >= 12)
    sv->res[11] = r41;
#endif
#if (SMALL_PRIMES >= 13)
    sv->res[12] = r43;
#endif
#if (SMALL_PRIMES >= 14)
    sv->res[13] = r47;
#endif
#if (SMALL_PRIMES >= 15)
    sv->res[14] = r53;
#endif
#if (SMALL_PRIMES >= 16)
    sv->res[15] = r59;
#endif
#if (SMALL_PRIMES >= 17)
    sv->res[16] = r61;
#endif
    i = SMALL_PRIMES;
  }
#else /* !SMALL_PRIMES */
  /* Set all bits of buf */
  memset(buf,-1,len*sizeof(unsigned long));
  i = 0;
#endif

  rnum = sv->res_num;
  for ( ; i < rnum; i++)
  {
    q = sieve_primes[i];
    
    for (r = sv->res[i]; r < bits; r += q)
      buf[r/ULONG_BIT] &= ~(1UL << r%ULONG_BIT);

    sv->res[i] = r-bits;
  }

  /* Extend residues if necessary. Never necessary for the GFN sieve. */
  rmax = sv->res_max;
  for (i = rnum; i < rmax; i++)
  {
    q = sieve_primes[i];

    if ((uint64_t)q*q >= p1)
      break;

    r = p0 % q;

    if (r != 0)
    {
      r = q-r;
      if (r % 2)
        r += q;
      r /= 2;
    }

    if (p0 > q)
      ;
    else
      r += q;

    for ( ; r < bits; r += q)
      buf[r/ULONG_BIT] &= ~(1UL << r%ULONG_BIT);

    sv->res[i] = r-bits;
  }
  sv->res_num = i;

  *pp = p0;

  /* Clear any bits corresponding to candidates >= p1 */
  if (len > 0 && bits%ULONG_BIT)
    buf[len-1] &= ((1UL << bits%ULONG_BIT)-1);

  return len;
}

/* Extend the list of sieve primes to include q.
 */
static void extend_sieve_primes(unsigned int q)
{
  sieve_data_t *sv;
  uint64_t p0;
  unsigned long buf[2048];
  unsigned int i, j, count, len;

  assert(q <= (1U<<31));

  if (q % 2)
    q++;

  if (q <= lim_sieve_primes)
    return;

  sv = create_sieve_data(lim_sieve_primes,q,0);
  for (count = 0; (len = sieve(sv,&p0,buf,sizeof(buf)/sizeof(unsigned long)));)
    for (i = 0; i < len; i++)
    {
#ifdef __GNUC__
      count += __builtin_popcountl(buf[i]);
#else
      unsigned int k;
      for (k = 0; k < ULONG_BIT; k++)
        if (buf[i] & (1UL<<k))
          count++;
#endif
    }

  free(sv);
  sieve_primes =
    xrealloc(sieve_primes,(num_sieve_primes+count+1)*sizeof(unsigned int));

  sv = create_sieve_data(lim_sieve_primes,q,0);
  for (j = num_sieve_primes; (len = sieve(sv,&p0,buf,2048)); )
    for (i = 0; i < len; i++)
    {
      unsigned int k;
      unsigned long u;

      for (k = 0, u = buf[i]; u != 0; u >>= 1, k++)
      {
#ifdef __GNUC__
        int h = __builtin_ctzl(u);
        u >>= h;
        k += h;
#else
        if ((u & 1UL) == 0)
          continue;
#endif
        sieve_primes[j++] = p0+2*(uint64_t)(i*ULONG_BIT+k);
      }
    }
  assert(j == num_sieve_primes+count);

  free(sv);
  sieve_primes[j] = -1U;
  num_sieve_primes = j;
  lim_sieve_primes = q;
}

/* Generate the list of sieve primes from 3 to at least q.
 */
void init_sieve_primes(unsigned int q)
{
  if (sieve_primes == NULL)
  {
    /* Generate all odd primes up to sqrt(2^31) */

    const unsigned char primes8[] = {
      3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
      71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
      149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211 };

    unsigned long buf[(46339-3)/(2*ULONG_BIT)];
    unsigned int i, j, q, r;

    memset(buf,-1,sizeof(buf));
    for (i = 0; i < sizeof(primes8); i++)
    {
      q = primes8[i];
      for (r = (q-3)/2+q; r < (46339-3)/2; r += q)
        buf[r/ULONG_BIT] &= ~(1UL << r%ULONG_BIT);
    }

    sieve_primes = xmalloc((4791+1)*sizeof(unsigned int));
    for (i = 0, j = 0; i < (46339-3)/2; i++)
      if (buf[i/ULONG_BIT] & (1UL << i%ULONG_BIT))
        sieve_primes[j++] = 3+2*i;
    assert(j == 4791);
    sieve_primes[4791] = -1U;
    num_sieve_primes = 4791;
    lim_sieve_primes = 46348;
  }

  extend_sieve_primes(q);

#ifndef NDEBUG
  printf("init_sieve_primes: Generated %u odd primes, max=%u, lim=%u\n",
         num_sieve_primes,sieve_primes[num_sieve_primes-1],lim_sieve_primes);
#endif
}

void free_sieve_primes(void)
{
  if (sieve_primes != NULL)
  {
    free(sieve_primes);
    sieve_primes = NULL;
  }
}


/* Create a new sieve object for odd candidates p in pmin <= p < pmax using
   sieve primes no larger than qmax. If qmax = 0 then use the greatest
   available sieve prime.
 */
sieve_t *create_sieve(uint64_t pmin,
                      uint64_t pmax,
                      unsigned int qmax,
                      unsigned int chunk_bytes,
                      unsigned int block_bytes,
                      unsigned int blocks)
{
  sieve_t *sv;
  unsigned int i, num, chunks;

  assert(pmin < pmax);
  assert(pmax > 3);
  assert(sieve_primes != NULL);
  assert(num_sieve_primes > 0);

  if (pmin < 2)
    pmin = 2;

  pmin |= 1;
  pmax |= 1;

  if (qmax == 0)
  {
    /* Restrict qmax to 2^31-1, the greatest possible sieve prime */
    if (pmax > (uint64_t)((1U<<31)-1)*((1U<<31)-1))
      qmax = (1U<<31)-1;
    else
      qmax = sqrt(pmax);
  }

#if SMALL_PRIMES
  if (qmax < sieve_primes[SMALL_PRIMES-1])
    qmax = sieve_primes[SMALL_PRIMES-1];
#endif

  /* Find the number of initialized sieve primes <= qmax */
  for (num = num_sieve_primes; num > 0; num--)
    if (sieve_primes[num-1] <= qmax)
      break;

  /* Allocate a new sieve object with room for num residues */
  sv = xmalloc(sizeof(sieve_t) + num*sizeof(unsigned int));

  sv->sieve_data.cand_next = pmin;
  sv->sieve_data.cand_max = pmax;
  sv->sieve_data.res_max = num;
#if SMALL_PRIMES
  sv->sieve_data.res_num = SMALL_PRIMES;
  init_residues(&sv->sieve_data);
#else
  sv->sieve_data.res_num = 0;
  /* Residues will be initialized during the first call to sieve(). */
#endif

#ifdef _WIN32
  InitializeCriticalSection(&sv->mutexA);
  InitializeCriticalSection(&sv->mutexB);
  sv->condC = CreateSemaphore(NULL,0,2147483647,NULL);
  sv->condC_waiting = 0;
#else
  pthread_mutex_init(&sv->mutexA,NULL);
  pthread_mutex_init(&sv->mutexB,NULL);
  pthread_cond_init(&sv->condC,NULL);
#endif

  if (blocks < 2)
    blocks = 2;
  else if (blocks > MAX_SIEVE_BLOCKS)
    blocks = MAX_SIEVE_BLOCKS;

  sv->cand_next = pmin;
  sv->free_blocks = blocks; /* Begin with all blocks free */
  sv->curr_block = blocks;  /* No blocks have ready chunks */
  sv->sieve_done = 0;
  sv->all_done = 0;

  block_bytes = block_bytes/sizeof(unsigned long)*sizeof(unsigned long);
  chunk_bytes = chunk_bytes/sizeof(unsigned long)*sizeof(unsigned long);

  if (chunk_bytes < 2*sizeof(unsigned long))
    chunk_bytes = 2*sizeof(unsigned long);

  if (block_bytes/2 < chunk_bytes)
    block_bytes = 2*chunk_bytes;
  else if (block_bytes > 1<<28) /* 2^31 bits */
    block_bytes = 1<<28;

  if (block_bytes/2 < chunk_bytes)
    chunks = 2;
  else
    chunks = block_bytes/chunk_bytes;

  sv->chunk_size = chunk_bytes/sizeof(unsigned long);
  sv->block_size = chunk_bytes/sizeof(unsigned long)*chunks;
  sv->chunk_bits = sv->chunk_size*ULONG_BIT;
  sv->block_bits = sv->block_size*ULONG_BIT;
  sv->num_chunks = chunks;
  sv->num_blocks = blocks;

  for (i = 0; i < blocks; i++)
  {
    sv->block[i].free_chunks = chunks;
    sv->block[i].next_chunk = chunks;
    sv->block[i].base = 0;
    /* TODO: Use aligned malloc */
    sv->block[i].bits = xmalloc(sv->block_size * sizeof(unsigned long));
  }

#ifndef NDEBUG
  printf("create_sieve: Using %u odd primes, max=%u, lim=%u\n",
         num,sieve_primes[num-1],qmax);
  printf("chunk_size=%u (%ub), block_size=%u (%ub), chunks=%u, blocks=%u\n",
         sv->chunk_size, sv->chunk_size*sizeof(unsigned long),
         sv->block_size, sv->block_size*sizeof(unsigned long),
         sv->num_chunks, sv->num_blocks);
#endif

  return sv;
}

/* Return b^n mod p. Doesn't have to be fast.
 */
static unsigned int powmod32(unsigned int b, unsigned int n, unsigned int p)
{
  unsigned int a;

  for (a = 1; n > 0; n >>= 1)
  {
    if (n & 1)
      a = (uint64_t)a*b%p;
    b = (uint64_t)b*b%p;
  }

  return a;
}

/* Return a-b mod p. Doesn't have to be fast.
 */
static unsigned int submod32(unsigned int a, unsigned int b, unsigned int p)
{
  if (a >= b)
    return a-b;
  else
    return a-b+p;
}

sieve_t *create_gfn_sieve(uint64_t kmin,
                          uint64_t kmax,
                          unsigned int n,
                          unsigned int qmax,
                          unsigned int chunk_bytes,
                          unsigned int block_bytes,
                          unsigned int blocks)
{
  sieve_t *sv;
  unsigned int i, num, chunks;

  assert(kmin < kmax);
  assert(kmax < UINT64_C(1)<<63);
  assert(n > 1);
  assert(sieve_primes != NULL);
  assert(num_sieve_primes > 0);

  kmin |= 1;
  kmax |= 1;

  if (qmax == 0)
    qmax = sieve_primes[num_sieve_primes-1];

#if SMALL_PRIMES
  if (qmax < sieve_primes[SMALL_PRIMES-1])
    qmax = sieve_primes[SMALL_PRIMES-1];
#endif

  /* Find the number of initialized sieve primes <= qmax */
  for (num = num_sieve_primes; num > 0; num--)
    if (sieve_primes[num-1] <= qmax)
      break;

  /* Allocate a new sieve object with room for num residues */
  sv = xmalloc(sizeof(sieve_t) + num*sizeof(unsigned int));

  sv->sieve_data.cand_next = kmin;
  sv->sieve_data.cand_max = kmax;
  sv->sieve_data.res_max = num;
  sv->sieve_data.res_num = num;

  /* Initialize all residues: For each prime q, find the smallest even r
     such that q divides (kmin+r)*2^n+1.
  */
  for (i = 0; i < num; i++)
  {
    unsigned int q, r;

    q = sieve_primes[i];
    r = powmod32((q+1)/2,n,q);  /* 1/2^n mod q, nonzero since q is odd */
    r = q - r;                  /* -1/2^n mod q */
    r = submod32(r,kmin%q,q);   /* so (kmin+r)*2^n+1 = 0 (mod q) */

    if (r % 2)
      r += q;
    r /= 2;

    /* If kmin*2^n+1 <= q then skip ahead to avoid kmin*2^n+1 = q. */
    if (n < 32 && ((kmin<<n)>>n)==kmin && (kmin<<n)+1 <= q)
      r += q;

    sv->sieve_data.res[i] = r;
  }

#ifdef _WIN32
  InitializeCriticalSection(&sv->mutexA);
  InitializeCriticalSection(&sv->mutexB);
  sv->condC = CreateSemaphore(NULL,0,2147483647,NULL);
  sv->condC_waiting = 0;
#else
  pthread_mutex_init(&sv->mutexA,NULL);
  pthread_mutex_init(&sv->mutexB,NULL);
  pthread_cond_init(&sv->condC,NULL);
#endif

  if (blocks < 2)
    blocks = 2;
  else if (blocks > MAX_SIEVE_BLOCKS)
    blocks = MAX_SIEVE_BLOCKS;

  sv->cand_next = kmin;
  sv->free_blocks = blocks; /* begin with all blocks free */
  sv->curr_block = blocks;
  sv->sieve_done = 0;
  sv->all_done = 0;

  block_bytes = block_bytes/sizeof(unsigned long)*sizeof(unsigned long);
  chunk_bytes = chunk_bytes/sizeof(unsigned long)*sizeof(unsigned long);

  if (chunk_bytes < 2*sizeof(unsigned long))
    chunk_bytes = 2*sizeof(unsigned long);

  if (block_bytes/2 < chunk_bytes)
    block_bytes = 2*chunk_bytes;
  else if (block_bytes > 1<<28) /* 2^31 bits */
    block_bytes = 1<<28;

  if (block_bytes/2 < chunk_bytes)
    chunks = 2;
  else
    chunks = block_bytes/chunk_bytes;

  sv->chunk_size = chunk_bytes/sizeof(unsigned long);
  sv->block_size = chunk_bytes/sizeof(unsigned long)*chunks;
  sv->chunk_bits = sv->chunk_size*ULONG_BIT;
  sv->block_bits = sv->block_size*ULONG_BIT;
  sv->num_chunks = chunks;
  sv->num_blocks = blocks;

  for (i = 0; i < blocks; i++)
  {
    sv->block[i].free_chunks = chunks;
    sv->block[i].next_chunk = chunks;
    sv->block[i].base = 0;
    /* TODO: Use aligned malloc */
    sv->block[i].bits = xmalloc(sv->block_size * sizeof(unsigned long));
  }

#ifndef NDEBUG
  printf("create_gfn_sieve: Using %u odd primes, max=%u, lim=%u\n",
         num,sieve_primes[num-1],qmax);
  printf("chunk_size=%u (%ub), block_size=%u (%ub), chunks=%u, blocks=%u\n",
         sv->chunk_size, sv->chunk_size*sizeof(unsigned long),
         sv->block_size, sv->block_size*sizeof(unsigned long),
         sv->num_chunks, sv->num_blocks);
#endif

  return sv;
}

/* Free resources allocated by create_sieve.
 */
void destroy_sieve(sieve_t *sv)
{
  unsigned int i;

  assert(sv != NULL);

#ifdef _WIN32
  DeleteCriticalSection(&sv->mutexA);
  DeleteCriticalSection(&sv->mutexB);
  CloseHandle(sv->condC);
#else
  pthread_mutex_destroy(&sv->mutexA);
  pthread_mutex_destroy(&sv->mutexB);
  pthread_cond_destroy(&sv->condC);
#endif

  for (i = sv->num_blocks; i > 0; i--)
    free(sv->block[i-1].bits);

  free(sv);
}


/* Get the next candidate to be handed out by get_chunk. sv is unchanged.
 */
uint64_t next_chunk(sieve_t *sv)
{
  uint64_t cand_next;

  assert(sv != NULL);

  /* sv->cand_next must be read atomically, which requires mutual exclusion
     if the register size is less than 64 bits. */

  if (sizeof(void *) < sizeof(uint64_t))
  {
#ifdef _WIN32
    EnterCriticalSection(&sv->mutexA);
#else
    pthread_mutex_lock(&sv->mutexA);
#endif
  }

  cand_next = sv->cand_next;

  if (sizeof(void *) < sizeof(uint64_t))
  {
#ifdef _WIN32
    LeaveCriticalSection(&sv->mutexA);
#else
    pthread_mutex_unlock(&sv->mutexA);
#endif
  }

  return cand_next;
}

#if TRACE
uint64_t get_chunk(int th, sieve_t *sv, unsigned long **bitmap)
#else
uint64_t get_chunk(sieve_t *sv, unsigned long **bitmap)
#endif
{
  uint64_t base;
  unsigned int i, j;

  assert(sv != NULL);
  assert(bitmap != NULL);

#ifdef _WIN32
  EnterCriticalSection(&sv->mutexA);
#else
  pthread_mutex_lock(&sv->mutexA);
#endif

#if TRACE
  printf("Thread %d: entering get_chunk()\n",th);
#endif

 start:
  if (sv->all_done == 0)
  {
    if (sv->sieve_done || sv->free_blocks == 0 ||
#ifdef _WIN32
        !TryEnterCriticalSection(&sv->mutexB)
#else
        pthread_mutex_trylock(&sv->mutexB)
#endif
        )
    {
      /* Either no more sieving is required, there are no free blocks, or
         another thread is sieving already */

#if TRACE
      printf("Thread %d: not sieving, curr_block=%u\n",th,sv->curr_block);
#endif

      if (sv->curr_block == sv->num_blocks)
      {
#if TRACE
        printf("Thread %d: waiting for more chunks\n",th);
#endif

        /* Block until more chunks are ready */
#ifdef _WIN32
        sv->condC_waiting++;
        LeaveCriticalSection(&sv->mutexA);
        WaitForSingleObject(sv->condC,INFINITE);
        EnterCriticalSection(&sv->mutexA);
#else
        pthread_cond_wait(&sv->condC,&sv->mutexA);
#endif
        goto start;
      }
    }
    else /* Sieving is needed, there is a free block, and we have the sieve */
    {
      /* Find a free block */
      for (i = 0; i < sv->num_blocks; i++)
        if (sv->block[i].free_chunks == sv->num_chunks) /* free block */
          break;

      assert(i < sv->num_blocks);

#if TRACE
      printf("Thread %d: sieving, using free block %u\n",th,i);
#endif

#ifdef _WIN32
      LeaveCriticalSection(&sv->mutexA);
#else
      pthread_mutex_unlock(&sv->mutexA);
#endif

      j = sieve(&sv->sieve_data,&sv->block[i].base,sv->block[i].bits,sv->block_size);
      assert(j != 0);

#ifdef _WIN32
      EnterCriticalSection(&sv->mutexA);
#else
      pthread_mutex_lock(&sv->mutexA);
#endif

      if (sv->sieve_data.cand_max <= sv->block[i].base + sv->block_bits*2)
      {
        unsigned int last_chunk, k;

#if TRACE
        printf("Thread %d: setting sieve_done\n",th);
#endif

        /* We have reached the end of the sieve */
        sv->sieve_done = 1;

        /* Clear bits of final chunk corresponding to candidates >= cand_max */
        last_chunk = (j-1)/sv->chunk_size;
        for (k = (last_chunk+1)*sv->chunk_size; j < k; j++)
          sv->block[i].bits[j] = 0;

        sv->block[i].free_chunks = (sv->num_chunks - last_chunk - 1);
        sv->block[i].next_chunk = 0;
      }
      else
      {
        assert(j == sv->block_size);
        sv->block[i].free_chunks = 0;
        sv->block[i].next_chunk = 0;
      }

      sv->free_blocks--;

      if (sv->curr_block == sv->num_blocks) /* no other chunks ready */
      {
#if TRACE
        printf("Thread %d: setting current block to %u\n",th,i);
#endif
        sv->curr_block = i;
      }
#ifdef _WIN32
      LeaveCriticalSection(&sv->mutexB);
#else
      pthread_mutex_unlock(&sv->mutexB);
#endif

#if TRACE
      printf("Thread %d: broadcasting\n",th);
#endif

      /* Unblock all threads waiting for a chunk */
#ifdef _WIN32
      if (sv->condC_waiting > 0)
      {
        ReleaseSemaphore(sv->condC,sv->condC_waiting,NULL);
        sv->condC_waiting = 0;
      }
#else
      pthread_cond_broadcast(&sv->condC);
#endif
    }

    /* Now mutexA is locked and there is at least one ready chunk. */
    i = sv->curr_block;
    assert(i < sv->num_blocks);
    assert(sv->block[i].free_chunks < sv->num_chunks);
    assert(sv->block[i].next_chunk < sv->num_chunks);
    assert(sv->cand_next ==
           sv->block[i].base+(uint64_t)sv->chunk_bits*sv->block[i].next_chunk*2);

    /* Return this chunk */
    base = sv->cand_next;
    *bitmap = &sv->block[i].bits[sv->block[i].next_chunk*sv->chunk_size];

#if TRACE
    printf("Thread %d: allocating chunk %u (%"PRIu64") from block %u\n",
           th,sv->block[i].next_chunk,base,i);
#endif

    /* Advance to next chunk */
    sv->block[i].next_chunk++;
    sv->cand_next += sv->chunk_bits*2;

    if (sv->cand_next >= sv->sieve_data.cand_max)
    {
      /* Last chunk for this sieve */
      assert(sv->sieve_done == 1);
#if TRACE
      printf("Thread %d: setting all_done\n",th);
#endif
      sv->all_done = 1;
    }
    else if (sv->block[i].next_chunk == sv->num_chunks)
    {
      /* Last chunk for this block. Locate next block */
      if ((j = i+1) == sv->num_blocks)
        j = 0;
      while (j != i)
      {
        if (sv->block[j].free_chunks < sv->num_chunks)
          if (sv->block[j].base == sv->cand_next)
            break;
        if (++j == sv->num_blocks)
          j = 0;
      }

      if (j == i)
      {
#if TRACE
        printf("Thread %d: no blocks with ready chunks\n",th);
#endif
        /* There are no blocks with ready chunks */
        sv->curr_block = sv->num_blocks;
      }
      else
      {
#if TRACE
        printf("Thread %d: next block with ready chunks is %u\n",th,j);
#endif

        assert(sv->block[j].next_chunk == 0);
        sv->curr_block = j;
      }
    }
  }
  else /* sv->all_done */
  {
#if TRACE
    printf("Thread %d: all_done set\n",th);
#endif
    base = sv->sieve_data.cand_max;
  }

#if TRACE
  printf("Thread %d: leaving get_chunk()\n",th);
#endif

#ifdef _WIN32
  LeaveCriticalSection(&sv->mutexA);
#else
  pthread_mutex_unlock(&sv->mutexA);
#endif

  return base;
}

#if TRACE
void free_chunk(int th,sieve_t *sv, uint64_t chunk)
#else
void free_chunk(sieve_t *sv, uint64_t chunk)
#endif
{
  unsigned int i;

  assert(sv != NULL);
  assert(chunk != 0);

#ifdef _WIN32
  EnterCriticalSection(&sv->mutexA);
#else
  pthread_mutex_lock(&sv->mutexA);
#endif

#if TRACE
  printf("Thread %d: entering free_chunk()\n",th);
#endif

  assert(sv->free_blocks < sv->num_blocks);

  /* Find which block this chunk belongs to and increment its free chunk
     count. */
  for (i = 0; i < sv->num_blocks; i++)
    if (sv->block[i].free_chunks < sv->num_chunks) /* not a free block */
      if (chunk >= sv->block[i].base)
        if (chunk < sv->block[i].base + (uint64_t)sv->block_bits*2)
        {
          assert((chunk - sv->block[i].base)%(sv->chunk_bits*2) == 0);
#if TRACE
          printf("Thread %d: freeing chunk %u (%"PRIu64") in block %u\n",th,(unsigned int)((chunk - sv->block[i].base)/(sv->chunk_bits*2)),chunk,i);
#endif
          if (++sv->block[i].free_chunks == sv->num_chunks)
          {
            /* This was the last non-free chunk so the block is now free */
            assert(sv->free_blocks < sv->num_blocks);
            if (0 == sv->free_blocks++)
            {
              /* It is possible, if there are fewer blocks than threads,
               that threads could be waiting for a free block. */
#ifdef _WIN32
              if (sv->condC_waiting > 0)
              {
                ReleaseSemaphore(sv->condC,sv->condC_waiting,NULL);
                sv->condC_waiting = 0;
              }
#else
              pthread_cond_broadcast(&sv->condC);
#endif
            }
#if TRACE
            printf("Thread %d: block %u freed\n",th,i);
#endif
          }

          break;
        }

  /* If no block contains this chunk then it was an invalid argument */
  assert(i < sv->num_blocks);

#if TRACE
  printf("Thread %d: leaving free_chunk()\n",th);
#endif

#ifdef _WIN32
  LeaveCriticalSection(&sv->mutexA);
#else
  pthread_mutex_unlock(&sv->mutexA);
#endif
}
