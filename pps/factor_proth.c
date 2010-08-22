/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* factor_proth.c -- (C) Ken Brazier August 2010.
   Factor a Proth number with small primes, and see if it breaks.
   To be used to test whether potential larger factors are useful.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include "inttypes.h"
#include "factor_proth.h"
#ifdef DEBUG64
#include "main.h"
#endif
#define MAX_SIEVE (32780/2)
// Global array (local to this file) of small primes:
static int small_primes[3514];

// Returns a list of small primes, 3-65537 inclusive.
// Doesn't have to be extremely fast; only done once.
void sieve_small_primes(int min) {
  char sieve[MAX_SIEVE];
  int i, j;

  // Initialize sieve to all prime, but 1/2==0.
  sieve[0] = 0;
  for(i=1; i < MAX_SIEVE; i++) sieve[i] = 1;

  // 181 is the last prime factoring in this range.
  for(j=3; j < 182; j += 2)
    if(sieve[j/2])
      for(i=(j*j)/2; i < MAX_SIEVE; i+=j)
        sieve[i] = 0;

  for(i=min/2,j=0; i < MAX_SIEVE; i++)
    if(sieve[i]) small_primes[j++] = i+i+1;
}

// find the log base 2 of a number.  Could use to be fast.
static int lg2(unsigned int v) {
#ifdef __GNUC__
  return 31-__builtin_clz(v);
#else
  /*
     int r = 0; // r will be lg(v)

     while (v >>= 1) // unroll for more speed...
     {
     r++;
     }
     */
  register unsigned int r; // result of log2(v) will go here
  register unsigned int shift;

  r =     (v > 0xFFFF) << 4; v >>= r;
  shift = (v > 0xFF  ) << 3; v >>= shift; r |= shift;
  shift = (v > 0xF   ) << 2; v >>= shift; r |= shift;
  shift = (v > 0x3   ) << 1; v >>= shift; r |= shift;
  r |= (v >> 1);

  return r;
#endif
}

// REDC stuff, for 16/32 bits, instead of 64/128 bits.
// Test method first:
#ifdef DEBUG64
static unsigned int inv2powmod(unsigned int N, unsigned int P) {
  unsigned int r = (P+1)/2;
  int bbits;

  bbits=0;
  while(N >> bbits) bbits++;
  bbits -= 2;

  // Now work through the other bits of N.
  for(; bbits >= 0; --bbits) {
    // Just keep squaring r.
    r = (r*r)%P;
    // If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod P).
    if(N & (1u << bbits)) {
      r += (r&1)?P:0;
      r >>= 1;
    }
  }
  return r;
}
#endif
static unsigned int invmod2pow_ul (const unsigned int n)
{
  unsigned int r;
  const unsigned short int in = (unsigned short int)n;

  //ASSERT (n % 2UL != 0UL);

  // Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
  // correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
  // Newton iterations are enough.
  r = (n+n+n) ^ ((unsigned int)2);
  // Newton iteration
  r += r - (unsigned short int) r * (unsigned short int) r * in;
  r += r - (unsigned short int) r * (unsigned short int) r * in;
  //r += r - (unsigned int) r * (unsigned int) r * in;
  r += r - r * r * n;

  return r;
}

static unsigned int mulmod_REDC (const unsigned int a, const unsigned int b, 
    const unsigned short int N, const unsigned short int Ns)
{
  unsigned int r;
  unsigned short int rax, rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  r = ((int)a)*((int)b);
  rcx = r >> 16;
  rax = (unsigned short int)r;
  //"movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
  //rcx = rdx;
  //"imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
  rax *= Ns;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx += (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = (unsigned short int)(((int)rax*(int)N)>>16);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  r = (unsigned int)rax + (unsigned int)rcx;
  rcx = r - N;
  rax = (r>(unsigned int)N)?rcx:r;

#ifdef DEBUG64
  if ((unsigned int)(rax<<16)%(unsigned int)N != (a*b)%(unsigned int)N)
  {
    fprintf (stderr, "Error, mulredc(%u,%u,%u) = %u\n", a, b, N, rax);
    bexit(1);
  }
#endif

  return rax;
}

/*
// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
static unsigned int onemod_REDC(const unsigned short int N, unsigned int rax) {
  unsigned short int rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //rcx = 0;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx = (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = (((int)rax*(int)N)>>16) + rcx;
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rcx = rax - N;
#ifdef DEBUG64
  if(rax > N && rcx > N) fprintf(stderr, "%u > %u and so is %u\n", rax, N, rcx);
#endif
  rax = (rax>N)?rcx:rax;


  return rax;
}

// Like mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
static unsigned int mod_REDC(const unsigned short int a, const unsigned short int N, const unsigned short int Ns) {
#ifndef DEBUG64
  return onemod_REDC(N, Ns*a);
#else
  const unsigned int r = onemod_REDC(N, ((int)Ns)*((int)a));

  if ((unsigned int)(r<<16)%(unsigned int)N != (a)%(unsigned int)N) {
    fprintf (stderr, "Error, redc(%lu,%lu) = %lu\n", a, N, r);
    bexit(1);
  }

  return r;
#endif
}
*/
// A Left-to-Right version of the powmod.  Calcualtes 2^-(first 4 bits), then just keeps squaring and dividing by 2 when needed.
static unsigned int invpowmod_REDClr (const unsigned int N, const unsigned short int P, const unsigned short int Ps, int bbits, unsigned short int r) {
  // First, make sure r < P.
  if(r >= P) r %= P;
  //if(r > P) fprintf(stderr, "original R %u > %u\n", r, P);

  // Now work through the other bits of N.
  for(; bbits >= 0; --bbits) {
    // Just keep squaring r.
    r = mulmod_REDC(r, r, P, Ps);
    //if(r > P) fprintf(stderr, "subsequent R %u > %u\n", r, P);
    // If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod P).
    if(N & (1u << bbits)) {
      r += (r&1)?P:0;
      r >>= 1;
    }
  }

  // Convert back to standard.
  r = mulmod_REDC (r, 1, P, Ps);
#ifdef DEBUG64
  if(r != inv2powmod(N, P)) {
    fprintf (stderr, "Error, inv2powmod(%u,%u) == %u, not %u\n", N, P, inv2powmod(N, P), r);
    bexit(1);
  }
#endif

  return r;
}

// Input: K>1, N>=16, sign (-1 or 1), p < 2^15 (32768)
static int try_factor(uint64_t K, unsigned int N, int sign, unsigned int p) {
  unsigned int kcalc;
  int bbits;
  unsigned int r0;

  // Prepare constants:
  bbits = lg2(N);
  //assert(r0 <= 32);
  /*
     if(bbits < 4) {
     fprintf(stderr, "Error: N too small at %d (must be at least 16).\n", N);
     bexit(1);
     }
     */
  // r = 2^-i * 2^16 (mod P)
  // If i is large there's a chance no mod is needed!
  r0 = ((unsigned int)1) << (16-(N >> (bbits-3)));

  bbits = bbits-4;

  kcalc = invpowmod_REDClr(N, p, -invmod2pow_ul(p), bbits, r0);
  if(sign > 0) kcalc = p-kcalc;
  return(kcalc == (unsigned int)(K%(uint64_t)p));
}

// Try factoring K and N with each prime. (In pairs, or quadruples, or even octuples with SSE2?)
int try_all_factors(uint64_t K, unsigned int N, int sign) {
  int *p;
  for(p=small_primes; *p < 32768; p++)
    if(try_factor(K, N, sign, *p))
      return *p;
  return 0;
}

#ifndef NDEBUG
int main(int argc, char **argv) {
  int f;
  uint64_t K, temp;
  unsigned int N;
  int sign = 1;

  // Sieve for small primes, < 2^15
  sieve_small_primes(0);

  // Temporary: Fiddle with the array to match what I'm doing with shifts and mods.
  //small_primes[8] = 31;
  //small_primes[9] = 40960;

  // Get K and N.
  //K = 7339; N=1173500; sign=1;
  //K = 355; N=54; sign=1;

  //for(K=3; K < 20000; K+=2) {
  while(fscanf(stdin, "%lu | %lu*2^%u+1\n", &temp, &K, &N) == 3) {
    f = try_all_factors(K, N, sign);
    if(f == 0) printf("%lu | %lu*2^%u%c1\n", temp, K, N, (sign<0)?'-':'+');
    else printf("sm %d | %lu*2^%u%c1\n", f, K, N, (sign<0)?'-':'+');
  }

  return 0;
}
#endif
