/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.cu -- (C) Ken Brazier August 2010.

   Proth Prime Search sieve OpenCL Kernel portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

//#define KERNEL_ONLY
//#include "appcl.h"
//#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable


/*** Kernel Helpers ***/
// Special thanks to Alex Kruppa for introducing me to Montgomery REDC math!
/* Compute a^{-1} (mod 2^(32 or 64)), according to machine's word size */

ulong
invmod2pow_ul (const ulong n)
{
  ulong r;
  //uint ir;
  const uint in = (uint)n;

  //ASSERT (n % 2UL != 0UL);

  // Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
  // correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
  // Newton iterations are enough.
  r = (n+n+n) ^ ((ulong)2);
  // Newton iteration
  r += r - (uint) r * (uint) r * in;
  r += r - (uint) r * (uint) r * in;
  r += r - (uint) r * (uint) r * in;
  r += r - r * r * n;

  return r;
}

ulong mulmod_REDC (const ulong a, const ulong b, 
    const ulong N, const ulong Ns)
{
  ulong rax, rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax = a*b;
  // No mad_hi here, because I want this independent of the rax computation.
  rcx = mul_hi(a,b);	// rcx = Th
  
  rax *= Ns;		// rax = (T*Ns) mod 2^64 = m
  rcx += (rax!=0)?1:0;	// if rax != 0, increase rcx
  rax = mad_hi(rax, N, rcx);
  // compute (rdx + rcx) mod N
  //rax += rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  /*
#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, b, N))
  {
    fprintf (stderr, "%sError, mulredc(%lu,%lu,%lu) = %lu\n", bmprefix(), a, b, N, rax);
    bexit(1);
  }
#endif
*/

  return rax;
}

// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
ulong onemod_REDC(const ulong N, ulong rax) {
  ulong rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //rcx = 0;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx = (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = mad_hi(rax, N, rcx);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  return rax;
}

// Like mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
ulong mod_REDC(const ulong a, const ulong N, const ulong Ns) {
//#ifndef DEBUG64
  return onemod_REDC(N, Ns*a);
  /*
//#else
  const ulong r = onemod_REDC(N, Ns*a);

  if (longmod (r, 0, N) != mulmod(a, 1, N)) {
    fprintf (stderr, "%sError, redc(%lu,%lu) = %lu\n", bmprefix(), a, N, r);
    bexit(1);
  }

  return r;
//#endif
*/
}

// Compute T=a<<s; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
// rax is passed in as a * Ns.
// rax's original value is destroyed, just to keep the register count down.
ulong shiftmod_REDC (const ulong a, 
    const ulong N, ulong rax)
{
  ulong rcx;

  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax <<= D_MONT_NSTEP; // So this is a*Ns*(1<<s) == (a<<s)*Ns.
  rcx = a >> D_NSTEP;
  //"movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
  //"imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
  //rax *= Ns;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx += (rax!=0)?1:0;
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = mad_hi(rax, N, rcx);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  /*
#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, ((ulong)1)<<D_MONT_NSTEP, N))
  {
    fprintf (stderr, "%sError, shiftredc(%lu,%u,%lu) = %lu\n", bmprefix(), a, D_MONT_NSTEP, N, rax);
    bexit(1);
  }
#endif
*/

  return rax;
}

// A Left-to-Right version of the powmod.  Calcualtes 2^-(first 6 bits), then just keeps squaring and dividing by 2 when needed.
ulong
invpowmod_REDClr (const ulong N, const ulong Ns, int bbits, ulong r) {
  // Now work through the other bits of nmin.
  for(; bbits >= 0; --bbits) {
    // Just keep squaring r.
    r = mulmod_REDC(r, r, N, Ns);
    // If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod N).
    if(D_NMIN & (1u << bbits)) {
      r += (r&1)?N:0;
      r >>= 1;
    }
  }

#ifdef DEBUG64
  //assert (mod_REDC (r, N, Ns) == invmod(powmod (D_NMIN, N), N));
#endif

  // Convert back to standard.
  //r = mod_REDC (r, N, Ns);

  return r;
}


// *** KERNELS ***

// Start checking N's.
__kernel void start_ns(__global ulong * P, __global ulong * Ps, __global ulong * K, __global uchar * factor_found_arr,
                                 // Device constants
                                 const ulong d_r0			// 4
                                 ) {
  uint n = D_NMIN; // = nmin;
  uint i;
  ulong k0;
  //uchar my_factor_found = 0;
  ulong my_P, my_Ps;
  i = get_global_id(0);
  my_P = P[i];

  // Better get this done before the first mulmod.
  my_Ps = -invmod2pow_ul (my_P); /* Ns = -N^{-1} % 2^64 */

  // Calculate k0, in Montgomery form.
  k0 = invpowmod_REDClr(my_P, my_Ps, D_BBITS, d_r0);

  //if(my_P == 42070000070587) printf("%lu^-1 = %lu (GPU)\n", my_P, my_Ps);

#ifdef D_SEARCH_PROTH
  k0 = my_P-k0;
#endif

  //my_factor_found = 0;
  //d_check_some_ns(my_P, my_Ps, k0, n, my_factor_found, i);

  //i = get_global_id(0);
  factor_found_arr[i] = 0;
  Ps[i] = my_Ps;
  K[i] = k0;
}

// Continue checking N's.
__kernel void check_more_ns(__global ulong * P, __global ulong * Psarr, __global ulong * K, __global uchar * factor_found_arr, const uint N
                                 // Device constants
#ifndef D_KMIN
                                 , const ulong d_kmin		// 5
#endif
#ifndef D_KMAX
                                 , const ulong d_kmax		// 6
#endif
                                 ) {
  uint i = get_global_id(0);
  const ulong my_P = P[i];
  const ulong Ps = Psarr[i];
  uint n = N;
  ulong k0 = K[i];
  uchar my_factor_found = factor_found_arr[i];
  ulong kpos, kPs;
  uint l_nmax = n + D_KERNEL_NSTEP;
  if(l_nmax > D_NMAX) l_nmax = D_NMAX;

#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (uint)k0, D_KERNEL_NSTEP);
#endif
  do { // Remaining steps are all of equal size nstep
    // Get K from the Montgomery form.
    // This is equivalent to mod_REDC(k, my_P, Ps), but the intermediate kPs value is kept for later.
    kPs = k0 * Ps;
    kpos = onemod_REDC(my_P, kPs);
    //i = __ffsll(kpos)-1;
    i = (uint)kpos;
    if(i != 0) {
      //i=(__float_as_int(__uint2float_rz(i & -i))>>23)-0x7f;
      i=31 - clz (i & -i);
    } else {
      i = (uint)(kpos>>32);
      i=63 - clz (i & -i);
    }

    kpos >>= i;
#ifdef D_KMAX
    if (kpos <= D_KMAX) {
#else
    if (kpos <= d_kmax) {
#endif
//#ifdef _DEVICEEMU
      //printf("%lu | %lu*2^%u+1 (P[%d])\n", my_P, kpos, n+i, get_global_id(0));
//#endif
      // Just flag this if kpos <= d_kmax.
#ifdef D_KMIN
      if(kpos >= D_KMIN)
#else
      if(kpos >= d_kmin)
#endif
        my_factor_found = 1;
    }

    // Proceed to the K for the next N.
    // kPs is destroyed, just to keep the register count down.
    k0 = shiftmod_REDC(k0, my_P, kPs);
    n += D_NSTEP;
  } while (n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (uint)k0);
#endif
  i = get_global_id(0);
  factor_found_arr[i] = my_factor_found;
  if(n < D_NMAX) {
    K[i] = k0;
  }
}
