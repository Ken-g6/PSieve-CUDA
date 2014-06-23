/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.cu -- (C) Ken Brazier August 2010.

   Proth Prime Search sieve OpenCL Kernel portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

//#define KERNEL_ONLY
#ifdef _DEVICEEMU
#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif
//#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
/*
#define D_MONT_NSTEP 35u
#define D_NSTEP 29u
#define D_NMIN 100u
#define D_NMAX 2000000u
#define D_BBITS 1u
#define D_KERNEL_NSTEP 4096u
#define D_KMIN 1200u
#define D_KMAX 9999u
#define VECSIZE 2
*/
#if VECSIZE == 1
#define VLONG ulong
#define VSIGNED_LONG long
#define VINT uint
#define VECTORIZE(x) (x)
#define V2VINT(x) ((uint)(x))
#define V2VLONG(x) ((ulong)(x))
#define VDOTX v
#else
#define VDOTX v.x
#define VECSIZEIT2(x,y) x##y
#define VECSIZEIT(x,y) VECSIZEIT2(x,y)
#define VLONG VECSIZEIT(ulong, VECSIZE)
#define VSIGNED_LONG VECSIZEIT(long, VECSIZE)
#define VINT VECSIZEIT(uint, VECSIZE)
#if VECSIZE == 2
#define VECTORIZE(x) ((x),(x))
#elif VECSIZE == 3
#define VECTORIZE(x) ((x),(x),(x))
#elif VECSIZE == 4
#define VECTORIZE(x) ((x),(x),(x),(x))
#endif
#define V2VINT(x) VECSIZEIT(convert_uint, VECSIZE)(x)
#define V2VLONG(x) VECSIZEIT(convert_ulong, VECSIZE)(x)
#endif


/*** Kernel Helpers ***/
// Special thanks to Alex Kruppa for introducing me to Montgomery REDC math!
/* Compute a^{-1} (mod 2^(32 or 64)), according to machine's word size */

VLONG
invmod2pow_ul (const VLONG n)
{
  VLONG r;
  //VINT ir;
  const VINT in = V2VINT(n);

  //ASSERT (n % 2UL != 0UL);

  // Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
  // correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
  // Newton iterations are enough.
  r = (n+n+n) ^ ((VLONG)2);
  // Newton iteration
  r += r - V2VLONG(V2VINT(r) * V2VINT(r) * in);
  r += r - V2VLONG(V2VINT(r) * V2VINT(r) * in);
  r += r - V2VLONG(V2VINT(r) * V2VINT(r) * in);
  r += r - r * r * n;

  return r;
}

VLONG mulmod_REDC (const VLONG a, const VLONG b, 
    const VLONG N, const VLONG Ns)
{
  VLONG rax, rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax = a*b;
  // No mad_hi here, because I want this independent of the rax computation.
  rcx = mul_hi(a,b);	// rcx = Th
  
  rax *= Ns;		// rax = (T*Ns) mod 2^64 = m
  rcx += (rax!=(VLONG)VECTORIZE(0))?((VLONG)VECTORIZE(1)):((VLONG)VECTORIZE(0));	// if rax != 0, increase rcx
  rax = mad_hi(rax, N, rcx);
  // compute (rdx + rcx) mod N
  //rax += rcx;
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  return rax;
}

// mulmod_REDC(1, 1, N, Ns)
// But note that mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
VLONG onemod_REDC(const VLONG N, VLONG rax) {
  VLONG rcx;

  // Akruppa's way, Compute T=a*b; m = (T*Ns)%2^64; T += m*N; if (T>N) T-= N;
  //rcx = 0;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  //rcx = (rax!=(VLONG)VECTORIZE(0))?((VLONG)VECTORIZE(1)):((VLONG)VECTORIZE(0));	// if rax != 0, increase rcx
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = mad_hi(rax, N, ((VLONG)VECTORIZE(1)));
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  return rax;
}

// Like mulmod_REDC(a, 1, N, Ns) == mulmod_REDC(1, 1, N, Ns*a).
VLONG mod_REDC(const VLONG a, const VLONG N, const VLONG Ns) {
//#ifndef DEBUG64
  return onemod_REDC(N, Ns*a);
  /*
//#else
  const VLONG r = onemod_REDC(N, Ns*a);

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
VLONG shiftmod_REDC (const VLONG a, 
    const VLONG N, VLONG rax)
{
  VLONG rcx;

  //( "mulq %[b]\n\t"           // rdx:rax = T 			Cycles 1-7
  rax = rax << D_MONT_NSTEP; // So this is a*Ns*(1<<s) == (a<<s)*Ns.
  rcx = a >> D_NSTEP;
  //"movq %%rdx,%%rcx\n\t"	// rcx = Th			Cycle  8
  //"imulq %[Ns], %%rax\n\t"  // rax = (T*Ns) mod 2^64 = m 	Cycles 8-12 
  //rax *= Ns;
  //"cmpq $1,%%rax \n\t"      // if rax != 0, increase rcx 	Cycle 13
  //"sbbq $-1,%%rcx\n\t"	//				Cycle 14-15
  rcx += (rax!=(VLONG)VECTORIZE(0))?((VLONG)VECTORIZE(1)):((VLONG)VECTORIZE(0));	// if rax != 0, increase rcx
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = mad_hi(rax, N, rcx);
  //"lea (%%rcx,%%rdx,1), %[r]\n\t" // compute (rdx + rcx) mod N  C 20 
  rcx = rax - N;
  rax = (rax>N)?rcx:rax;

  /*
#ifdef DEBUG64
  if (longmod (rax, 0, N) != mulmod(a, ((VLONG)1)<<D_MONT_NSTEP, N))
  {
    fprintf (stderr, "%sError, shiftredc(%lu,%u,%lu) = %lu\n", bmprefix(), a, D_MONT_NSTEP, N, rax);
    bexit(1);
  }
#endif
*/

  return rax;
}

// Same function, for a constant NSTEP.
// Multiply two 32-bit integers to get a 64-bit result.
VLONG mul_wide_u32 (const VINT a, const VINT b) {
  return upsample(mul_hi(a, b), a*b);
}

#if D_NSTEP == 32
// Same function, for a constant NSTEP of 32.
VLONG shiftmod_REDC32 (VLONG rcx, const VLONG N, const VINT rax)
{
  //unsigned int temp;
  //uint64_t rcx;

  rcx >>= 32;
  //temp = ((unsigned int)(N>>32));
  // This isn't in an asm, but should be compiled as mad.hi.u32. Four cycles.
  // We know this can fit an unsigned int because (N-1)*(N-1) = N^2-2N-1, so adding the equivalent of 1N is OK.
  rcx += V2VLONG(mad_hi(rax, V2VINT(N), (rax!=(VINT)VECTORIZE(0))?((VINT)VECTORIZE(1)):((VINT)VECTORIZE(0))));
  // A wide multiply should take ~8 cycles here, and the add can't be combined.  But that's better than 16.
  rcx += mul_wide_u32(rax, V2VINT(N>>32));
  // Two cycles for this one add!
  //rcx = ((((uint64_t)__umulhi((unsigned int)rax,(unsigned int)(N>>32))) << 32) | (((unsigned int)rax)*((unsigned int)(N>>32)))) + rcx;
  //rax = ((uint64_t)((unsigned int)rax))*((uint64_t)((unsigned int)(N>>32))) + temp;

  // And the rest is normal, but squashed.
  rcx = (rcx>N)?(rcx-N):rcx;
  return rcx;
}
#endif
#if D_NSTEP < 32
// Same function for nstep < 32. (SMall.)
// Third argument must be passed in as only the low register, as we're effectively left-shifting 32 plus a small number.
VLONG shiftmod_REDCsm (VLONG rcx, const VLONG N, VINT rax)
{
  //unsigned int temp;
  //uint64_t rcx;

  rax <<= (D_MONT_NSTEP-32);
  rcx >>= D_NSTEP;
  //temp = ((unsigned int)(N>>32));
  // This isn't in an asm, but should be compiled as mad.hi.u32. Four cycles.
  // We know this can fit an unsigned int because (N-1)*(N-1) = N^2-2N-1, so adding the equivalent of 1N is OK.
  rcx += V2VLONG(mad_hi(rax, V2VINT(N), (rax!=(VINT)VECTORIZE(0))?((VINT)VECTORIZE(1)):((VINT)VECTORIZE(0))));
  // A wide multiply should take ~8 cycles here, and the add can't be combined.  But that's better than 16.
  rcx += mul_wide_u32(rax, V2VINT(N>>32));
  // Two cycles for this one add!
  //rcx = ((((uint64_t)__umulhi((unsigned int)rax,(unsigned int)(N>>32))) << 32) | (((unsigned int)rax)*((unsigned int)(N>>32)))) + rcx;
  //rax = ((uint64_t)((unsigned int)rax))*((uint64_t)((unsigned int)(N>>32))) + temp;

  // And the rest is normal, but squashed.
  rcx = (rcx>N)?(rcx-N):rcx;
  return rcx;
}
#endif

// Same function, for a custom NSTEP.  step should be no bigger than NSTEP.
VLONG shiftmod_REDCX (const VLONG a, const VLONG N, VLONG rax, const uint step) { \
  VLONG rcx; \
  rax = rax << (64-step); \
  rcx = a >> step; \
  rcx += (rax!=(VLONG)VECTORIZE(0))?((VLONG)VECTORIZE(1)):((VLONG)VECTORIZE(0)); \
  rax = mad_hi(rax, N, rcx); \
  rcx = rax - N; \
  rax = (rax>N)?rcx:rax; \
  return rax; \
}

// A Left-to-Right version of the powmod.  Calcualtes 2^-(first 6 bits), then just keeps squaring and dividing by 2 when needed.
VLONG
invpowmod_REDClr (const VLONG N, const VLONG Ns, int bbits, const ulong r0) {
  VLONG r = r0;
  // Now work through the other bits of nmin.
  for(; bbits >= 0; --bbits) {
    // Just keep squaring r.
    r = mulmod_REDC(r, r, N, Ns);
    // If there's a one bit here, multiply r by 2^-1 (aka divide it by 2 mod N).
    if(D_NMIN & (1u << bbits)) {
      r += ((r&((VLONG)VECTORIZE(1))) == (VLONG)VECTORIZE(1))?N:((VLONG)VECTORIZE(0));
      r = r >> 1u;
    }
  }

#ifdef DEBUG64
  //assert (mod_REDC (r, N, Ns) == invmod(powmod (D_NMIN, N), N));
#endif

  // Convert back to standard.
  r = mod_REDC (r, N, Ns);

  return r;
}


// *** KERNELS ***
#if VECSIZE == 1
 #define VLOAD(X, Y) Y[X]
 #define VSTORE(X, Y, Z) Z[Y] = (X)
#else
 #define VLOAD VECSIZEIT(vload,VECSIZE)
 #define VSTORE VECSIZEIT(vstore,VECSIZE)
#endif
// Set up to check N's by getting in position with division only.  Doesn't use the CPU for a powmod, but should only be used for small nmin.
__kernel void start_ns(__global ulong * P, __global ulong * Ps, __global ulong * K, __global uint * factor_found_arr,
                                 // Device constants
                                 const ulong d_r0			// 4
                                 ) {
  uint n = D_NMIN; // = nmin;
  uint i;
  VLONG k0;
  VLONG my_P, my_Ps;
  i = get_global_id(0);
  my_P = VLOAD(i, P); //P[i];

  my_Ps = -invmod2pow_ul (my_P); // Ns = -N^{-1} % 2^64
  // Might need to load this from the CPU, in case the bug is there.
  //my_Ps = VLOAD(i, Ps)/* Ns = -N^{-1} % 2^64 */

  // Calculate k0, in Montgomery form.
  //k0 = invpowmod_REDClr(my_P, my_Ps, D_BBITS, d_r0);
  while(n >= 128) {
    k0 = mod_REDC(k0, my_P, my_Ps);
    n -= 64;
  }
  if(n > 64) {
    k0 = shiftmod_REDCX(k0, my_P, k0*my_Ps, n/2);
    n -= n/2;
  }
  if(n > 0) {
    k0 = shiftmod_REDCX(k0, my_P, k0*my_Ps, n);
  }

  //if(my_P == 42070000070587) printf("%lu^-1 = %lu (GPU)\n", my_P, my_Ps);

#ifdef D_SEARCH_PROTH
  k0 = my_P-k0;
#endif

  //my_factor_found = 0;
  //d_check_some_ns(my_P, my_Ps, k0, n, my_factor_found, i);

  //i = get_global_id(0);
  //factor_found_arr[i] = 0;
  VSTORE((VINT)0, i, factor_found_arr);
  //Ps[i] = my_Ps;
  VSTORE(my_Ps, i, Ps);
  //K[i] = k0;
  VSTORE(k0, i, K);
}
/*
// Start checking N's.
__kernel void start_ns(__global ulong * P, __global ulong * Ps, __global ulong * K, __global uint * factor_found_arr,
                                 // Device constants
                                 const ulong d_r0			// 4
                                 ) {
  uint n = D_NMIN; // = nmin;
  uint i;
  VLONG k0;
  VLONG my_P, my_Ps;
  i = get_global_id(0);
  my_P = VLOAD(i, P); //P[i];

  // Better get this done before the first mulmod.
  my_Ps = -invmod2pow_ul (my_P); // Ns = -N^{-1} % 2^64

  // Calculate k0, in Montgomery form.
  k0 = invpowmod_REDClr(my_P, my_Ps, D_BBITS, d_r0);

  //if(my_P == 42070000070587) printf("%lu^-1 = %lu (GPU)\n", my_P, my_Ps);

#ifdef D_SEARCH_PROTH
  k0 = my_P-k0;
#endif

  //my_factor_found = 0;
  //d_check_some_ns(my_P, my_Ps, k0, n, my_factor_found, i);

  //i = get_global_id(0);
  //factor_found_arr[i] = 0;
  VSTORE(VINT(0), i, factor_found_arr);
  //Ps[i] = my_Ps;
  VSTORE(my_Ps, i, Ps);
  //K[i] = k0;
  VSTORE(k0, i, K);
}
 */

// Continue checking N's.
//i=(__float_as_int(__VINT2float_rz(i & -i))>>23)-0x7f;
#define VEC_CTZLL(I,X) \
      if(I.X != 0) { \
        I.X=31u - clz (I.X & -I.X); \
      } else { \
        I.X = (uint)(kpos.X>>32); \
        I.X=63u - clz (I.X & -I.X); \
      }

#ifdef SEARCH_TWIN
#define GPUNAME "TPS GPU"
#else
#define GPUNAME "PPS GPU"
#endif
#ifdef _DEVICEEMU
#define DEBUG_PRINT_RESULT2(X) printf("%lu | %lu*2^%u(+/-)1 (%s)\n", my_P.X, (kpos.X >> v.X), n+v.X, GPUNAME);
#define DEBUG_PRINT_RESULT(X) DEBUG_PRINT_RESULT2(X)
//#define DEBUG_PRINT_RESULT(X)
#else
#define DEBUG_PRINT_RESULT(X)
#endif
#ifdef SEARCH_TWIN
#define NSTEP_COMP <=
#else
#define NSTEP_COMP <
#endif
    // Just flag this if kpos <= d_kmax.
#ifdef D_KMAX
#ifdef D_KMIN
// This test is still useful if it may reduce reported candidates by, say, half.
#define MFFTEST if((kpos.X >> v.X) >= D_KMIN)
// Removing these extra tests may result in a few more reported candidates, but should speed up the testing slightly.
//&& v.X NSTEP_COMP D_NSTEP && n+v.X NSTEP_COMP l_nmax
#else
#define MFFTEST
#endif
#define VEC_FLAG_TEST(X) \
    if ((kpos.X >> v.X) <= D_KMAX) { \
      DEBUG_PRINT_RESULT(X) \
      MFFTEST my_factor_found.X |= 1; \
    }
#define VEC_FAST_FLAG_TEST(X) \
  if ((((uint)(kpos.X >> 32))>>v.X) == 0) { \
    if(((uint)(kpos.X >> v.X)) <= D_KMAX) { \
      DEBUG_PRINT_RESULT(X) \
      MFFTEST my_factor_found.X |= 1; \
    } \
  }
#else
#ifdef D_KMIN
#define THE_KMIN D_KMIN
#else
#define THE_KMIN d_kmin
#endif
#define VEC_FLAG_TEST(X) \
    if ((kpos.X >> v.X) <= d_kmax) { \
      DEBUG_PRINT_RESULT(X) \
      if((kpos.X >> v.X) >= THE_KMIN) \
        my_factor_found.X |= 1; \
    }
#define VEC_FAST_FLAG_TEST(X) VEC_FLAG_TEST(X)
#endif
#if(VECSIZE == 1)
// At VECSIZE 1, if this is called, we already know v == 0.
#define ALL_CTZLL \
        v = (uint)(kpos>>32); \
        v=63u - clz (v & -v);

/*
      if(v != 0) { \
        v=31u - clz (v & -v); \
      } else { \
      }
      */
#ifdef D_KMAX
#ifdef D_KMIN
// This test is still useful if it may reduce reported candidates by, say, half.
#define MFFONETEST if((kpos >> v) >= D_KMIN)
// Removing these extra tests may result in a few more reported candidates, but should speed up the testing slightly.
// && v NSTEP_COMP D_NSTEP && n+v NSTEP_COMP l_nmax)
#else
#define MFFONETEST
#endif
#define ALL_FLAG_TEST \
    if ((kpos >> v) <= D_KMAX) { \
      MFFONETEST my_factor_found |= 1; \
    }
#define ALL_FAST_FLAG_TEST \
  if ((((uint)(kpos >> 32))>>v) == 0) { \
    if(((uint)(kpos >> v)) <= D_KMAX) { \
      MFFONETEST my_factor_found |= 1; \
    } \
  }
#else
#define ALL_FLAG_TEST \
    if ((kpos >> v) <= d_kmax) { \
      if((kpos >> v) >= THE_KMIN) \
        my_factor_found |= 1; \
    }
#endif
#define VEC_IF if(v != 0)
#elif(VECSIZE == 2)
#define ALL_CTZLL \
      VEC_CTZLL(v,x) \
      VEC_CTZLL(v,y)
#define ALL_FLAG_TEST \
      VEC_FLAG_TEST(x) \
      VEC_FLAG_TEST(y)
#define ALL_FAST_FLAG_TEST \
      VEC_FAST_FLAG_TEST(x) \
      VEC_FAST_FLAG_TEST(y)
#define VEC_IF if(v.x != 0 && v.y != 0)
#elif(VECSIZE == 4)
#define ALL_CTZLL \
      VEC_CTZLL(v,x) \
      VEC_CTZLL(v,y) \
      VEC_CTZLL(v,z) \
      VEC_CTZLL(v,w)
#define ALL_FLAG_TEST \
      VEC_FLAG_TEST(x) \
      VEC_FLAG_TEST(y) \
      VEC_FLAG_TEST(z) \
      VEC_FLAG_TEST(w)
#define ALL_FAST_FLAG_TEST \
      VEC_FAST_FLAG_TEST(x) \
      VEC_FAST_FLAG_TEST(y) \
      VEC_FAST_FLAG_TEST(z) \
      VEC_FAST_FLAG_TEST(w)
#define VEC_IF if(v.x != 0 && v.y != 0 && v.z != 0 && v.w != 0)
#else
#error "Invalid vecsize" #VECSIZE
#endif
      // ALL_CTZLL happens, but my calculator can't compute how rarely!
      // About 1 in every 2^30 times.
#define ALL_IF_CTZLL \
    v = V2VINT(kpos); \
    VEC_IF { \
      v=31u - clz (v & -v); \
      ALL_FAST_FLAG_TEST \
    } else { \
      ALL_CTZLL \
      ALL_FLAG_TEST \
    }

#ifdef SEARCH_TWIN
// Select the even one here, so as to use the zero count and shift.
// The other side (whether positive or negative) is odd then, with no zeroes on the right.
#if(VECSIZE == 1)
#define TWIN_CHOOSE_EVEN_K0 kpos = (((unsigned int)k0) & 1)?(my_P - k0):k0;
#else
#define TWIN_CHOOSE_EVEN_K0 kpos = (((k0) & (VLONG)VECTORIZE(1)) == (VLONG)VECTORIZE(1))?(my_P - k0):k0;
#endif
//#define TWIN_CHOOSE_EVEN kpos = (((kpos) & (VLONG)VECTORIZE(1)) == (VLONG)VECTORIZE(1))?(my_P - kpos):kpos;
#else
#define TWIN_CHOOSE_EVEN_K0 kpos = k0;
//#define TWIN_CHOOSE_EVEN
#endif

__kernel void check_more_ns(__global ulong * P, __global ulong * Psarr, __global ulong * K, __global uint * factor_found_arr, const uint N, const uint shift
                                 // Device constants
#ifndef D_KMAX
#ifndef D_KMIN
                                 , const ulong d_kmin		// 6
#endif
                                 , const ulong d_kmax		// 7
#endif
                                 ) {
  VINT v;
  VLONG my_P;// = P[VDOTX];
  VLONG Ps;// = Psarr[VDOTX];
  uint n = N;
  VLONG k0;// = K[VDOTX];
  VINT my_factor_found;// = factor_found_arr[VDOTX];
  VLONG kpos;
  uint l_nmax = n + D_KERNEL_NSTEP;
  if(l_nmax > D_NMAX) l_nmax = D_NMAX;

  VDOTX = get_global_id(0);
  Ps = VLOAD(VDOTX, Psarr);
  k0 = VLOAD(VDOTX, K);
  my_P = VLOAD(VDOTX, P);
  if(n == D_NMIN) {
    my_factor_found = (VINT)VECTORIZE(0);
  } else {
    my_factor_found = VLOAD(VDOTX, factor_found_arr);
    my_factor_found <<= shift;
  }
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (VINT)k0, D_KERNEL_NSTEP);
#endif
  do { // Remaining steps are all of equal size nstep
    // The following macros use kpos.
    TWIN_CHOOSE_EVEN_K0
    //i = __ffsll(kpos)-1;
    // Just flag this if kpos <= d_kmax.
    ALL_IF_CTZLL

    n += D_NSTEP;
#if D_NSTEP <= 32
#if D_NSTEP == 32
    k0 = shiftmod_REDC32(k0, my_P, V2VINT(k0)*V2VINT(Ps));
#else // D_NSTEP < 32
    k0 = shiftmod_REDCsm(k0, my_P, V2VINT(k0)*V2VINT(Ps));
#endif
#else
    // Proceed to the K for the next N.
    // kpos is destroyed, just to keep the register count down.
    kpos = k0 * Ps;
    k0 = shiftmod_REDC(k0, my_P, kpos);
#endif
  } while (n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (VINT)k0);
#endif
  VDOTX = get_global_id(0);
  //factor_found_arr[VDOTX] = my_factor_found;
  if(n < D_NMAX) {
    //K[VDOTX] = k0;
    VSTORE(k0, VDOTX, K);
  } else {
    // Prepend a checksum of the K result to my_factor_found.
    my_factor_found |= V2VINT(k0)<<8;
  }
  VSTORE(my_factor_found, VDOTX, factor_found_arr);
}
