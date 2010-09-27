/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.cu -- (C) Ken Brazier August 2010.

   Proth Prime Search sieve OpenCL Kernel portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

//#define KERNEL_ONLY
//#pragma OPENCL EXTENSION cl_amd_printf : enable
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
  rcx = (rax!=(VLONG)VECTORIZE(0))?((VLONG)VECTORIZE(1)):((VLONG)VECTORIZE(0));	// if rax != 0, increase rcx
  //"mulq %[N]\n\t"           // rdx:rax = m * N 		Cycle 13?-19?
  rax = mad_hi(rax, N, rcx);
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
#define VLOAD VECSIZEIT(vload,VECSIZE)
#define VSTORE VECSIZEIT(vstore,VECSIZE)
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
  //factor_found_arr[i] = 0;
  VSTORE(0, i, factor_found_arr);
  //Ps[i] = my_Ps;
  VSTORE(my_Ps, i, Ps);
  //K[i] = k0;
  VSTORE(k0, i, K);
}

// Continue checking N's.
//i=(__float_as_int(__VINT2float_rz(i & -i))>>23)-0x7f;
#define VEC_CTZLL(I,X) \
        if(I.X != 0) { \
          I.X=31u - clz (I.X & -I.X); \
      } else { \
        I.X = (uint)(kpos.X>>32); \
        I.X=63 - clz (I.X & -I.X); \
      }
    // Just flag this if kpos <= d_kmax.
#ifdef D_KMAX
#define VEC_FLAG_TEST(X) \
    if ((((uint)(kpos.X >> 32))>>v.X) == 0) { \
     if(((uint)(kpos.X >> v.X)) <= D_KMAX) { \
      if((kpos.X >> v.X) >= D_KMIN) \
        my_factor_found.X = 1; \
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
      if((kpos.X >> v.X) >= THE_KMIN) \
        my_factor_found.X = 1; \
    }
#endif
#if(VECSIZE == 2)
#define ALL_CTZLL \
      VEC_CTZLL(v,x) \
      VEC_CTZLL(v,y)
#define ALL_FLAG_TEST \
      VEC_FLAG_TEST(x) \
      VEC_FLAG_TEST(y)
#define VEC_IF if(v.x == 0 || v.y == 0)
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
#define VEC_IF if(v.x == 0 || v.y == 0 || v.z == 0 || v.w == 0)
#else
#error "Invalid vecsize" #VECSIZE
#endif
      // ALL_CTZLL happens, but my calculator can't compute how rarely!
      // About 1 in every 2^30 times.
#define ALL_IF_CTZLL \
    v = V2VINT(kpos); \
    VEC_IF { \
      ALL_CTZLL \
    } else { \
      v=31u - clz (v & -v); \
    }



__kernel void check_more_ns(__global ulong * P, __global ulong * Psarr, __global ulong * K, __global uint * factor_found_arr, const uint N
                                 // Device constants
#ifndef D_KMIN
                                 , const ulong d_kmin		// 5
#endif
#ifndef D_KMAX
                                 , const ulong d_kmax		// 6
#endif
                                 ) {
  VINT v;
  VLONG my_P;// = P[v.x];
  VLONG Ps;// = Psarr[v.x];
  uint n = N;
  VLONG k0;// = K[v.x];
  VINT my_factor_found;// = factor_found_arr[v.x];
  VLONG kpos; //, kPs;
  uint l_nmax = n + D_KERNEL_NSTEP;
  if(l_nmax > D_NMAX) l_nmax = D_NMAX;

  v.x = get_global_id(0);
  Ps = VLOAD(v.x, Psarr);
  k0 = VLOAD(v.x, K);
  my_P = VLOAD(v.x, P);
  my_factor_found = VLOAD(v.x, factor_found_arr);
  
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Started at n=%u, k=%u; running %u n's (GPU)\n", n, (VINT)k0, D_KERNEL_NSTEP);
#endif
  do { // Remaining steps are all of equal size nstep
    // Get K from the Montgomery form.
    // This is equivalent to mod_REDC(k, my_P, Ps), but the intermediate kPs value is kept for later.
    kpos = k0;
    //i = __ffsll(kpos)-1;
    ALL_IF_CTZLL

    // Just flag this if kpos <= d_kmax.
    ALL_FLAG_TEST

    // Proceed to the K for the next N.
    // kpos is destroyed, just to keep the register count down.
    kpos = k0 * Ps;
    k0 = shiftmod_REDC(k0, my_P, kpos);
    n += D_NSTEP;
  } while (n < l_nmax);
#ifdef _DEVICEEMU
  //if(my_P == 42070000070587) printf("Stopped at n=%u, k=%u (GPU)\n", n, (VINT)k0);
#endif
  v.x = get_global_id(0);
  //factor_found_arr[v.x] = my_factor_found;
  VSTORE(my_factor_found, v.x, factor_found_arr);
  if(n < D_NMAX) {
    //K[v.x] = k0;
    VSTORE(k0, v.x, K);
  }
}
