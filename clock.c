/* clock.c -- (C) Geoffrey Reynolds, March 2009.

   Elapsed and processor times in microsecond units.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifdef _WIN32
# include <windows.h>
#include "stdint.h"
#else
#include <stdint.h>
# include <sys/resource.h>
# include <sys/time.h>
#endif
#include "clock.h"


/* Time elapsed since some fixed base time. */
uint64_t elapsed_usec(void)
{
#ifdef _WIN32
  FILETIME ft;
  ULARGE_INTEGER ns100;
  GetSystemTimeAsFileTime(&ft);
  ns100.u.LowPart = ft.dwLowDateTime;
  ns100.u.HighPart = ft.dwHighDateTime;
  return (uint64_t)ns100.QuadPart/10;
#else
  struct timeval t;
  gettimeofday(&t,0);
  return (uint64_t)t.tv_sec*1000000 + t.tv_usec;
#endif
}

/* Total processor time consumed by all threads of the current process. */
uint64_t processor_usec(void)
{
#ifdef _WIN32
  FILETIME ft_create, ft_exit, ft_kernel, ft_user;
  ULARGE_INTEGER ns100_kernel, ns100_user;
  GetProcessTimes(GetCurrentProcess(),&ft_create,&ft_exit,&ft_kernel,&ft_user);
  ns100_kernel.u.LowPart = ft_kernel.dwLowDateTime;
  ns100_kernel.u.HighPart = ft_kernel.dwHighDateTime;
  ns100_user.u.LowPart = ft_user.dwLowDateTime;
  ns100_user.u.HighPart = ft_user.dwHighDateTime;
  return (uint64_t)(ns100_kernel.QuadPart + ns100_user.QuadPart)/10;
#else
  struct rusage r;
  getrusage(RUSAGE_SELF,&r);
  return (uint64_t)(r.ru_utime.tv_sec + r.ru_stime.tv_sec)*1000000
    + (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
#endif
}

/* CPU cycles */
uint64_t processor_cycles(void)
{
  uint64_t cycles;

#if defined(__GNUC__) && defined(__x86_64__)
  asm volatile ("rdtsc"                 "\n\t"
                "shl     $32, %%rdx"    "\n\t"
                "or      %%rdx, %%rax"
                : "=a" (cycles) : : "%rdx", "cc" );

#elif defined(_MSC_VER) && defined(_M_X64)
  __asm { rdtsc
          shl     rdx, 32
          or      rax, rdx
          mov     qword ptr cycles, rax }

#elif defined(__GNUC__) && defined(__i386__)
  asm volatile ("rdtsc"
                : "=A" (cycles) : : "cc" );

#elif defined(_MSC_VER) && defined(_M_IX86)
  __asm { rdtsc
          mov     dword ptr cycles, eax
          mov     dword ptr cycles+4, edx }

#elif defined(__ppc64__) || defined(__powerpc64__)
  /* TODO: how to scale result to get clock cycles? */
  asm volatile ("mftb %0"
                : "=r" (cycles) );

#else
  cycles = 0;

#endif

  return cycles;
}
