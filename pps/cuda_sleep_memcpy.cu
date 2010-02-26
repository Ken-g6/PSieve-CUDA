/* cuda_sleep_memcpy.cu -- (C) Ken Brazier February 2010.

   Helps a consistent CUDA kernel not busy-wait so long.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifndef NDEBUG
#include <stdio.h>
#endif
#include <unistd.h>
#include "clock.h"

// Allows one to perform a CUDA memcpy without (as much) busy-wait.
// For best performance, both the CUDA kernel and any CPU code after its launch
// should take consistent amounts of time.
//
// * Delay is time to sleep in usec, passed by reference.
//   It is dynamically adjusted by the function.
//
// * Overlap is time that CUDA should busy-wait, in usec.
//   Larger overlaps recover from extra-long kernel calls more quickly,
//   but they use more CPU all other times.
cudaError_t cudaSleepMemcpy(void *dst, const void *src, size_t count, enum cudaMemcpyKind kind, int *delay, const int overlap) {
	cudaError_t ret;
#ifndef BUSYWAIT
	// Timer variables.
	uint64_t start_t, busy_wait_t;

	// First, sleep as long as seemed a good idea last time.
#ifndef NDEBUG
	fprintf(stderr, "Sleeping %d usec.\n", *delay);
#endif
	usleep(*delay);

	// Now, time the busy-wait.
	start_t = elapsed_usec();
#endif
	ret = cudaMemcpy(dst, src, count, kind);
#ifndef BUSYWAIT
	busy_wait_t = elapsed_usec()-start_t;
	
	// Add that to the delay.
	*delay += busy_wait_t;

	// Subtract the expected overlap.
	*delay -= overlap;

	// Can't sleep less than 0 seconds.
	if(*delay < 0) *delay = 0;

#ifndef NDEBUG
	fprintf(stderr, "Will sleep %d usec next time.\n", *delay);
#endif
#endif
	return ret;
}

