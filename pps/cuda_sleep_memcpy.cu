/* cuda_sleep_memcpy.cu -- (C) Ken Brazier February 2010.

   Helps a consistent CUDA kernel not busy-wait so long.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifdef TRACE
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
#ifdef TRACE
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

#ifdef TRACE
	fprintf(stderr, "Will sleep %d usec next time.\n", *delay);
#endif
#endif
	return ret;
}

// Allows one to perform a CUDA memcpy without (as much) busy-wait.
// This version expects elapsed_usec() to be called immediately
// after the kernel call; then passed in here later.
// For best performance, the CUDA kernel should take a consistent
// amounts of time, but the CPU time doesn't matter.
//
// * Delay is time to sleep in usec, passed by reference.
//   It is dynamically adjusted by the function.
//
// * Overlap is time that CUDA should busy-wait, in usec.
//   Larger overlaps recover from extra-long kernel calls more quickly,
//   but they use more CPU all other times.
//
// * Start_t is the time the kernel was started, retrieved from
//   elapsed_usec().

cudaError_t cudaSleepMemcpyFromTime(void *dst, const void *src, size_t count, enum cudaMemcpyKind kind, int *delay, const int overlap, const uint64_t start_t) {
	cudaError_t ret;
#ifndef BUSYWAIT
	// Timer variables.
	int busy_wait_t;

	// First, sleep as long as seemed a good idea last time.
#ifdef TRACE
	fprintf(stderr, "Sleeping %d usec.\n", *delay);
#endif
	// Figure out how much time is left before a result is expected.
	if(*delay > 0) {
		busy_wait_t = *delay - (int)(elapsed_usec()-start_t);
		// Sleep that long.
		if(busy_wait_t > 0) usleep(busy_wait_t);
	}

	// Now, time the busy-wait.
#endif
	ret = cudaMemcpy(dst, src, count, kind);
#ifndef BUSYWAIT
	// Set that to the delay.
	*delay = elapsed_usec()-start_t;

	// Subtract the expected overlap.
	*delay -= overlap;

	// Can't sleep less than 0 seconds.
	if(*delay < 0) *delay = 0;

#ifdef TRACE
	fprintf(stderr, "Will sleep %d usec next time.\n", *delay);
#endif
#endif
	return ret;
}

