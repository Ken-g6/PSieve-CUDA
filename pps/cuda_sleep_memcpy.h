/* cuda_sleep_memcpy.h -- (C) Ken Brazier February 2010.

   Helps a consistent CUDA kernel not busy-wait so long.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifndef _CUDA_SLEEP_MEMCPY_H
#define CUDA_SLEEP_MEMCPY_H 1

cudaError_t cudaSleepMemcpy(void *dst, const void *src, size_t count, enum cudaMemcpyKind kind, int *delay, const int overlap);
cudaError_t cudaSleepMemcpyFromTime(void *dst, const void *src, size_t count, enum cudaMemcpyKind kind, int *delay, const int overlap, const uint64_t start_t);
cudaError_t cudaSleepWait(cudaEvent_t &stop, int *delay, const int overlap, const uint64_t start_t);
#endif
