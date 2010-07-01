/* clock.h -- (C) Geoffrey Reynolds, March 2009.

   Elapsed and processor times in microsecond units.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CLOCK_H
#define _CLOCK_H 1

#include "stdint.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Time elapsed since some fixed base time. */
uint64_t elapsed_usec(void);

/* Total processor time consumed by all threads of the current process. */
uint64_t processor_usec(void);

/* CPU cycles */
uint64_t processor_cycles(void);
#ifdef __cplusplus
}
#endif
#endif /* _CLOCK_H */
