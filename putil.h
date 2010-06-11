/* util.h -- (C) Geoffrey Reynolds, March 2009.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _UTIL_H
#define _UTIL_H 1

#include <stdlib.h>
#include "stdint.h"


/* Stringification macros */
#define XSTR(ARG) STR(ARG)
#define STR(ARG) #ARG

/* malloc-or-die functions */
void *xmalloc(size_t size);
void *xrealloc(void *mem, size_t size);
char *xstrdup(const char *str);

/* Allocating string functions */
#ifndef _GNU_SOURCE
int asprintf(char **out, const char *fmt, const char *str);
#endif
char* astrcpy(char **dest, const char *src);

int parse_uint(unsigned int *result, const char *str,
               unsigned int lo, unsigned int hi);
int parse_uint64(uint64_t *result, const char *str,
                 uint64_t lo, uint64_t hi);

#ifndef __GNUC__
unsigned long __builtin_ctzll(unsigned __int64 i);
#endif
#endif /* _UTIL_H */
