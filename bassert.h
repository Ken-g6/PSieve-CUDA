// BOINC assert

#include "main.h"
#ifndef assert

#define assert(exp)  if(!(exp)) { fprintf(stderr, "Assertion failed: %s\n", #exp); bexit(-193); }

#endif
