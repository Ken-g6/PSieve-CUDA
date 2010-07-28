/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* do_init_boinc.h -- (C) Ken Brazier July 2010

   Initialize BOINC in C++, for the C code elsewhere.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifdef USE_BOINC
#ifndef _DO_BOINC_INIT
#define _DO_BOINC_INIT

#ifdef __cplusplus
extern "C" {
#endif
int do_boinc_init();
#ifdef __cplusplus
}
#endif

#endif
#endif
