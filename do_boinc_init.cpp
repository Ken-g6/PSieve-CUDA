/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* do_init_boinc.c -- (C) Ken Brazier July 2010

   Initialize BOINC in C++, for the C code elsewhere.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
#ifdef USE_BOINC
#if defined(__APPLE__) || defined(_WIN32)
#include "filesys.h"
#endif
#include "boinc_api.h"
#include "diagnostics.h"     // boinc_init_diagnostics()
#include "main.h"
#include "do_boinc_init.h"

int do_boinc_init() {
  int retval;
  BOINC_OPTIONS options;
  if (!diagnostics_is_initialized()) {
    retval = boinc_init_diagnostics(BOINC_DIAG_REDIRECTSTDERR|
        BOINC_DIAG_MEMORYLEAKCHECKENABLED|
        BOINC_DIAG_DUMPCALLSTACKENABLED| 
        BOINC_DIAG_TRACETOSTDERR);
    if (retval) return retval;
  }
  boinc_options_defaults(options);
#ifdef _WIN32
  if(priority_opt < 14) options.normal_thread_priority = 1;
#endif
  return boinc_init_options(&options);
}
#endif
