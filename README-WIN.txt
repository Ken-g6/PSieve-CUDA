Building TPSieve-cl on Win32 with Visual Studio Express Edition 2008:

0. Install the AMD APP SDK v2.9 32-bit.  (Yes, even on a 64-bit computer.)
1. Extract or Git the source code to a new folder.
2. (BOINC only) Create a subfolder of that folder called "boinc"
3. (BOINC only) Place the BOINC source code there.  You may also need libssl.
4. (BOINC only) Build libboinc_staticcrt and libboincapi_staticcrt under win_build, with Visual Studio.  Be sure to build the Release versions.
5. In the pps subfolder of the main folder, run convert_visual_c.bat.  This moves files to their correct locations, and loads the Solution.
6. Build the TPS-Release version for a standalone app, or the TPS-BOINC-Release version.
