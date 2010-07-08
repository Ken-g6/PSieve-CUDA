appname=primegrid_ppsieve_1.00_x86_64-apple-darwin__cuda23
BOINC_DIR=../../../boinc/boinc-6.7.4
BOINC_API_DIR=$BOINC_DIR/api
BOINC_LIB_DIR=$BOINC_DIR/lib
BOINC_LOAD_LIBS="-I$BOINC_DIR -I$BOINC_LIB_DIR -I$BOINC_API_DIR -L$BOINC_DIR -L$BOINC_LIB_DIR -L$BOINC_API_DIR"

MACOSX_DEPLOYMENT_TARGET=10.5 nvcc -arch=sm_11 --ptxas-options=-v -O3 -m64 -Xcompiler "-arch x86_64 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5" -Xlinker "-syslibroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5" -D_GNU_SOURCE -DUSE_BOINC -DNDEBUG -D_REENTRANT -I. -I.. -o $appname $BOINC_LOAD_LIBS ../main.c ../sieve.c ../clock.c ../putil.c cuda_sleep_memcpy.cu appcu.cu app.c -lm -lpthread -lboinc_api -lboinc -L/usr/local/cuda/lib -lcuda
