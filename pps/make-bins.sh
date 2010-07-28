appname=ppsieve-cuda
#DOEMU="-deviceemu -D_DEVICEEMU -g"
# Test: ./ppsieve-cuda-x86_64-linux -p42070e9 -P42070003e6 -k 1201 -K 9999 -N 2000000 -z normal
#cleanvars="-fomit-frame-pointer -s"
#cleanvars="-s"
# To build the BOINC Linux 64 version, get and make BOINC from subversion.  Place its directory location below.
#BOINC_DIR=/downloads/distributed/boinc610/server_stable
#BOINC_DIR=../../ppsb/boinc
BOINC_API_DIR=$BOINC_DIR/api
BOINC_LIB_DIR=$BOINC_DIR/lib
arch=`uname -m`
link='-lm -lpthread' # -lcuda'
#if [ "$DOEMU" == "" ] ; then link="$link -lcuda" ; fi

if [ "$BOINC_DIR" != "" ] ; then
	BOINC_LOAD_LIBS="-I$BOINC_DIR -I$BOINC_LIB_DIR -I$BOINC_API_DIR -L$BOINC_DIR -L$BOINC_LIB_DIR -L$BOINC_API_DIR"
else
	BOINC_LOAD_LIBS=""
fi
if uname -a | grep " 2.6.28[^ ]* [^ a-zA-Z]*Ubuntu" > /dev/null ; then 
	kernel=2.6.15
	#cleanvars="-s -fno-stack-protector"
fi

# 32-bit Linux (BOINC or non-BOINC)
if [ "$kernel" != "" ] ; then export LD_ASSUME_KERNEL=$kernel ; fi
if [ "$1" != "boinc" ] ; then
	nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DNDEBUG -D_REENTRANT -m32 -I. -I.. -o $appname-x86-linux ../main.c ../sieve.c ../clock.c ../putil.c cuda_sleep_memcpy.cu appcu.cu app.c $link
#nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DNDEBUG -D_REENTRANT -march=pentium-m -msse2 -I. -I.. -o $appname-x86-linux-sse2 ../main.c ../sieve.c ../clock.c ../putil.c app.c $link
else
	if [ "$arch" == "i686" ] ; then
		echo Making i686 BOINC version.
		nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DUSE_BOINC -DNDEBUG -D_REENTRANT -m32 -I. -I.. -o $appname-boinc-x86-linux $BOINC_LOAD_LIBS ../main.c ../sieve.c ../clock.c ../putil.c ../do_boinc_init.cpp cuda_sleep_memcpy.cu appcu.cu app.c $link -lboinc_api -lboinc `g++ -print-file-name=libstdc++.a` -DAPP_GRAPHICS
	fi
fi
if [ "$kernel" != "" ] ; then unset LD_ASSUME_KERNEL ; fi
#if [ -f /multimedia/mingw/bin/gcc.exe ] ; then
	#wine /multimedia/mingw/bin/gcc.exe -O3 -fomit-frame-pointer -DNDEBUG -D_REENTRANT -m32 -march=i686 -mtune=pentium3 -I. -I.. -s -o $appname-x86-windows.exe ../main.c ../sieve.c ../clock.c ../putil.c app.c -lm
	#wine /multimedia/mingw/bin/gcc.exe -O3 -fomit-frame-pointer -DNDEBUG -D_REENTRANT -m32 -march=i686 -mtune=core2 -msse2 -I. -I.. -s -o $appname-x86-windows-sse2.exe ../main.c ../sieve.c ../clock.c ../putil.c app.c -lm
#fi
#
if [ "$kernel" != "" ] ; then export LD_ASSUME_KERNEL=$kernel ; fi
#gcc -O3 $cleanvars -DNDEBUG -D_REENTRANT -m64 -march=k8 -mno-3dnow -mtune=core2 -I. -I.. -o $appname-x86_64-linux ../main.c ../sieve.c ../clock.c ../putil.c app.c -lm -lpthread
# --ptxas-options=-v : Print the number of registers used by each kernel.
# -deviceemu : compile for the CPU emulator.
if [ "$1" != "boinc" ] ; then
nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DNDEBUG -D_REENTRANT -m64 -I. -I.. -o $appname-x86_64-linux ../main.c ../sieve.c ../clock.c ../putil.c cuda_sleep_memcpy.cu appcu.cu app.c $link
else
	if [ "$arch" == "x86_64" ] ; then
		echo Making x86_64 BOINC version.
		nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DUSE_BOINC -DNDEBUG -D_REENTRANT -m64 -I. -I.. -o $appname-boinc-x86_64-linux $BOINC_LOAD_LIBS ../main.c ../sieve.c ../clock.c ../putil.c ../do_boinc_init.cpp cuda_sleep_memcpy.cu appcu.cu app.c $link -lboinc_api -lboinc `g++ -print-file-name=libstdc++.a` -DAPP_GRAPHICS
	fi
fi
if [ "$kernel" != "" ] ; then unset LD_ASSUME_KERNEL ; fi
#
#if [ -f /multimedia/mingw64/bin/x86_64-w64-mingw32-gcc.exe ] ; then
	#wine /multimedia/mingw64/bin/x86_64-w64-mingw32-gcc.exe -O3 -fomit-frame-pointer -DNDEBUG -D_REENTRANT -m64 -march=k8 -mno-3dnow -mtune=core2 -I. -I.. -s -o $appname-x86_64-windows.exe ../main.c ../sieve.c ../clock.c ../putil.c app.c -lm
#fi
