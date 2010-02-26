appname=ppsieve-cuda
DOEMU=-deviceemu
#cleanvars="-fomit-frame-pointer -s"
#cleanvars="-s"
if uname -a | grep " 2.6.28[^ ]* [^ a-zA-Z]*Ubuntu" > /dev/null ; then 
	kernel=2.6.15
	#cleanvars="-s -fno-stack-protector"
fi
if [ "$kernel" != "" ] ; then export LD_ASSUME_KERNEL=$kernel ; fi
#
nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DNDEBUG -D_REENTRANT -m32 -I. -I.. -o $appname-x86-linux ../main.c ../sieve.c ../clock.c ../util.c cuda_sleep_memcpy.cu appcu.cu app.c -lm -lpthread
#nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DNDEBUG -D_REENTRANT -march=pentium-m -msse2 -I. -I.. -o $appname-x86-linux-sse2 ../main.c ../sieve.c ../clock.c ../util.c app.c -lm -lpthread
if [ "$kernel" != "" ] ; then unset LD_ASSUME_KERNEL ; fi
#if [ -f /multimedia/mingw/bin/gcc.exe ] ; then
	#wine /multimedia/mingw/bin/gcc.exe -O3 -fomit-frame-pointer -DNDEBUG -D_REENTRANT -m32 -march=i686 -mtune=pentium3 -I. -I.. -s -o $appname-x86-windows.exe ../main.c ../sieve.c ../clock.c ../util.c app.c -lm
	#wine /multimedia/mingw/bin/gcc.exe -O3 -fomit-frame-pointer -DNDEBUG -D_REENTRANT -m32 -march=i686 -mtune=core2 -msse2 -I. -I.. -s -o $appname-x86-windows-sse2.exe ../main.c ../sieve.c ../clock.c ../util.c app.c -lm
#fi
#
if [ "$kernel" != "" ] ; then export LD_ASSUME_KERNEL=$kernel ; fi
#gcc -O3 $cleanvars -DNDEBUG -D_REENTRANT -m64 -march=k8 -mno-3dnow -mtune=core2 -I. -I.. -o $appname-x86_64-linux ../main.c ../sieve.c ../clock.c ../util.c app.c -lm -lpthread
# --ptxas-options=-v : Print the number of registers used by each kernel.
# -deviceemu : compile for the CPU emulator.
nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DBUSYWAIT -DBITSATATIME=3 -DNDEBUG -D_REENTRANT -m64 -I. -I.. -o $appname-24bit-x86_64-linux ../main.c ../sieve.c ../clock.c ../util.c cuda_sleep_memcpy.cu appcu.cu app.c -lm -lpthread
nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DBUSYWAIT -DNDEBUG -D_REENTRANT -m64 -I. -I.. -o $appname-64bit-x86_64-linux ../main.c ../sieve.c ../clock.c ../util.c cuda_sleep_memcpy.cu appcu.cu app.c -lm -lpthread
nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DNDEBUG -D_REENTRANT -m64 -I. -I.. -o $appname-x86_64-linux ../main.c ../sieve.c ../clock.c ../util.c cuda_sleep_memcpy.cu appcu.cu app.c -lm -lpthread
if [ "$kernel" != "" ] ; then unset LD_ASSUME_KERNEL ; fi
#
#if [ -f /multimedia/mingw64/bin/x86_64-w64-mingw32-gcc.exe ] ; then
	#wine /multimedia/mingw64/bin/x86_64-w64-mingw32-gcc.exe -O3 -fomit-frame-pointer -DNDEBUG -D_REENTRANT -m64 -march=k8 -mno-3dnow -mtune=core2 -I. -I.. -s -o $appname-x86_64-windows.exe ../main.c ../sieve.c ../clock.c ../util.c app.c -lm
#fi
#
#rm -f $appname-$appver-bin.zip
#zip -9 $appname-$appver-bin.zip README.txt CHANGES.txt ppconfig.txt $appname-x86-linux $appname-x86-linux-sse2 $appname-x86_64-linux $appname-x86-windows.exe $appname-x86-windows-sse2.exe $appname-x86_64-windows.exe make-bins.sh make-bins.bat app.c app.h ../src-0.2.4.zip
