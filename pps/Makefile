APPNAME = ppsieve-cuda
COMMON_OBJS = ../clock.o cuda_sleep_memcpy.o factor_proth.o
OBJS = ../main.o ../putil.o ../sieve.o app.o appcu.o $(COMMON_OBJS)
BOINC_OBJS = ../do_boinc_init.o ../main_boinc.o ../putil_boinc.o ../sieve_boinc.o app_boinc.o appcu_boinc.o $(COMMON_OBJS)
COMMON64_OBJS = ../clock64.o cuda_sleep_memcpy64.o factor_proth64.o
OBJS64 = ../main64.o ../putil64.o ../sieve64.o app64.o appcu64.o $(COMMON64_OBJS)
BOINC64_OBJS = ../do_boinc_init64.o ../main_boinc64.o ../putil_boinc64.o ../sieve_boinc64.o app_boinc64.o appcu_boinc64.o $(COMMON64_OBJS)
BOINC_BINS = $(APPNAME)-boinc-x86-linux $(APPNAME)-boinc-x86_64-linux
NON_BOINC_BINS = $(APPNAME)-x86-linux $(APPNAME)-x86_64-linux
BINS = $(NON_BOINC_BINS) $(BOINC_BINS)
CC = gcc
CC32 = $(CC) -c -m32 -march=i586 -mtune=core2
CC64 = $(CC) -c -m64 -march=k8 -mno-3dnow -mtune=core2
NVCC = nvcc

# Flags:
# Uncomment to use device emulation:
#DOEMU = -deviceemu -D_DEVICEEMU -g
COMMON_FLAGS = -O3 $(DOEMU) -DNDEBUG -D_REENTRANT -I. -I..
CFLAGS = -fomit-frame-pointer -s $(COMMON_FLAGS)
NVCFLAGS = --ptxas-options=-v $(COMMON_FLAGS)
LINK = -lm -lpthread # -lcuda
BOINC64_LINK = $(LINK) -lboinc_api -lboinc `g++ -m64 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
BOINC_LINK = $(LINK) -lboinc_api -lboinc `g++ -m32 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
#nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DUSE_BOINC -DNDEBUG -D_REENTRANT -m32 -I. -I.. -o $appname-boinc-x86-linux $BOINC_LOAD_LIBS ../main.c ../sieve.c ../clock.c ../putil.c ../do_boinc_init.cpp cuda_sleep_memcpy.cu appcu.cu app.c factor_proth.c $link -lboinc_api -lboinc `g++ -print-file-name=libstdc++.a` -DAPP_GRAPHICS

# Make-able targets:
all: $(BINS)

# Make bins before zips, but making zips is too complicated for this Makefile (so far.)
# Should also have $(APPNAME)-x86-windows.exe $(APPNAME)-boinc-x86-windows.exe, but those can't be made on Linux. (Yet.)
zips: $(BINS)
	bash ./make-zips.sh

non-boinc: $(NON_BOINC_BINS)

boinc: $(BOINC_BINS)

clean:
	rm *.o ../*.o $(BINS)

$(APPNAME)-x86-linux: $(OBJS)
	$(NVCC) -m32 -o $@ $(NVCFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86-linux: $(BOINC_OBJS)
	$(NVCC) -m32 -o $@ $(NVCFLAGS) $^ $(BOINC_LINK)

$(APPNAME)-x86_64-linux: $(OBJS64)
	$(NVCC) -m64 -o $@ $(NVCFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86_64-linux: $(BOINC64_OBJS)
	$(NVCC) -m64 -o $@ $(NVCFLAGS) $^ $(BOINC64_LINK)

# OBJS
../main.o: ../main.c ../main.h
	$(CC32) -o $@ $(CFLAGS) $<

../putil.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

../sieve.o: ../sieve.c ../sieve.h
	$(CC32) -o $@ $(CFLAGS) $<

app.o: app.c app.h
	$(CC32) -o $@ $(CFLAGS) $<

appcu.o: appcu.cu appcu.h
	$(NVCC) -c -m32 -o $@ $(NVCFLAGS) $<


# BOINC_OBJS
../do_boinc_init.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(NVCC) -c -m32 -o $@ $(NVCFLAGS) $<

../main_boinc.o: ../main.c ../main.h
	$(CC32) -o $@ $(CFLAGS) $<

../putil_boinc.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

../sieve_boinc.o: ../sieve.c ../sieve.h
	$(CC32) -o $@ $(CFLAGS) $<

app_boinc.o: app.c app.h
	$(CC32) -o $@ $(CFLAGS) $<

appcu_boinc.o: appcu.cu appcu.h
	$(NVCC) -c -m32 -o $@ $(NVCFLAGS) $<

#  COMMON_OBJS
../clock.o: ../clock.c ../clock.h
	$(CC32) -o $@ $(CFLAGS) $<

factor_proth.o: factor_proth.c factor_proth.h
	$(CC32) -o $@ $(CFLAGS) $<

cuda_sleep_memcpy.o: cuda_sleep_memcpy.cu cuda_sleep_memcpy.h
	$(NVCC) -c -m32 -o $@ $(NVCFLAGS) $<

# OBJS64
../main64.o: ../main.c ../main.h
	$(CC64) -o $@ $(CFLAGS) $<

../putil64.o: ../putil.c ../putil.h
	$(CC64) -o $@ $(CFLAGS) $<

../sieve64.o: ../sieve.c ../sieve.h
	$(CC64) -o $@ $(CFLAGS) $<

app64.o: app.c app.h
	$(CC64) -o $@ $(CFLAGS) $<

appcu64.o: appcu.cu appcu.h
	$(NVCC) -c -m64 -o $@ $(NVCFLAGS) $<


# BOINC64_OBJS
../do_boinc_init64.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(NVCC) -c -m64 -o $@ $(NVCFLAGS) $<

../main_boinc64.o: ../main.c ../main.h
	$(CC64) -o $@ $(CFLAGS) $<

../putil_boinc64.o: ../putil.c ../putil.h
	$(CC64) -o $@ $(CFLAGS) $<

../sieve_boinc64.o: ../sieve.c ../sieve.h
	$(CC64) -o $@ $(CFLAGS) $<

app_boinc64.o: app.c app.h
	$(CC64) -o $@ $(CFLAGS) $<

appcu_boinc64.o: appcu.cu appcu.h
	$(NVCC) -c -m64 -o $@ $(NVCFLAGS) $<


#  COMMON64_OBJS
../clock64.o: ../clock.c ../clock.h
	$(CC64) -o $@ $(CFLAGS) $<

factor_proth64.o: factor_proth.c factor_proth.h
	$(CC64) -o $@ $(CFLAGS) $<

cuda_sleep_memcpy64.o: cuda_sleep_memcpy.cu cuda_sleep_memcpy.h
	$(NVCC) -c -m64 -o $@ $(NVCFLAGS) $<
