APPNAME = ppsieve-cuda
CL_APPNAME = ppsieve-cl
APPPATH = pps
COMMON_OBJS = ../clock.o factor_proth.o
OBJS = ../main.o ../putil.o ../sieve.o $(COMMON_OBJS)
CUDA_OBJS = $(OBJS) app.o cuda_sleep_memcpy.o appcu.o
CL_OBJS = $(OBJS) cl_app.o appcl.o
BOINC_OBJS = ../do_boinc_init.o ../main_boinc.o ../putil_boinc.o ../sieve_boinc.o $(COMMON_OBJS)
CUDA_BOINC_OBJS = $(BOINC_OBJS) app_boinc.o cuda_sleep_memcpy.o appcu_boinc.o
CL_BOINC_OBJS = $(BOINC_OBJS) cl_app_boinc.o appcl.o
COMMON64_OBJS = ../clock64.o factor_proth64.o
OBJS64 = ../main64.o ../putil64.o ../sieve64.o $(COMMON64_OBJS)
CUDA_OBJS64 = $(OBJS64) app64.o cuda_sleep_memcpy64.o appcu64.o
CL_OBJS64 = $(OBJS64) appcl64.o cl_app64.o
BOINC_OBJS64 = ../do_boinc_init64.o ../main_boinc64.o ../putil_boinc64.o ../sieve_boinc64.o $(COMMON64_OBJS)
CUDA_BOINC_OBJS64 = $(BOINC_OBJS64) app_boinc64.o cuda_sleep_memcpy64.o appcu_boinc64.o
CL_BOINC_OBJS64 = $(BOINC_OBJS64) cl_app_boinc64.o appcl64.o
BOINC_BINS = $(APPNAME)-boinc-x86-linux $(APPNAME)-boinc-x86_64-linux
NON_BOINC_BINS = $(APPNAME)-x86-linux $(APPNAME)-x86_64-linux
CL_BOINC_BINS = $(CL_APPNAME)-boinc-x86-linux $(CL_APPNAME)-boinc-x86_64-linux
CL_NON_BOINC_BINS = $(CL_APPNAME)-x86-linux $(CL_APPNAME)-x86_64-linux
# BINS should also have $(APPNAME)-x86-windows.exe $(APPNAME)-boinc-x86-windows.exe, but those can't be made on Linux. (Yet.)
BINS = $(NON_BOINC_BINS) $(BOINC_BINS)
CL_BINS = $(CL_NON_BOINC_BINS) $(CL_BOINC_BINS)

CC = gcc
CC32 = $(CC) -c -m32 -march=i586 -mtune=core2
CC64 = $(CC) -c -m64 -march=k8 -mno-3dnow -mtune=core2
NVCC = nvcc
CL_DEFINE = -DUSE_OPENCL
CL_INCLUDE = -I $(ATISTREAMSDKROOT)/include
CL32 = $(CXX) -c -m32 -march=i586 -mtune=core2 $(CL_INCLUDE)
CL64 = $(CXX) -c -m64 -march=k8 -mno-3dnow -mtune=core2 $(CL_INCLUDE)

# Flags:
# Uncomment to use device emulation:
#DOEMU = -deviceemu -D_DEVICEEMU -g
COMMON_FLAGS = -O3 -DNDEBUG -D_REENTRANT -I. -I..
CFLAGS = -Wall -fomit-frame-pointer -s $(COMMON_FLAGS)
NVCFLAGS = --ptxas-options=-v $(DOEMU) $(COMMON_FLAGS)
LINK = -lm -lpthread # -lcuda
CL_LINK = -lOpenCL -L$(ATISTREAMSDKROOT)/lib/x86
CL_LINK64 = -lOpenCL -L$(ATISTREAMSDKROOT)/lib/x86_64

# Uncomment this line if BOINC is in a user-specified directory.
#BOINC_DIR = /downloads/distributed/boinc610/server_stable
ifdef BOINC_DIR
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_LOAD_LIBS = -I$(BOINC_DIR) -I$(BOINC_LIB_DIR) -I$(BOINC_API_DIR)
BOINC_LINK_LIBS = -L$(BOINC_DIR) -L$(BOINC_LIB_DIR) -L$(BOINC_API_DIR)
else
# This line includes where some repositories put the BOINC includes:
BOINC_LOAD_LIBS = -I/usr/include/BOINC
endif
BOINC_CFLAGS = -DUSE_BOINC $(CFLAGS) $(BOINC_LOAD_LIBS)
BOINC_NVCFLAGS = -DUSE_BOINC $(NVCFLAGS) $(BOINC_LOAD_LIBS)

BOINC_LINK64 = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc `g++ -m64 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
BOINC_LINK = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc `g++ -m32 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
#nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DUSE_BOINC -DNDEBUG -D_REENTRANT -m32 -I. -I.. -o $appname-boinc-x86-linux $BOINC_LOAD_LIBS ../main.c ../sieve.c ../clock.c ../putil.c ../do_boinc_init.cpp cuda_sleep_memcpy.cu appcu.cu app.c factor_proth.c $link -lboinc_api -lboinc `g++ -print-file-name=libstdc++.a` -DAPP_GRAPHICS

# Make-able targets:
all: $(BINS)

# Making zipfiles (requires binaries):
zips: $(APPNAME).zip $(APPNAME)-src.zip

$(APPNAME).zip: $(BINS) gitstatus changes
	-rm $@ ; zip -9 $(APPNAME).zip $(APPNAME)-x86-linux $(APPNAME)-boinc-x86-linux $(APPNAME)-x86_64-linux $(APPNAME)-boinc-x86_64-linux $(APPNAME)-boinc-x86-windows.exe $(APPNAME)-x86-windows.exe $(APPNAME)-boinc-x86_64-windows.exe $(APPNAME)-x86_64-windows.exe README.txt CHANGES.txt $(APPNAME).sh $(APPNAME).bat *config.txt license.txt

$(APPNAME)-src.zip: $(BINS) gitstatus changes
	-rm $@ ; cd .. ; zip -9 $(APPPATH)/$(APPNAME)-src.zip Makefile $(APPPATH)/Makefile *.txt *.[ch] *.cpp *.sln $(APPPATH)/*.[ch] $(APPPATH)/*.cu $(APPPATH)/README.txt $(APPPATH)/CHANGES.txt $(APPPATH)/*config.txt $(APPPATH)/make-*.sh $(APPPATH)/$(APPNAME).sh $(APPPATH)/*.bat vc/*

gitstatus:
	git status | grep nothing > /dev/null || git status

changes:
	git log | grep "^\(    \|Date:\)" > CHANGES.txt

# Binaries:
non-boinc: $(NON_BOINC_BINS)

boinc: $(BOINC_BINS)

clean:
	-rm *.o ../*.o $(BINS) $(CL_BINS) $(APPNAME).zip $(APPNAME)-src.zip

$(APPNAME)-x86-linux: $(CUDA_OBJS)
	$(NVCC) -m32 -o $@ $(NVCFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86-linux: $(CUDA_BOINC_OBJS)
	$(NVCC) -m32 -o $@ $(BOINC_NVCFLAGS) $^ $(BOINC_LINK)

$(APPNAME)-x86_64-linux: $(CUDA_OBJS64)
	$(NVCC) -m64 -o $@ $(NVCFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86_64-linux: $(CUDA_BOINC_OBJS64)
	$(NVCC) -m64 -o $@ $(BOINC_NVCFLAGS) $^ $(BOINC_LINK64)

$(CL_APPNAME)-x86-linux: $(CL_OBJS)
	$(CXX) -m32 -o $@ $(CFLAGS) $^ $(CL_LINK) $(LINK)

$(CL_APPNAME)-boinc-x86-linux: $(CL_BOINC_OBJS)
	$(CXX) -m32 -o $@ $(BOINC_CFLAGS) $^ $(CL_LINK) $(BOINC_LINK)

$(CL_APPNAME)-x86_64-linux: $(CL_OBJS64)
	$(CXX) -m64 -o $@ $(CFLAGS) $^ $(CL_LINK64) $(LINK)

$(CL_APPNAME)-boinc-x86_64-linux: $(CL_BOINC_OBJS64)
	$(CXX) -m64 -o $@ $(BOINC_CFLAGS) $^ $(CL_LINK64) $(BOINC_LINK64)

# OBJS
../main.o: ../main.c ../main.h
	$(CC32) -o $@ $(CFLAGS) $<

../putil.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

../sieve.o: ../sieve.c ../sieve.h
	$(CC32) -o $@ $(CFLAGS) $<

app.o: app.c app.h
	$(CC32) -o $@ $(CFLAGS) $<

cl_app.o: app.c app.h
	$(CC32) -o $@ $(CFLAGS) $(CL_DEFINE) $<

appcu.o: appcu.cu appcu.h
	$(NVCC) -c -m32 -o $@ $(NVCFLAGS) $<


# BOINC_OBJS
../do_boinc_init.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(NVCC) -c -m32 -o $@ $(BOINC_NVCFLAGS) $<

../main_boinc.o: ../main.c ../main.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

../putil_boinc.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

../sieve_boinc.o: ../sieve.c ../sieve.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

app_boinc.o: app.c app.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

cl_app_boinc.o: app.c app.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $(CL_DEFINE) $<

appcu_boinc.o: appcu.cu appcu.h
	$(NVCC) -c -m32 -o $@ $(BOINC_NVCFLAGS) $<

#  COMMON_OBJS
../clock.o: ../clock.c ../clock.h
	$(CC32) -o $@ $(CFLAGS) $<

factor_proth.o: factor_proth.c factor_proth.h
	$(CC32) -o $@ $(CFLAGS) $<

cuda_sleep_memcpy.o: cuda_sleep_memcpy.cu cuda_sleep_memcpy.h
	$(NVCC) -c -m32 -o $@ $(NVCFLAGS) $<

appcl.o: appcl.cpp appcl.h
	$(CL32) -o $@ $(CFLAGS) $<

# OBJS64
../main64.o: ../main.c ../main.h
	$(CC64) -o $@ $(CFLAGS) $<

../putil64.o: ../putil.c ../putil.h
	$(CC64) -o $@ $(CFLAGS) $<

../sieve64.o: ../sieve.c ../sieve.h
	$(CC64) -o $@ $(CFLAGS) $<

app64.o: app.c app.h
	$(CC64) -o $@ $(CFLAGS) $<

cl_app64.o: app.c app.h
	$(CC64) -o $@ $(CFLAGS) $(CL_DEFINE) $<

appcu64.o: appcu.cu appcu.h
	$(NVCC) -c -m64 -o $@ $(NVCFLAGS) $<


# BOINC_OBJS64
../do_boinc_init64.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(NVCC) -c -m64 -o $@ $(BOINC_NVCFLAGS) $<

../main_boinc64.o: ../main.c ../main.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $<

../putil_boinc64.o: ../putil.c ../putil.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $<

../sieve_boinc64.o: ../sieve.c ../sieve.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $<

app_boinc64.o: app.c app.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $<

cl_app_boinc64.o: app.c app.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $(CL_DEFINE) $<

appcu_boinc64.o: appcu.cu appcu.h
	$(NVCC) -c -m64 -o $@ $(BOINC_NVCFLAGS) $<

#  COMMON64_OBJS
../clock64.o: ../clock.c ../clock.h
	$(CC64) -o $@ $(CFLAGS) $<

factor_proth64.o: factor_proth.c factor_proth.h
	$(CC64) -o $@ $(CFLAGS) $<

cuda_sleep_memcpy64.o: cuda_sleep_memcpy.cu cuda_sleep_memcpy.h
	$(NVCC) -c -m64 -o $@ $(NVCFLAGS) $<

appcl64.o: appcl.cpp appcl.h
	$(CL64) -o $@ $(CFLAGS) $<

# Generate the appcl.h file from the .h part and the appcl.cl part.
appcl.h: appcl_start.h appcl.cl cltoh.pl
	perl cltoh.pl appcl.cl | cat appcl_start.h - > appcl.h
