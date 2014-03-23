# Build Linux apps by default, set PLATFORM=mac for Mac build
PLATFORM ?= linux
CL_APPNAME = ppsieve-cl
CL_TPS_APPNAME = tpsieve-cl
APPPATH = pps
COMMON_OBJS = ../clock.o factor_proth.o

# Non-BOINC 32-bit objects.
OBJS = ../main.o ../putil.o ../sieve.o $(COMMON_OBJS)
CL_OBJS = $(OBJS) cl_app.o appcl.o
CL_TPSOBJS = ../tpmain.o ../putil.o ../sieve.o cl_tpapp.o tpappcl.o $(COMMON_OBJS)

# BOINC 32-bit objects.
BOINC_OBJS = ../do_boinc_init.o ../main_boinc.o ../putil_boinc.o ../sieve_boinc.o $(COMMON_OBJS)
CL_BOINC_OBJS = $(BOINC_OBJS) cl_app_boinc.o appcl.o
BOINC_TPSOBJS = ../do_boinc_init.o ../tpmain_boinc.o ../putil_boinc.o ../sieve_boinc.o tpapp_boinc.o $(COMMON_OBJS)
CL_BOINC_TPSOBJS = ../do_boinc_init.o ../tpmain_boinc.o ../putil_boinc.o ../sieve_boinc.o cl_tpapp_boinc.o tpappcl.o $(COMMON_OBJS)

# 64-bit objects.
COMMON64_OBJS = ../clock64.o factor_proth64.o

# non-BOINC 64-bit objects.
OBJS64 = ../main64.o ../putil64.o ../sieve64.o $(COMMON64_OBJS)
CL_OBJS64 = $(OBJS64) appcl64.o cl_app64.o
TPSOBJS64 = ../tpmain64.o ../putil64.o ../sieve64.o tpapp64.o $(COMMON64_OBJS)
CL_TPSOBJS64 = ../tpmain64.o ../putil64.o ../sieve64.o cl_tpapp64.o tpappcl64.o $(COMMON64_OBJS)

# BOINC 64-bit objects.
BOINC_OBJS64 = ../do_boinc_init64.o ../main_boinc64.o ../putil_boinc64.o ../sieve_boinc64.o $(COMMON64_OBJS)
CL_BOINC_OBJS64 = $(BOINC_OBJS64) cl_app_boinc64.o appcl64.o
BOINC_TPSOBJS64 = ../do_boinc_init64.o ../tpmain_boinc64.o ../putil_boinc64.o ../sieve_boinc64.o tpapp_boinc64.o $(COMMON64_OBJS)
CL_BOINC_TPSOBJS64 = ../do_boinc_init64.o ../tpmain_boinc64.o ../putil_boinc64.o ../sieve_boinc64.o cl_tpapp_boinc64.o tpappcl64.o $(COMMON64_OBJS)

# Bins
BOINC_BINS = $(APPNAME)-boinc-x86-$(PLATFORM) $(APPNAME)-boinc-x86_64-$(PLATFORM)
NON_BOINC_BINS = $(APPNAME)-x86-$(PLATFORM) $(APPNAME)-x86_64-$(PLATFORM)
CL_BOINC_BINS = $(CL_APPNAME)-boinc-x86-$(PLATFORM) $(CL_APPNAME)-boinc-x86_64-$(PLATFORM) $(CL_TPS_APPNAME)-boinc-x86-$(PLATFORM) $(CL_TPS_APPNAME)-boinc-x86_64-$(PLATFORM)
CL_NON_BOINC_BINS = $(CL_APPNAME)-x86-$(PLATFORM) $(CL_APPNAME)-x86_64-$(PLATFORM) $(CL_TPS_APPNAME)-x86-$(PLATFORM) $(CL_TPS_APPNAME)-x86_64-$(PLATFORM)
# BINS should also have $(APPNAME)-x86-windows.exe $(APPNAME)-boinc-x86-windows.exe, but those can't be made on Linux. (Yet.)
BINS = $(NON_BOINC_BINS) $(BOINC_BINS)
CL_BINS = $(CL_NON_BOINC_BINS) $(CL_BOINC_BINS)
TPS_BINS = $(TPS_APPNAME)-x86-$(PLATFORM) $(TPS_APPNAME)-x86_64-$(PLATFORM) 
TPS_BOINC_BINS = $(TPS_APPNAME)-boinc-x86-$(PLATFORM) $(TPS_APPNAME)-boinc-x86_64-$(PLATFORM) 

ifeq ($(PLATFORM),linux)
  CC = gcc
  CC32 = $(CC) -c -m32 -march=i586 -mtune=core2
  CC64 = $(CC) -c -m64 -march=k8 -mno-3dnow -mtune=core2
  NVCC = g++
  CL_DEFINE = -DUSE_OPENCL
  CL_INCLUDE = -I $(ATISTREAMSDKROOT)/include
  CL32 = $(CXX) -c -m32 -march=i586 -mtune=core2 $(CL_INCLUDE)
  CL64 = $(CXX) -c -m64 -march=k8 -mno-3dnow -mtune=core2 $(CL_INCLUDE)
else ifeq ($(PLATFORM),mac)
  TARGET_OS = MACOSX_DEPLOYMENT_TARGET=10.6
  MACOS_MIN_VERSION = -mmacosx-version-min=10.6
  ISYSROOT = -isysroot /Developer/SDKs/MacOSX10.6.sdk
  CC = $(TARGET_OS) gcc
  CC32 = $(CC) -c -arch i386
  CC64 = $(CC) -c -arch x86_64
  NVCC = $(TARGET_OS) nvcc
  LD32 = $(CC) -arch i386
  LD64 = $(CC) -arch x86_64
  CL_DEFINE = -DUSE_OPENCL
  CL_INCLUDE = -framework OpenCL
  CL32 = $(CXX) -c -arch i386 $(CL_INCLUDE)
  CL64 = $(CXX) -c -arch x86_64 $(CL_INCLUDE)
endif
 

# Flags:
# Uncomment to use device emulation:
#DOEMU = -deviceemu -D_DEVICEEMU -g
COMMON_FLAGS = -O3 -DNDEBUG -D_REENTRANT -I. -I..
# Uncomment -D_DEVICEEMU below to use ATI OpenCL device emulation.
CFLAGS = -Wall -fomit-frame-pointer -s $(COMMON_FLAGS) #-D_DEVICEEMU
NVCFLAGS = --ptxas-options=-v $(DOEMU) $(COMMON_FLAGS)
LINK = -lm -lpthread # -lcuda
CL_LINK = -lOpenCL -L$(ATISTREAMSDKROOT)/lib/x86
CL_LINK64 = -lOpenCL -L$(ATISTREAMSDKROOT)/lib/x86_64

ifeq ($(PLATFORM),mac)
  COMMON_FLAGS += -D_GNU_SOURCE
  CFLAGS += $(MACOS_MIN_VERSION) $(ISYSROOT)
  NVCFLAGS += -arch=sm_11
  LINK += $(ISYSROOT) -lstdc++
  CL_LINK = -framework OpenCL
  CL_LINK64 = $(CL_LINK)
endif

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

ifeq ($(PLATFORM),linux)
  BOINC_LINK64 = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc `g++ -m64 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
  BOINC_LINK = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc `g++ -m32 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
else ifeq ($(PLATFORM),mac)
  BOINC_LINK64 = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc -DAPP_GRAPHICS
  BOINC_LINK = $(BOINC_LINK64)
endif


# Make-able targets:
all: $(BINS) cl tps

cl: appcl.h $(CL_BINS)

# Making zipfiles (requires binaries):
zips: $(APPNAME).zip $(TPS_APPNAME).zip $(APPNAME)-src.zip

$(TPS_APPNAME).zip: $(TPS_BINS) $(TPS_BOINC_BINS) gitstatus changes
	-rm $@ ; zip -9 $(TPS_APPNAME).zip $(TPS_APPNAME)-x86-$(PLATFORM) $(TPS_APPNAME)-boinc-x86-$(PLATFORM) $(TPS_APPNAME)-x86_64-$(PLATFORM) $(TPS_APPNAME)-boinc-x86_64-$(PLATFORM) $(TPS_APPNAME)-boinc-x86-windows.exe $(TPS_APPNAME)-x86-windows.exe $(TPS_APPNAME)-boinc-x86_64-windows.exe $(TPS_APPNAME)-x86_64-windows.exe README.txt CHANGES.txt $(TPS_APPNAME).sh $(TPS_APPNAME).bat *config.txt license.txt

$(APPNAME).zip: $(BINS) gitstatus changes
	-rm $@ ; zip -9 $(APPNAME).zip $(APPNAME)-x86-$(PLATFORM) $(APPNAME)-boinc-x86-$(PLATFORM) $(APPNAME)-x86_64-$(PLATFORM) $(APPNAME)-boinc-x86_64-$(PLATFORM) $(APPNAME)-boinc-x86-windows.exe $(APPNAME)-x86-windows.exe $(APPNAME)-boinc-x86_64-windows.exe $(APPNAME)-x86_64-windows.exe README.txt CHANGES.txt $(APPNAME).sh $(APPNAME).bat *config.txt license.txt

$(APPNAME)-src.zip: $(BINS) gitstatus changes
	-rm $@ ; cd .. ; zip -9 $(APPPATH)/$(APPNAME)-src.zip Makefile $(APPPATH)/Makefile *.txt *.[ch] *.cpp *.sln $(APPPATH)/*.[ch] $(APPPATH)/README.txt $(APPPATH)/CHANGES.txt $(APPPATH)/*config.txt $(APPPATH)/make-*.sh $(APPPATH)/$(APPNAME).sh $(APPPATH)/*.bat vc/*

gitstatus:
	git status | grep nothing > /dev/null || git status

changes:
	git log | grep "^\(    \|Date:\)" > CHANGES.txt

# Binaries:
tps: $(TPS_BINS) $(TPS_BOINC_BINS)
	
# $(CL_TPS_BINS) $(CL_TPS_BOINC_BINS)

non-boinc: $(NON_BOINC_BINS) $(TPS_BINS)

boinc: $(BOINC_BINS) $(TPS_BOINC_BINS)

clean:
	-rm *.o ../*.o appcl.h $(BINS) $(CL_BINS) $(TPS_BINS) $(TPS_BOINC_BINS) $(TPS_APPNAME).zip $(APPNAME).zip $(APPNAME)-src.zip

$(APPNAME)-x86-$(PLATFORM): $(CUDA_OBJS)
	$(NVCC) -m32 -o $@ $(NVCFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86-$(PLATFORM): $(CUDA_BOINC_OBJS)
	$(NVCC) -m32 -o $@ $(BOINC_NVCFLAGS) $^ $(BOINC_LINK)

$(TPS_APPNAME)-x86-$(PLATFORM): $(TPSOBJS)
	$(NVCC) -m32 -o $@ $(NVCFLAGS) $^ $(LINK) -DSEARCH_TWIN

$(TPS_APPNAME)-boinc-x86-$(PLATFORM): $(BOINC_TPSOBJS)
	$(NVCC) -m32 -o $@ $(BOINC_NVCFLAGS) $^ $(BOINC_LINK) -DSEARCH_TWIN

$(APPNAME)-x86_64-$(PLATFORM): $(CUDA_OBJS64)
	$(NVCC) -m64 -o $@ $(NVCFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86_64-$(PLATFORM): $(CUDA_BOINC_OBJS64)
	$(NVCC) -m64 -o $@ $(BOINC_NVCFLAGS) $^ $(BOINC_LINK64)

$(TPS_APPNAME)-x86_64-$(PLATFORM): $(TPSOBJS64)
	$(NVCC) -m64 -o $@ $(NVCFLAGS) $^ $(LINK) -DSEARCH_TWIN

$(TPS_APPNAME)-boinc-x86_64-$(PLATFORM): $(BOINC_TPSOBJS64)
	$(NVCC) -m64 -o $@ $(BOINC_NVCFLAGS) $^ $(BOINC64_LINK) -DSEARCH_TWIN

$(CL_APPNAME)-x86-$(PLATFORM): $(CL_OBJS)
	$(CXX) -m32 -o $@ $(CFLAGS) $^ $(CL_LINK) $(LINK)

$(CL_APPNAME)-boinc-x86-$(PLATFORM): $(CL_BOINC_OBJS)
	$(CXX) -m32 -o $@ $(BOINC_CFLAGS) $^ $(CL_LINK) $(BOINC_LINK)

$(CL_TPS_APPNAME)-x86-$(PLATFORM): $(CL_TPSOBJS)
	$(CXX) -m32 -o $@ $(CFLAGS) $^ $(CL_LINK) $(LINK) -DSEARCH_TWIN

$(CL_TPS_APPNAME)-boinc-x86-$(PLATFORM): $(CL_BOINC_TPSOBJS)
	$(CXX) -m32 -o $@ $(BOINC_CFLAGS) $^ $(CL_LINK) $(BOINC_LINK) -DSEARCH_TWIN

$(CL_APPNAME)-x86_64-$(PLATFORM): $(CL_OBJS64)
	$(CXX) -m64 -o $@ $(CFLAGS) $^ $(CL_LINK64) $(LINK)

$(CL_APPNAME)-boinc-x86_64-$(PLATFORM): $(CL_BOINC_OBJS64)
	$(CXX) -m64 -o $@ $(BOINC_CFLAGS) $^ $(CL_LINK64) $(BOINC_LINK64)

$(CL_TPS_APPNAME)-x86_64-$(PLATFORM): $(CL_TPSOBJS64)
	$(CXX) -m64 -o $@ $(CFLAGS) $^ $(CL_LINK64) $(LINK) -DSEARCH_TWIN

$(CL_TPS_APPNAME)-boinc-x86_64-$(PLATFORM): $(CL_BOINC_TPSOBJS64)
	$(CXX) -m64 -o $@ $(BOINC_CFLAGS) $^ $(CL_LINK64) $(BOINC_LINK64) -DSEARCH_TWIN

# OBJS
../putil.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

../sieve.o: ../sieve.c ../sieve.h ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

cl_app.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(CFLAGS) $(CL_DEFINE) $<

cl_tpapp.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(CFLAGS) $(CL_DEFINE) $< -DSEARCH_TWIN


# BOINC_OBJS
../do_boinc_init.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(CXX) -c -m32 -o $@ $(BOINC_CFLAGS) $<

../putil_boinc.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

../sieve_boinc.o: ../sieve.c ../sieve.h ../putil.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

cl_app_boinc.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $(CL_DEFINE) $<

cl_tpapp_boinc.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $(CL_DEFINE) $< -DSEARCH_TWIN

#  COMMON_OBJS
../clock.o: ../clock.c ../clock.h
	$(CC32) -o $@ $(CFLAGS) $<

factor_proth.o: factor_proth.c factor_proth.h
	$(CC32) -o $@ $(CFLAGS) $<

appcl.o: appcl.cpp appcl.h
	$(CL32) -o $@ $(CFLAGS) $<

tpappcl.o: appcl.cpp appcl.h
	$(CL32) -o $@ $(CFLAGS) $< -DSEARCH_TWIN

# OBJS64
../putil64.o: ../putil.c ../putil.h
	$(CC64) -o $@ $(CFLAGS) $<

../sieve64.o: ../sieve.c ../sieve.h ../putil.h
	$(CC64) -o $@ $(CFLAGS) $<

cl_app64.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC64) -o $@ $(CFLAGS) $(CL_DEFINE) $<

cl_tpapp64.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC64) -o $@ $(CFLAGS) $(CL_DEFINE) $< -DSEARCH_TWIN


# BOINC_OBJS64
../do_boinc_init64.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(CXX) -c -m64 -o $@ $(BOINC_CFLAGS) $<

../putil_boinc64.o: ../putil.c ../putil.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $<

../sieve_boinc64.o: ../sieve.c ../sieve.h ../putil.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $<

cl_app_boinc64.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $(CL_DEFINE) $<

cl_tpapp_boinc64.o: app.c app.h ../putil.h appcl.h ../clock.h factor_proth.h ../main.h
	$(CC64) -o $@ $(BOINC_CFLAGS) $(CL_DEFINE) $< -DSEARCH_TWIN


#  COMMON64_OBJS
../clock64.o: ../clock.c ../clock.h
	$(CC64) -o $@ $(CFLAGS) $<

factor_proth64.o: factor_proth.c factor_proth.h
	$(CC64) -o $@ $(CFLAGS) $<

appcl64.o: appcl.cpp appcl.h
	$(CL64) -o $@ $(CFLAGS) $<

tpappcl64.o: appcl.cpp appcl.h
	$(CL64) -o $@ $(CFLAGS) $< -DSEARCH_TWIN

# Generate the appcl.h file from the .h part and the appcl.cl part.
appcl.h: appcl_start.h appcl.cl cltoh.pl
	perl cltoh.pl appcl.cl | cat appcl_start.h - > appcl.h
