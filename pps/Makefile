# Build Linux apps by default, set PLATFORM=mac for Mac build
PLATFORM ?= linux
TPS_APPNAME = tpsieve
APPNAME = ppsieve
APPPATH = pps
COMMON_OBJS = ../clock.o factor_proth.o
OBJS = ../main.o ../putil.o ../sieve.o app.o $(COMMON_OBJS)
TPSOBJS = ../tpmain.o ../putil.o ../sieve.o tpapp.o $(COMMON_OBJS)
BOINC_TPSOBJS = ../do_boinc_init.o ../tpmain_boinc.o ../putil_boinc.o ../sieve_boinc.o tpapp_boinc.o $(COMMON_OBJS)
BOINC_OBJS = ../do_boinc_init.o ../main_boinc.o ../putil_boinc.o ../sieve_boinc.o app_boinc.o $(COMMON_OBJS)
BOINC_BINS = $(APPNAME)-boinc-x86-$(PLATFORM)
NON_BOINC_BINS = $(APPNAME)-x86-$(PLATFORM)
# BINS should also have $(APPNAME)-x86-windows.exe $(APPNAME)-boinc-x86-windows.exe, but those can't be made on Linux. (Yet.)
BINS = $(NON_BOINC_BINS) $(BOINC_BINS)
TPS_BINS = $(TPS_APPNAME)-x86-$(PLATFORM)
TPS_BOINC_BINS = $(TPS_APPNAME)-boinc-x86-$(PLATFORM)

ifeq ($(PLATFORM),linux)
  CC = gcc
  CC32 = $(CC) -c -m32 -march=i586 -mtune=core2
  NVCC = g++
  LD32 = $(CC) -m32
else ifeq ($(PLATFORM), mac)
  TARGET_OS = MACOSX_DEPLOYMENT_TARGET=10.5
  MACOS_MIN_VERSION = -mmacosx-version-min=10.5
  ISYSROOT = -isysroot /Developer/SDKs/MacOSX10.5.sdk
  CC = $(TARGET_OS) gcc
  CC32 = $(CC) -c -arch i386
  NVCC = $(TARGET_OS) g++
  LD32 = $(CC) -arch i386
endif

# Flags:
# Uncomment to use device emulation:
#DOEMU = -deviceemu -D_DEVICEEMU -g
COMMON_FLAGS = -O3 -DNDEBUG -D_REENTRANT -I. -I..
CFLAGS = -fomit-frame-pointer -s $(COMMON_FLAGS)
NVCFLAGS = $(DOEMU) $(COMMON_FLAGS)
LINK = -lm -lpthread # -lcuda

ifeq ($(PLATFORM),mac)
  COMMON_FLAGS += -D_GNU_SOURCE
  CFLAGS += $(MACOS_MIN_VERSION) $(ISYSROOT)
  NVCFLAGS += -arch=sm_11
  LINK += $(ISYSROOT) -lstdc++
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
  LDFLAGS = $(NVCFLAGS)
  BOINC_LDFLAGS = $(BOINC_NVCFLAGS)
  BOINC_LINK = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc `g++ -m32 -print-file-name=libstdc++.a` -DAPP_GRAPHICS
else ifeq($(PLATFORM),mac)
  BOINC_LINK = $(BOINC_LINK_LIBS) $(LINK) -lboinc_api -lboinc  -DAPP_GRAPHICS
endif

#nvcc $DOEMU --ptxas-options=-v -O3 $cleanvars -DUSE_BOINC -DNDEBUG -D_REENTRANT -m32 -I. -I.. -o $appname-boinc-x86-$(PLATFORM) $BOINC_LOAD_LIBS ../main.c ../sieve.c ../clock.c ../putil.c ../do_boinc_init.cpp cuda_sleep_memcpy.cu appcu.cu app.c factor_proth.c $link -lboinc_api -lboinc `g++ -print-file-name=libstdc++.a` -DAPP_GRAPHICS

# Make-able targets:
all: $(BINS) $(TPS_BINS) $(TPS_BOINC_BINS)

# Making zipfiles (requires binaries):
zips: $(APPNAME).zip $(TPS_APPNAME).zip $(APPNAME)-src.zip

$(TPS_APPNAME).zip: $(TPS_BINS) $(TPS_BOINC_BINS) gitstatus changes
	-rm $@ ; zip -9 $(TPS_APPNAME).zip $(TPS_APPNAME)-x86-$(PLATFORM) $(TPS_APPNAME)-boinc-x86-$(PLATFORM) $(TPS_APPNAME)-boinc-x86-windows.exe $(TPS_APPNAME)-x86-windows.exe README.txt CHANGES.txt $(TPS_APPNAME).sh $(TPS_APPNAME).bat *config.txt license.txt

$(APPNAME).zip: $(BINS) gitstatus changes
	-rm $@ ; zip -9 $(APPNAME).zip $(APPNAME)-x86-$(PLATFORM) $(APPNAME)-boinc-x86-$(PLATFORM) $(APPNAME)-boinc-x86-windows.exe $(APPNAME)-x86-windows.exe README.txt CHANGES.txt $(APPNAME).sh $(APPNAME).bat *config.txt license.txt

$(APPNAME)-src.zip: $(BINS) gitstatus changes
	-rm $@ ; cd .. ; zip -9 $(APPPATH)/$(APPNAME)-src.zip Makefile $(APPPATH)/Makefile *.txt *.[ch] *.cpp *.sln $(APPPATH)/*.[ch] $(APPPATH)/README.txt $(APPPATH)/CHANGES.txt $(APPPATH)/*config.txt $(APPPATH)/make-*.sh $(APPPATH)/$(APPNAME).sh $(APPPATH)/*.bat vc/*

gitstatus:
	git status | grep nothing > /dev/null || git status

changes:
	git log | grep "^\(    \|Date:\)" > CHANGES.txt

# Binaries:
tps: $(TPS_BINS) $(TPS_BOINC_BINS)

non-boinc: $(NON_BOINC_BINS) $(TPS_BINS)

boinc: $(BOINC_BINS) $(TPS_BOINC_BINS)

clean:
	-rm *.o ../*.o $(BINS) $(TPS_BINS) $(TPS_BOINC_BINS) $(TPS_APPNAME).zip $(APPNAME).zip $(APPNAME)-src.zip

$(APPNAME)-x86-$(PLATFORM): $(OBJS)
	$(LD32) -o $@ $(LDFLAGS) $^ $(LINK)

$(APPNAME)-boinc-x86-$(PLATFORM): $(BOINC_OBJS)
	$(LD32) -o $@ $(BOINC_LDFLAGS) $^ $(BOINC_LINK)

$(TPS_APPNAME)-x86-$(PLATFORM): $(TPSOBJS)
	$(LD32) -o $@ $(LDFLAGS) $^ $(LINK) -DSEARCH_TWIN

$(TPS_APPNAME)-boinc-x86-$(PLATFORM): $(BOINC_TPSOBJS)
	$(LD32) -o $@ $(BOINC_LDFLAGS) $^ $(BOINC_LINK) -DSEARCH_TWIN

# OBJS
../main.o: ../main.c ../main.h ../clock.h ../sieve.h ../putil.h app.h
	$(CC32) -o $@ $(CFLAGS) $<

../putil.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

../sieve.o: ../sieve.c ../sieve.h ../putil.h
	$(CC32) -o $@ $(CFLAGS) $<

app.o: app.c app.h ../putil.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(CFLAGS) $<

../tpmain.o: ../main.c ../main.h ../clock.h ../sieve.h ../putil.h app.h
	$(CC32) -o $@ $(CFLAGS) $< -DSEARCH_TWIN

tpapp.o: app.c app.h ../putil.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(CFLAGS) $< -DSEARCH_TWIN


# BOINC_OBJS
../do_boinc_init.o: ../do_boinc_init.cpp ../do_boinc_init.h
	$(NVCC) -c -m32 -o $@ $(BOINC_NVCFLAGS) $<

../main_boinc.o: ../main.c ../main.h ../clock.h ../sieve.h ../putil.h app.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

../putil_boinc.o: ../putil.c ../putil.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

../sieve_boinc.o: ../sieve.c ../sieve.h ../putil.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

app_boinc.o: app.c app.h ../putil.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $<

../tpmain_boinc.o: ../main.c ../main.h ../clock.h ../sieve.h ../putil.h app.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $< -DSEARCH_TWIN

tpapp_boinc.o: app.c app.h ../putil.h ../clock.h factor_proth.h ../main.h
	$(CC32) -o $@ $(BOINC_CFLAGS) $< -DSEARCH_TWIN

#  COMMON_OBJS
../clock.o: ../clock.c ../clock.h
	$(CC32) -o $@ $(CFLAGS) $<

factor_proth.o: factor_proth.c factor_proth.h
	$(CC32) -o $@ $(CFLAGS) $<
