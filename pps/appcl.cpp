/* ex: set softtabstop=2 shiftwidth=2 expandtab: */
/* app.cu -- (C) Ken Brazier, August 2010.

   Proth Prime Search sieve OpenCL calling portion (for many K and many N).

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <utility>
#define __NO_STD_VECTOR
#include <CL/cl.hpp>
#include "main.h"
#include "putil.h"
#include "app.h"
#include "appcl.h"

#include <unistd.h>

#define INLINE static inline
/*
#ifndef BITSATATIME
#define BITSATATIME 4
#endif
#define BITSMASK ((1<<BITSATATIME)-1)*/
/*
#if(BITSATATIME == 3)
#define SHIFT_CAST unsigned int
#elif(BITSATATIME == 4)
#define SHIFT_CAST uint64_t
#else
#error "Invalid BITSATATIME."
#endif
 */
// Extern vars in appcl.hpp:
unsigned int ld_nstep;
int ld_bbits;
uint64_t ld_r0;

// Device constants
//__constant__ unsigned int d_bitsatatime;
//__constant__ unsigned int d_len;//=(1<<bitsatatime); 
//__constant__ unsigned int d_halflen;//=(1<<bitsatatime)/2; 
/*__constant__ uint64_t d_kmin;
__constant__ uint64_t d_kmax;
__constant__ unsigned int d_nmin;
__constant__ unsigned int d_nmax;
__constant__ unsigned int d_nstep;
__constant__ unsigned int d_kernel_nstep;
__constant__ unsigned int d_search_proth;

__constant__ int d_bbits;
__constant__ unsigned int d_mont_nstep;
__constant__ uint64_t d_r0;*/
// Device arrays
// Output:
cl_mem d_factor_found;
// Input (read-only):
cl_mem d_P;
// Local storage, between kernels:
cl_mem d_Ps, d_K;
//unsigned int *d_N;

// Timing variables:
//const int setup_ps_overlap = 5000;
//const int check_ns_overlap = 50000;

static unsigned int ld_kernel_nstep;
static size_t global_cthread_count[1];
static cl_event comp_done_event;

//static bool blocking_sync_ok=true;
const char *source;
//globals for OpenCL
/* problem size for a 2D matrix. */
// Note: we will handle the problem as a 1D matrix.
cl_uint width;
cl_uint height;

cl_context          context;
cl_device_id        *devices;
cl_command_queue    commandQueue;

cl_program program;

/* This program uses two kernels */
cl_kernel  start_ns_kernel;
cl_kernel  check_more_ns_kernel;

// find the log base 2 of a number.  Need not be fast; only done once.
int lg2(uint64_t v) {
  int r = 0; // r will be lg(v)

  while (v >>= 1) // unroll for more speed...
  {
    r++;
  }
  return r;
}

void cuda_finalize(void) {
  cl_int status;

  status = clReleaseKernel(start_ns_kernel);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseKernel (start_ns_kernel)\n");
  }

  status = clReleaseKernel(check_more_ns_kernel);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseKernel (check_more_ns_kernel)\n");
  }

  status = clReleaseProgram(program);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseProgram\n");
  }

  status = clReleaseMemObject(d_P);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseMemObject (inputBuffer)\n");
  }

  status = clReleaseMemObject(d_Ps);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseMemObject (inputBuffer)\n");
  }

  status = clReleaseMemObject(d_K);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseMemObject (inputBuffer)\n");
  }

  status = clReleaseMemObject(d_factor_found);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseMemObject (inputBuffer)\n");
  }

  status = clReleaseCommandQueue(commandQueue);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseCommandQueue\n");
  }

  status = clReleaseContext(context);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: In clReleaseContext\n");
  }
}

static const char* printable_cl_error(cl_int err) {
  switch (err) {
    case CL_SUCCESS:                            return "Success!";
    case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
    case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
    case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
    case CL_OUT_OF_RESOURCES:                   return "Out of resources";
    case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
    case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
    case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
    case CL_MAP_FAILURE:                        return "Map failure";
    case CL_INVALID_VALUE:                      return "Invalid value";
    case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
    case CL_INVALID_PLATFORM:                   return "Invalid platform";
    case CL_INVALID_DEVICE:                     return "Invalid device";
    case CL_INVALID_CONTEXT:                    return "Invalid context";
    case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
    case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
    case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
    case CL_INVALID_SAMPLER:                    return "Invalid sampler";
    case CL_INVALID_BINARY:                     return "Invalid binary";
    case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
    case CL_INVALID_PROGRAM:                    return "Invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
    case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
    case CL_INVALID_KERNEL:                     return "Invalid kernel";
    case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
    case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
    case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
    case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
    case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
    case CL_INVALID_EVENT:                      return "Invalid event";
    case CL_INVALID_OPERATION:                  return "Invalid operation";
    case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
    case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
    case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
    default: return "Unknown";
  }
}

static void checkCUDAErr(cl_int err, const char* msg) {
  if(err != CL_SUCCESS) {
    fprintf(stderr, "Cuda error: %s: %s\n", msg, printable_cl_error(err));
    cuda_finalize();
    bexit(EXIT_FAILURE);
  }
}


/*
 * Converts the contents of a file into a string
 * To read the .cl file.
 */
/*
char * convert_to_string(const char *fileName) {
    int count=0;
    char *s;
    char c;
    int i=0;

    // look for "atiopencl_kernels.cl" in "boinc/samples/atiopencl/debug" or
    // in "boinc/samples/atiopencl/release". Note that "atiopencl_kernels.cl"
    // is automatically copied to these directories along the building process.
    FILE *infile=fopen(fileName,"r");
    if (!infile) { //not found. This typically happens on Linux or Mac.
        //look for "atiopencl_kernels.cl" in "boinc/sample/atiopencl/" instead.
        infile = fopen(KERNELS_FILEPATH,"r");
        if (!infile) {
            bmsg("appcl.cl file open Error!");
            bexit(ERR_FOPEN);
        }
    }
    fseek(infile,0,SEEK_SET);
    while (fgetc(infile)!=EOF) count++;
    s=(char *) malloc(sizeof(char)*(count+1)); //add 1 for string terminator.
    fseek(infile,0,SEEK_SET);   
    while ((c=fgetc(infile))!=EOF) {
        s[i++]=c;
    }
    s[i]='\0';
    return s;
}*/

// Just return the string from the .h file instead.
const char *convert_to_string(const char *fileName) {
  return appcl;
}

/*
 * \brief OpenCL related initialization
 *        Create Context, Device list, Command Queue
 *        Load CL file, compile, link CL source
 *                Build program and kernel objects
 */

 // Note: OpenCL memory buffer objects will be created in calling
 //       function before kernel calls are made.
static int initialize_cl(int deviceno, unsigned int *cthread_count) {
  cl_int status = 0;
  size_t deviceListSize;

  //localThreads[0]  = LOCAL_WORK_SIZE;
  //globalThreads[0] = GLOBAL_WORK_SIZE;

  /*
   * Have a look at the available platforms and pick either
   * the AMD one if available or a reasonable default.
   */

  cl_uint numPlatforms;
  cl_platform_id platform = NULL;
  status = clGetPlatformIDs(0, NULL, &numPlatforms);
  if(status != CL_SUCCESS) {
    fprintf(stderr, "Error: Getting Platforms. (clGetPlatformsIDs)\n");
    return 1;
  }

  if (numPlatforms > 0) {
    cl_platform_id* platforms = (cl_platform_id *)
      malloc(sizeof(cl_platform_id)*numPlatforms);
    status = clGetPlatformIDs(numPlatforms, platforms, NULL);
    if (status != CL_SUCCESS) {
      fprintf(stderr, "Error: Getting Platform Ids. (clGetPlatformsIDs)\n");
      return 1;
    }
    for (unsigned int i=0; i < numPlatforms; ++i) {
      char pbuff[100];
      status = clGetPlatformInfo(platforms[i],
          CL_PLATFORM_VENDOR,
          sizeof(pbuff),
          pbuff,
          NULL);
      if (status != CL_SUCCESS) {
        fprintf(stderr, "Error: Getting Platform Info.(clGetPlatformInfo)\n");
        return 1;
      }
      platform = platforms[i];
      if (!strcmp(pbuff, "Advanced Micro Devices, Inc.")) {
        break;
      }
    }
    delete platforms;
  }

  if(NULL == platform) {
    fprintf(stderr, "NULL platform found so Exiting Application.\n");
    return 1;
  }

  /*
   * If we could find our platform, use it. Otherwise use just available platform.
   */
  cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM,
    (cl_context_properties)platform,
    0
  };

  /////////////////////////////////////////////////////////////////
  // Create an OpenCL context
  /////////////////////////////////////////////////////////////////
  context = clCreateContextFromType(cps, CL_DEVICE_TYPE_ALL, NULL, NULL, &status);
  if (status != CL_SUCCESS) { 
    fprintf(stderr, "Error: Creating Context. (clCreateContextFromType)\n");
    return 1;
  }

  /* First, get the size of device list data */
  status = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize);
  if (status != CL_SUCCESS) { 
    fprintf(stderr, "Error: Getting Context Info (device list size, clGetContextInfo)\n");
    return 1;
  }
  if(deviceListSize <= (size_t)deviceno) {
    fprintf(stderr, "Error: Device %d not found.\n", deviceno);
    return 1;
  }

  /////////////////////////////////////////////////////////////////
  // Detect OpenCL devices
  /////////////////////////////////////////////////////////////////
  devices = (cl_device_id *)malloc(deviceListSize);
  if (devices == 0) {
    fprintf(stderr, "Error: No devices found.\n");
    return 1;
  }

  /* Now, get the device list data */
  status = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceListSize, devices, NULL);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: Getting Context Info (device list, clGetContextInfo)\n");
    return 1;
  }

  /////////////////////////////////////////////////////////////////
  // Create an OpenCL command queue
  /////////////////////////////////////////////////////////////////
  commandQueue = clCreateCommandQueue(context, devices[deviceno], 0, &status);
  if(status != CL_SUCCESS) {
    fprintf(stderr, "Creating Command Queue. (clCreateCommandQueue)\n");
    return 1;
  }

  /////////////////////////////////////////////////////////////////
  // Load CL file, build CL program object, create CL kernel object
  /////////////////////////////////////////////////////////////////
  source = convert_to_string(KERNELS_FILENAME);
  size_t sourceSize[]    = { strlen(source) };
  program = clCreateProgramWithSource(context, 1, &source, sourceSize, &status);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)\n");
    return 1;
  }

  /* create a cl program executable for all the devices specified */
  status = clBuildProgram(program, 1, &devices[deviceno], "-g", NULL, NULL);
  if (status != CL_SUCCESS)  {
    fprintf(stderr, "Error: Building Program (clBuildProgram)\n");
    return 1;
  }

  /* get a kernel object handle for a kernel with the given name */
  start_ns_kernel = clCreateKernel(program, "start_ns", &status);
  if (status != CL_SUCCESS) { 
    fprintf(stderr, "Error: clCreateKernel (start_ns)\n");
    return 1;
  }

  check_more_ns_kernel = clCreateKernel(program, "check_more_ns", &status);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "Error: clCreateKernel (check_more_ns)\n");
    return 1;
  }

  // Knowing the device, we need to calculate how many threads to give it.
  cl_uint compute_units;
  clGetDeviceInfo(devices[deviceno],
      CL_DEVICE_MAX_COMPUTE_UNITS,
      sizeof(cl_uint),
      &compute_units,
      NULL);

  fprintf(stderr, "%sDetected %d multiprocessors (%d SPUs) on device %d.\n",
      bmprefix(), compute_units, compute_units*5, deviceno);
  // 7 wavefronts per SIMD
  // Double this if using ulong2.
  *cthread_count = compute_units * (7 * BLOCKSIZE);

  // N's to search each time a kernel is run:
  ld_kernel_nstep = ITERATIONS_PER_KERNEL;
  // Adjust for differing block sizes.
  ld_kernel_nstep *= 384;
  ld_kernel_nstep /= (*cthread_count/compute_units);
  ld_kernel_nstep *= ld_nstep;

  return 0;
}

/* This function is called once before any threads are started.
 */
unsigned int cuda_app_init(int gpuno, unsigned int cthread_count)
{
  unsigned int i;
  cl_int status = 0;
  //unsigned int ld_bitsatatime = 0;
  //unsigned int ld_halflen=(1<<bitsatatime)/2; 
  //unsigned int ld_bitsmask;
  //unsigned int ld_bpernstep;

  if(initialize_cl(gpuno, &cthread_count)) {
    bexit(ERR_NOT_IMPLEMENTED);
  }
  /* Assume N >= 2^32. */
  if(pmin <= ((uint64_t)1)<<32) {
    bmsg("Error: PMin is too small, <= 2^32!\n");
    bexit(1);
  }

  /*
  if(gpuprop.totalGlobalMem < cthread_count*13) {
    fprintf(stderr, "%sInsufficient GPU memory: %u bytes.\n", bmprefix(), (unsigned int)(gpuprop.totalGlobalMem));
#ifdef USE_BOINC
    bexit(1);
#else
    return 0;
#endif
  }
  // Calculate ld_bitsatatime given memory constraints, and possibly nmin-nmax via nstep vs. 2^ld_bitsatatime
  // Things change if nmax-nmin < 1000000 or so, but for now let's go with a constant maximum of ld_bitsatatime<=13.
  i = gpuprop.totalGlobalMem/sizeof(uint64_t); // Total number of 64-bit numbers that can be stored.
  //ld_bitsatatime = BITSATATIME;
  //ld_bitsmask = BITSMASK+1;
  */

  // Allocate device arrays:
  //printf("cthread_count = %u\n", cthread_count);
  d_P = clCreateBuffer(context, CL_MEM_READ_ONLY, cthread_count*sizeof(cl_ulong), NULL, &status);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "%sInsufficient available memory on GPU %d.\n", bmprefix(), gpuno);
    checkCUDAErr(status, "Creating buffer d_P");
    bexit(ERR_INSUFFICIENT_RESOURCE);
  }
  CL_SET_BUF_ARG(start_ns_kernel, d_P);
  CL_SET_BUF_ARG(check_more_ns_kernel, d_P);

  d_Ps = clCreateBuffer(context, CL_MEM_READ_WRITE, cthread_count*sizeof(cl_ulong), NULL, &status);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "%sInsufficient available memory on GPU %d.\n", bmprefix(), gpuno);
    checkCUDAErr(status, "Creating buffer d_Ps");
    bexit(ERR_INSUFFICIENT_RESOURCE);
  }
  CL_SET_BUF_ARG(start_ns_kernel, d_Ps);
  CL_SET_BUF_ARG(check_more_ns_kernel, d_Ps);

  d_K = clCreateBuffer(context, CL_MEM_READ_WRITE, cthread_count*sizeof(cl_ulong), NULL, &status);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "%sInsufficient available memory on GPU %d.\n", bmprefix(), gpuno);
    checkCUDAErr(status, "Creating buffer d_K");
    bexit(ERR_INSUFFICIENT_RESOURCE);
  }
  CL_SET_BUF_ARG(start_ns_kernel, d_K);
  CL_SET_BUF_ARG(check_more_ns_kernel, d_K);

  d_factor_found = clCreateBuffer(context, CL_MEM_READ_WRITE, cthread_count*sizeof(cl_uchar), NULL, &status);
  if (status != CL_SUCCESS) {
    fprintf(stderr, "%sInsufficient available memory on GPU %d.\n", bmprefix(), gpuno);
    checkCUDAErr(status, "Creating buffer d_factor_found");
    bexit(ERR_INSUFFICIENT_RESOURCE);
  }
  CL_SET_BUF_ARG(start_ns_kernel, d_factor_found);
  CL_SET_BUF_ARG(check_more_ns_kernel, d_factor_found);

  //ld_bitsmask--; // Finalize bitsmask

  if (ld_nstep > (nmax-nmin+1))
    ld_nstep = (nmax-nmin+1);

  //assert((1ul << (64-nstep)) < pmin);
  if((((uint64_t)1) << (64-ld_nstep)) > pmin) {
    bmsg("Error: pmin is not large enough (or nmax is close to nmin).\n");
    bexit(ERR_INVALID_PARAM);
  }
  // Set the constants.
  //CL_MEMCPY_TO_SYMBOL(d_bitsatatime, &ld_bitsatatime, sizeof(ld_bitsatatime));

  // Prepare constants:
  ld_bbits = lg2(nmin);
  //assert(d_r0 <= 32);
  if(ld_bbits < 6) {
    fprintf(stderr, "%sError: nmin too small at %d (must be at least 64).\n", bmprefix(), nmin);
    bexit(ERR_INVALID_PARAM);
  }
  // r = 2^-i * 2^64 (mod N), something that can be done in a uint64_t!
  // If i is large (and it should be at least >= 32), there's a very good chance no mod is needed!
  ld_r0 = ((uint64_t)1) << (64-(nmin >> (ld_bbits-5)));

  ld_bbits = ld_bbits-6;
  CL_MEMCPY_TO_SYMBOL(d_bbits, &ld_bbits, sizeof(ld_bbits));
  // d_mont_nstep is the montgomerized version of nstep.
  i = 64-ld_nstep;
  CL_MEMCPY_TO_SYMBOL(d_mont_nstep, &i, sizeof(i));
  CL_MEMCPY_TO_SYMBOL(d_r0, &ld_r0, sizeof(ld_r0));
  CL_MEMCPY_TO_SYMBOL(d_kernel_nstep, &ld_kernel_nstep, sizeof(ld_kernel_nstep));
  CL_MEMCPY_TO_SYMBOL(d_kmax, &kmax, sizeof(kmax));
  CL_MEMCPY_TO_SYMBOL(d_kmin, &kmin, sizeof(kmin));
  CL_MEMCPY_TO_SYMBOL(d_nmin, &nmin, sizeof(nmin));
  CL_MEMCPY_TO_SYMBOL(d_nmax, &nmax, sizeof(nmax));
  CL_MEMCPY_TO_SYMBOL(d_nstep, &ld_nstep, sizeof(ld_nstep));
  i = (search_proth == 1)?1:0;	// search_proth is 1 or -1, not 0.
  CL_MEMCPY_TO_SYMBOL(d_search_proth, &i, sizeof(i));

  global_cthread_count[0] = cthread_count;

  return cthread_count;
}


// *** Host Kernel-calling functions ***

// Pass the arguments to the CUDA device, run the code, and get the results.
void check_ns(const uint64_t *P, const unsigned int cthread_count) {
  //const unsigned int cblockcount = cthread_count/BLOCKSIZE;
  unsigned int n;
  // timing variables:

  // Pass P.
  //printf("Starting tests on P's starting with %lu\n", P[0]);
  //cudaMemcpy(d_P, P, cthread_count*sizeof(uint64_t), cudaMemcpyHostToDevice);
  checkCUDAErr(clEnqueueWriteBuffer(commandQueue,
        d_P,
        CL_TRUE,
        0,
        cthread_count*sizeof(cl_ulong),
        P,
        0,
        NULL,
        &comp_done_event)
      ,"cudaMemcpy");

  /* Wait for the write buffer to finish execution */
  checkCUDAErr(clWaitForEvents(1, &comp_done_event), "Waiting for write buffer call to finish. (clWaitForEvents)");
  // Let the loop release the event.
  //checkCUDAErr(clReleaseEvent(dev_write_event), "Release event object. (clReleaseEvent)");

#ifndef NDEBUG
  bmsg("Setup successful...\n");
#endif
  //printf("cthread_count now = %u\n", cthread_count);
  checkCUDAErr(clEnqueueNDRangeKernel(commandQueue,
        start_ns_kernel,
        1,
        NULL,
        global_cthread_count,
        NULL,
        0,
        NULL,
        NULL),
      "kernel invocation");
#ifndef NDEBUG
  bmsg("Main kernel successful...\n");
#endif
  // Continue checking until nmax is reached.
  for(n = nmin; n < nmax; n += ld_kernel_nstep) {
    checkCUDAErr(clReleaseEvent(comp_done_event), "Release event object. (clReleaseEvent)");
    // Write the N parameter.
    checkCUDAErr(clSetKernelArg(check_more_ns_kernel, ARGNO_N, sizeof(n), (void *)&n),
      "Error: Setting kernel argument. (N)\n");

    // Note that comp_done_event is set many times.  Only the last one is kept after the loop.
    checkCUDAErr(clEnqueueNDRangeKernel(commandQueue,
          check_more_ns_kernel,
          1,
          NULL,
          global_cthread_count,
          NULL,
          0,
          NULL,
          &comp_done_event),
        "kernel 2 invocation");
  }
}

void get_factors_found(unsigned char *factor_found, const unsigned int cthread_count, const uint64_t start_t, int *check_ns_delay) {
  cl_event dev_read_event;
  /* Wait for the computation to finish execution */
  checkCUDAErr(clWaitForEvents(1, &comp_done_event), "Waiting for computation to finish. (clWaitForEvents)");
  checkCUDAErr(clReleaseEvent(comp_done_event), "Release event object 2. (clReleaseEvent)");

  // Get d_factor_found, into the thread'th factor_found array.
  checkCUDAErr(clEnqueueReadBuffer(commandQueue,
        d_factor_found,
        CL_TRUE,
        0,
        cthread_count*sizeof(cl_uchar),
        factor_found,
        0,
        NULL,
        &dev_read_event), "Retrieving results");

  checkCUDAErr(clWaitForEvents(1, &dev_read_event), "Waiting for results read. (clWaitForEvents)");
  checkCUDAErr(clReleaseEvent(dev_read_event), "Release event object 3. (clReleaseEvent)");

  //  cudaMemcpy(factor_found, d_factor_found, cthread_count*sizeof(unsigned char), cudaMemcpyDeviceToHost);
#ifndef NDEBUG
  printf("Retrieve successful...\n");
#endif
}
