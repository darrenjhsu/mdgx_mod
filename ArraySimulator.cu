#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <nvml.h>
#include <curand_kernel.h>

#include "MultiSimDS.h"
#include "KernelMacros.h"

__device__ __constant__ gpuMultiSim cGms;

//-----------------------------------------------------------------------------
// GetGpuSpecs: obtain specs for the GPU in the system.  This is here to have
//              all of the CUDA-specific data structures in this one unit,
//              while returning the information that's really needed in a
//              custom struct full of POD types to be incorporated into the
//              main peptide control data structure in Peptide.c.
//
// Arguments:
//   reckless:   taken from the trajectory data structure and fed in from the
//               command line, lets mdgx go ahead and run on GPUs that are
//               already in use
//-----------------------------------------------------------------------------
extern "C" gpuSpecs GetGpuSpecs(int reckless)
{
  int i, ndev, gpucount, seldev;
  unsigned int nvmlItemCount;
  int* validGpus;
  int* gpuList;
  cudaError_t stt;
  nvmlReturn_t sttNV;
  cudaDeviceProp devPRP;
  gpuSpecs devspc;
  nvmlProcessInfo_t nvmlInfo[32];
  nvmlDevice_t ntdev;
  
  // Test that there is a GPU in the system
  stt = cudaGetDeviceCount(&gpucount);
  if (gpucount == 0) {
    printf("mdgx >> Error.  No CUDA-capable devices were found.\n");
    cudaDeviceReset();
    exit(1);
  }

  // Activate zero-copy
  cudaSetDeviceFlags(cudaDeviceMapHost);

  // Initialize the NVIDIA Management Library
  nvmlInit();
  
  // Get device properties
  validGpus = (int*)malloc(gpucount * sizeof(int));
  gpuList = (int*)malloc(gpucount * sizeof(int));
  ndev = 0;
  for (i = 0; i < gpucount; i++) {
    cudaGetDeviceProperties(&devPRP, i);
    if (devPRP.major >= 3) {
      nvmlDeviceGetHandleByIndex(i, &ntdev);
      nvmlItemCount = 0;
      sttNV = nvmlDeviceGetComputeRunningProcesses(ntdev, &nvmlItemCount,
                                                   nvmlInfo);
      if (sttNV != NVML_SUCCESS && sttNV != NVML_ERROR_INSUFFICIENT_SIZE) {
	printf("mdgx >> Warning.  Unable to monitor activity on GPU %d "
	       "[error %u]\n", i, sttNV);
      }
      if (nvmlItemCount == 0 || reckless == 1) {
        validGpus[i] = 1;
        gpuList[ndev] = i;
        ndev++;
      }
    }
  }
  if (ndev == 0 && reckless == 0) {
    printf("mdgx >> All GPUs are unavailable, or assisting other customers.  "
           "If you believe\nmdgx >> you have received this message in error, "
	   "or if you know people are\nmdgx >> using the GPUs and just don't "
	   "care, you may re-run mdgx with the\nmdgx >> -Reckless flag.  What "
           "your colleagues say at the water cooler is\nmdgx >> not the "
	   "responsibility of Amber developers.\n");
    exit(1);
  }

  // Shut down the NVIDIA Management Lbirary
  nvmlShutdown();
  
  // Select a device from the list
  stt = cudaSetValidDevices(gpuList, ndev);
  if (stt != cudaSuccess) {
    printf("mdgx >> Error searching for CUDA-compatible GPU.\n");
    cudaDeviceReset();
    exit(1);
  }

  // Establish the CUDA context
  stt = cudaFree(0);
  if (stt != cudaSuccess) {
    printf("mdgx >> Error selecting compatible GPU.\n");
    cudaDeviceReset();
    exit(1);
  }

  // Get the device
  stt = cudaGetDevice(&seldev);
  if (stt != cudaSuccess) {
    printf("mdgx >> Error setting cuda device.\n");
    cudaDeviceReset();
    exit(1);
  }
  cudaDeviceSynchronize();
  cudaGetDeviceProperties(&devPRP, seldev);

  // Copy the relevant information for shipment back to the calling function
  devspc.major          = devPRP.major;
  devspc.minor          = devPRP.minor;
  devspc.MPcount        = devPRP.multiProcessorCount;
  devspc.maxThrPerMP    = devPRP.maxThreadsPerMultiProcessor;
  devspc.maxThrPerBlock = devPRP.maxThreadsPerBlock;
  devspc.cardMemory     = devPRP.totalGlobalMem;
  i = strlen(devPRP.name);
  if (i > 127) {
    i = 127;
  }
  strncpy(devspc.name, devPRP.name, i);
  devspc.name[i] = '\0';
  
  // Free allocated memory
  free(gpuList);
  free(validGpus);

  return devspc;
}

//-----------------------------------------------------------------------------
// kGpuPRNGSetup: kernel for initializing GPU random number generators
//-----------------------------------------------------------------------------
__global__ void kGpuPRNGSetup(curandState_t *states, int igseed)
{
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  curand_init(igseed, tid, 0, &states[tid]);
}

//-----------------------------------------------------------------------------
// InitGpuPRNG: initialize pseudo-random number generators on the GPU.
//
// Arguments:
//   gms:       the repository for all parameters amd coordinates
//   igseed:    the random number generator seed
//   nblocks:   the number of blocks that the main dynamics kernels will run
//   blockDim:  dimension of the main dynamics kernel blocks
//-----------------------------------------------------------------------------
extern "C" void InitGpuPRNG(gpuMultiSim *gms, int igseed, int nblocks,
			    int blockdim)
{
  cudaMalloc((void **)&gms->prngStates,
             nblocks * blockdim * sizeof(curandState));
  kGpuPRNGSetup<<<nblocks, blockdim>>>((curandState_t*)gms->prngStates,
                                       igseed);
}

//-----------------------------------------------------------------------------
// SetGmsImage: function to establish a GPU Multi-Simulator on the device,
//              with pointers to all of the device-allocated memory as well as
//              constants describing the simulation conditions.
//
// Arguments:
//   gms:    the repository for all parameters amd coordinates
//-----------------------------------------------------------------------------
extern "C" void SetGmsImage(gpuMultiSim *gms)
{
  cudaError_t status;

  status = cudaMemcpyToSymbol(cGms, gms, sizeof(gpuMultiSim));
  if (status != cudaSuccess) {
    printf("SetGmsImage >> Unable to copy gpuMultiSim struct to the "
           "device (error %d).\n", (int)status);
    exit(1);
  }
}

//----------------------------------------------------------------------------
// kSetSystemCounters: set counters to guide the blocks as they step through
//                     systems during each segment of dynamics, in between
//                     coordinate writes.
//-----------------------------------------------------------------------------
__global__ void kSetSystemCounters(int blocks)
{
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < 2 * cGms.nsgmdout) {
    cGms.DVCsystemPos[tid] = blocks;
    tid += gridDim.x * blockDim.x;
  }
}

//-----------------------------------------------------------------------------
// Dynamics kernels with RATTLE to compute forces (kDynLoop) or force and
// energies (kEStep).  The concept is that kEStep will launch to get forces and
// energies for the first step, then kDynLoop will fire off for (ntpr - 1)
// steps, then kEStep will launch again, and so on until the maximum number of
// steps has been reached.
//-----------------------------------------------------------------------------
#define GO_RATTLE
#define ATOM_LIMIT SM_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 256
#else
#  define THREAD_COUNT 288
#endif
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kEStepSmallRtt(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kDynLoopSmallRtt(int sgc)
#include "kDynamics.h"
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kEStepSmallGBRtt(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kDynLoopSmallGBRtt(int sgc)
#include "kDynamics.h"
#undef GBSOLVENT
#undef THREAD_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT MD_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 512
#else
#  define THREAD_COUNT 576
#endif
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kEStepMedRtt(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kDynLoopMedRtt(int sgc)
#include "kDynamics.h"
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kEStepMedGBRtt(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kDynLoopMedGBRtt(int sgc)
#include "kDynamics.h"
#undef GBSOLVENT
#undef THREAD_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT LG_ATOM_COUNT
#define THREAD_COUNT 1024
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kEStepLargeRtt(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kDynLoopLargeRtt(int sgc)
#include "kDynamics.h"
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kEStepLargeGBRtt(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kDynLoopLargeGBRtt(int sgc)
#include "kDynamics.h"
#undef GBSOLVENT
#undef THREAD_COUNT
#undef ATOM_LIMIT
#undef GO_RATTLE

//-----------------------------------------------------------------------------
// Kernels without RATTLE.  The concept behind DynLoop and EStep kernels is the
// same as above, but without the register burden of RATTLE the kernels can
// engage additional threads.
//-----------------------------------------------------------------------------
#define ATOM_LIMIT SM_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 256
#else
#  define THREAD_COUNT 320
#endif
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kEStepSmall(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kDynLoopSmall(int sgc)
#include "kDynamics.h"
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kEStepSmallGB(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 4) kDynLoopSmallGB(int sgc)
#include "kDynamics.h"
#undef GBSOLVENT
#undef THREAD_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT MD_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 512
#else
#  define THREAD_COUNT 640
#endif
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kEStepMed(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kDynLoopMed(int sgc)
#include "kDynamics.h"
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kEStepMedGB(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 2) kDynLoopMedGB(int sgc)
#include "kDynamics.h"
#undef GBSOLVENT
#undef THREAD_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT LG_ATOM_COUNT
#define THREAD_COUNT 1024
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kEStepLarge(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kDynLoopLarge(int sgc)
#include "kDynamics.h"
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kEStepLargeGB(int sgc)
#include "kDynamics.h"
#undef COMPUTE_ENERGY
__global__ void __launch_bounds__(THREAD_COUNT, 1) kDynLoopLargeGB(int sgc)
#include "kDynamics.h"
#undef GBSOLVENT
#undef THREAD_COUNT
#undef ATOM_LIMIT
  
//-----------------------------------------------------------------------------
// LaunchDynamics: launch the appropriate kernels for energy and forces.  As
//                 is done in the pmemd code, this function in the CUDA unit
//                 encapsualtes the launch so that the .c libraries can be
//                 built with a standard C compiler.
//
// Arguments:
//   gms:       the repository for all parameters amd coordinates
//   blockDim:  the block size to use, determined by PlanGpuUtilization in
//              Peptide.c
//   devspc:    device specifications
//-----------------------------------------------------------------------------
extern "C" void LaunchDynamics(gpuMultiSim *gms, int blockDim, int nblocks,
			       gpuSpecs *devspc)
{
  int i;

  // Initialize system counters for this portion of dynamics
  kSetSystemCounters<<<nblocks, blockDim>>>(nblocks);
  
  // Vacuum-phase dynamics
  if (gms->igb == 6) {
    if (gms->rattle == 0) {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < 512) {
          kEStepSmall<<<nblocks, blockDim>>>(i);
          kDynLoopSmall<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < 1024) {
          kEStepMed<<<nblocks, blockDim>>>(i);
          kDynLoopMed<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLarge<<<nblocks, blockDim>>>(i);
          kDynLoopLarge<<<nblocks, blockDim>>>(i);
        }
      }
    }
    else {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < 512) {
          kEStepSmallRtt<<<nblocks, blockDim>>>(i);
          kDynLoopSmallRtt<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < 1024) {
          kEStepMedRtt<<<nblocks, blockDim>>>(i);
          kDynLoopMedRtt<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLargeRtt<<<nblocks, blockDim>>>(i);
          kDynLoopLargeRtt<<<nblocks, blockDim>>>(i);
        }
      }
    }
  }

  // Dynamics in Generalized Born solvent
  else {
    if (gms->rattle == 0) {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < 512) {
          kEStepSmallGB<<<nblocks, blockDim>>>(i);
          kDynLoopSmallGB<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < 1024) {
          kEStepMedGB<<<nblocks, blockDim>>>(i);
          kDynLoopMedGB<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLargeGB<<<nblocks, blockDim>>>(i);
          kDynLoopLargeGB<<<nblocks, blockDim>>>(i);
	}
      }
    }
    else {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < 512) {
          kEStepSmallGBRtt<<<nblocks, blockDim>>>(i);
          kDynLoopSmallGBRtt<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < 1024) {
          kEStepMedGBRtt<<<nblocks, blockDim>>>(i);
          kDynLoopMedGBRtt<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLargeGBRtt<<<nblocks, blockDim>>>(i);
          kDynLoopLargeGBRtt<<<nblocks, blockDim>>>(i);
	}
      }
    }
  }
}
