//-----------------------------------------------------------------------------
// Useful macros for kernel operations on the GPU
//-----------------------------------------------------------------------------
#ifndef KernelMacroFramework
#define KernelMacroFramework

#if defined(CUDA_VERSION) && (CUDA_VERSION < 9000)
#  define __shfl_down_sync(a, b, c, d)  __shfl_down(b, c, d)
#  define __shfl_xor_sync(a, b, c, d)   __shfl_xor(b, c, d)
#  define __shfl_up_sync(a, b, c, d)    __shfl_up(b, c, d)
#  define __shfl_sync(a, b, c, d)       __shfl(b, c, d)
#  define __ballot_sync(a, b)           __ballot(a, b)
#  define __syncwarp(a)                 
#endif

#  include "GpuDS.h"

//-----------------------------------------------------------------------------
// WarpREDUCE: warp reduction to sum all values into the first lane
//-----------------------------------------------------------------------------
#define WarpREDUCE(var) \
{ \
  var += __shfl_down_sync(0xffffffff, var, 16, 32); \
  var += __shfl_down_sync(0xffffffff, var,  8, 32); \
  var += __shfl_down_sync(0xffffffff, var,  4, 32); \
  var += __shfl_down_sync(0xffffffff, var,  2, 32); \
  var += __shfl_down_sync(0xffffffff, var,  1, 32); \
}

//-----------------------------------------------------------------------------
// DvcCROSSPf: compute the cross product vA x vB = vC
//-----------------------------------------------------------------------------
#define DvcCROSSPf(vA, vB, vC) \
{ \
  vC[0] = vA[1]*vB[2] - vA[2]*vB[1]; \
  vC[1] = vA[2]*vB[0] - vA[0]*vB[2]; \
  vC[2] = vA[0]*vB[1] - vA[1]*vB[0]; \
}

#define ISNAN(x) ((__float_as_uint(x) & 0x7f000000) == 0x7f000000 &&	\
		  (__float_as_uint(x) & 0xffffff) > 0)

#endif
