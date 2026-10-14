/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// This function copies positions (in float) and type.
__global__ void kernel_build_posf(real* __restrict__ x,
                                  real* __restrict__ y,
                                  real* __restrict__ z,
                                  int*  __restrict__ type,
                                  int*  __restrict__ particle_index,
				  float4* __restrict__ posf) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;
  
  int j = particle_index[i];

  float zj;
  if (dim == 3)
    zj = (float) z[j];
  posf[i] = make_float4((float) x[j], (float) y[j], zj, __int_as_float(type[j]));
  // The type can be recovered with __float_as_int function
}
