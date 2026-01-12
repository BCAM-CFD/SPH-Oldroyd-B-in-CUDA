#include "kernel_functions.h"
#include "config.h"

// Forces are set to zero
__global__ void kernel_forces_to_zero(real* __restrict__ fx,
				      real* __restrict__ fy,
				      real* __restrict__ fz)  {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  //-------- Forces are set to zero ---------
  fx[i] = 0.0;
  fy[i] = 0.0;
  if (dim == 3)
    fz[i] = 0.0;
}
