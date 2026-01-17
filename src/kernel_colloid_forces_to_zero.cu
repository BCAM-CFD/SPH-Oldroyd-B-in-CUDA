/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// Forces are set to zero
__global__ void kernel_colloid_forces_to_zero(real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid)  {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;
  
  //-------- Forces are set to zero ---------
  fx_colloid[i] = 0.0;
  fy_colloid[i] = 0.0;
  tz_colloid[i] = 0.0;
  if (dim == 3) {
    fz_colloid[i] = 0.0;
    tx_colloid[i] = 0.0;
    ty_colloid[i] = 0.0;
  }
  
}
