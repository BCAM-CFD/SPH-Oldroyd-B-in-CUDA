/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// Forces and densities are set to zero
__global__ void kernel_zero(real* __restrict__ fx,
			    real* __restrict__ fy,
			    real* __restrict__ fz,
			    real* __restrict__ dens)  
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  //-------- Forces are set to zero ---------
  fx[i] = 0.0;
  fy[i] = 0.0;
  if (dim == 3)
    fz[i] = 0.0;

  //------ Densities are set to zero ------
  dens[i] = 0.0;
}
