/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// Densities are set to zero
__global__ void kernel_density_to_zero(real* __restrict__ dens)  
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  //------ Densities are set to zero ------
  dens[i] = 0.0;
}
