/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

//Function to set total kinetic energy and total momentum to zero
__global__ void kernel_macro_vars_to_zero(real* __restrict__ kin_energy) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  *kin_energy = 0.0;
}
