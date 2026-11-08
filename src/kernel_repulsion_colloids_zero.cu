/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"
#include <math.h>

#include <stdio.h>
// Repulsion force is initialized to zero
__global__ void kernel_repulsion_colloids_zero(real* __restrict__ fx_rep_colloid,
					       real* __restrict__ fy_rep_colloid,
					       real* __restrict__ fz_rep_colloid) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;

  //-- Forces are initialized to zero --
  fx_rep_colloid[i] = 0;
  fy_rep_colloid[i] = 0;
  fz_rep_colloid[i] = 0;  
    
}
