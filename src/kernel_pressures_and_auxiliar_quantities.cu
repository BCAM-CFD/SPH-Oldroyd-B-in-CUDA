/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// Calculation of the pressure with an equation of state
__global__ void kernel_pressures_and_auxiliar_quantities(real* __restrict__ press,
							 real* __restrict__ dens,
							 real* __restrict__ mass,
							 real* __restrict__ inv_dens) {
  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  real dens_i = dens[i];
  
  press[i] = csq * (mass[i] * dens_i - rho0);
  inv_dens[i] = 1.0/dens_i;
  
}
