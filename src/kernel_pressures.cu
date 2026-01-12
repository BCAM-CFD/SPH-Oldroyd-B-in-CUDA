#include "kernel_functions.h"
#include "config.h"

// Calculation of the pressure with an equation of state
__global__ void kernel_pressures(real* __restrict__ press,
				 real* __restrict__ dens,
				 real* __restrict__ mass) {
  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;
  
  
  press[i] = csq * (mass[i] * dens[i] - rho0);  
  
}
