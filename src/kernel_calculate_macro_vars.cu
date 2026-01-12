#include "kernel_functions.h"
#include "config.h"

//Function to calculate total kinetic energy and total momentum of the system
__global__ void kernel_calculate_macro_vars(real* __restrict__ mass,
					    real* __restrict__ vx,
					    real* __restrict__ vy,
					    real* __restrict__ vz,
					    real* __restrict__ kin_energy) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  real part_sum;

  if (dim == 2)
    part_sum = 0.5 * mass[i] * (vx[i] * vx[i] + vy[i] * vy[i]);
  else //--- dim = 3 ---
    part_sum = 0.5 * mass[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  atomicAdd(kin_energy, part_sum);

}
