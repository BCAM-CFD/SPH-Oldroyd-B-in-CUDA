#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>

// Function to collect the forces on walls
__global__ void kernel_collect_force_walls(real* __restrict__ fx_wall,
					   real* __restrict__ fy_wall,
					   real* __restrict__ fz_wall,
					   real* __restrict__ fx,
					   real* __restrict__ fy,
					   real* __restrict__ fz,				    
				   int*  __restrict__ walls_list,
				   int*  __restrict__ walls_start) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= Nwalls) return;
  // i is the identity of the wall

  real fx_wall_i = 0.0;
  real fy_wall_i = 0.0;
  real fz_wall_i = 0.0;
  int walls_end_i;
  if (i < Nwalls - 1)
    walls_end_i = walls_start[i+1] - 1;
  else  // i == Nwalls
    walls_end_i = Nlist_walls - 1;
    
  for (int j = walls_start[i]; j <= walls_end_i; ++j) {
    int part = walls_list[j];
    fx_wall_i = fx_wall_i + fx[part];
    fy_wall_i = fy_wall_i + fy[part];
    if (dim == 3)
      fz_wall_i = fz_wall_i + fz[part];
  }

  // The result is storaged in the GPU
  fx_wall[i] = fx_wall_i;
  fy_wall[i] = fy_wall_i;
  if (dim == 3)
    fz_wall[i] = fz_wall_i;

}

