#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>

// Function to collect forces and torques of computational particles of colloids
__global__ void kernel_collect_force_colloids(real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid,
					      real* __restrict__ fx,
					      real* __restrict__ fy,
					      real* __restrict__ fz,
					      real* __restrict__ x_center,
					      real* __restrict__ y_center,
					      real* __restrict__ z_center,
					      int*  __restrict__ colloids_list,
					      int*  __restrict__ colloids_start) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;
  // i is the identity of the colloid

  real fx_colloid_i = fx_colloid[i];
  real fy_colloid_i = fy_colloid[i];
  real fz_colloid_i = fz_colloid[i];
  real tx_colloid_i = tx_colloid[i];
  real ty_colloid_i = ty_colloid[i];
  real tz_colloid_i = tz_colloid[i];  
  int colloids_end_i;
  if (i < N_colloids - 1)
    colloids_end_i = colloids_start[i+1] - 1;
  else  // i == Nwalls
    colloids_end_i = Nlist_colloids - 1;

  for (int j = colloids_start[i]; j <= colloids_end_i; ++j) {
    int part = colloids_list[j];
    
    //--- Forces ----
    fx_colloid_i = fx_colloid_i + fx[part];
    fy_colloid_i = fy_colloid_i + fy[part];
    if (dim == 3)
      fz_colloid_i = fz_colloid_i + fz[part];
    //--- Torques ---
    tz_colloid_i = tz_colloid_i + fy[part] * x_center[j] - fx[part] * y_center[j];
    if (dim == 3) {
      tx_colloid_i = tx_colloid_i + fz[part] * y_center[j] - fy[part] * z_center[j];
      ty_colloid_i = ty_colloid_i + fx[part] * z_center[j] - fz[part] * x_center[j];
    }
  }

  // The result is storaged in the GPU
  fx_colloid[i] = fx_colloid_i;
  fy_colloid[i] = fy_colloid_i;
  if (dim == 3)
    fz_colloid[i] = fz_colloid_i;

  tz_colloid[i] = tz_colloid_i;
  if (dim == 3) {
    tx_colloid[i] = tx_colloid_i;
    ty_colloid[i] = ty_colloid_i;
  }

}

