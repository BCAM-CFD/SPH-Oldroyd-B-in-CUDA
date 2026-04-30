/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

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
					      real* __restrict__ x,
					      real* __restrict__ y,
					      real* __restrict__ z,
					      real* __restrict__ fx,
					      real* __restrict__ fy,
					      real* __restrict__ fz,
					      real* __restrict__ coll_x,
					      real* __restrict__ coll_y,
					      real* __restrict__ coll_z,
					      int*  __restrict__ colloids_list,
					      int*  __restrict__ colloids_start) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;
  // i is the identity of the colloid

  real fx_colloid_i = fx_colloid[i];
  real fy_colloid_i = fy_colloid[i];
  real tx_colloid_i = tx_colloid[i];
  real fz_colloid_i;
  real ty_colloid_i;
  real tz_colloid_i;  
  if (dim == 3) {
    fz_colloid_i = fz_colloid[i];
    ty_colloid_i = ty_colloid[i];
    tz_colloid_i = tz_colloid[i];
  }
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
    // Relative position of boundary particle with respect to colloid center ---
    //    x_center, y_center,... cannot be used, because they are not changed during the sim
    real rx = x[part] - coll_x[i];
    real ry = y[part] - coll_y[i];
    real rz;
    if (dim == 3)
      rz = z[part] - coll_z[i];
    // Minimum-image correction in periodic directions
    if (rx >  0.5 * L[0]) rx -= L[0];
    if (rx < -0.5 * L[0]) rx += L[0];
    if (wall == 0) {
      if (ry >  0.5 * L[1]) ry -= L[1];
      if (ry < -0.5 * L[1]) ry += L[1];
    }
    if (dim == 3) {
      rz = z[part] - coll_z[i];
      if (rz >  0.5 * L[2]) rz -= L[2];
      if (rz < -0.5 * L[2]) rz += L[2];
    }
    // Torques: T = r x F 
    tz_colloid_i = tz_colloid_i + rx * fy[part] - ry * fx[part];
    if (dim == 3) {
      tx_colloid_i = tx_colloid_i + ry * fz[part] - rz * fy[part];
      ty_colloid_i = ty_colloid_i + rz * fx[part] - rx * fz[part];        
    
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

