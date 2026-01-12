#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>
//Part 2 of the function to move colloids with a velocity Verlet method with lambda = 0.5.
__global__ void kernel_move_colloids_VV_part2(real* __restrict__ coll_vx,
					      real* __restrict__ coll_vy,
					      real* __restrict__ coll_vz,
					      real* __restrict__ coll_omegax,
					      real* __restrict__ coll_omegay,
					      real* __restrict__ coll_omegaz,
					      real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;

  //--- Some initializations ---
  real half_dt_over_mass = 0.5 * dt / coll_mass;
  real half_dt_over_I    = 0.5 * dt / coll_I;

  // Velocity at t+dt
  coll_vx[i] += half_dt_over_mass * fx_colloid[i];
  coll_vy[i] += half_dt_over_mass * fy_colloid[i];
  if (dim == 3)
    coll_vz[i] += half_dt_over_mass * fz_colloid[i];

  // Angular velocity at t+dt
  coll_omegaz[i] += half_dt_over_I * tz_colloid[i];          
  if (dim == 3) {
    coll_omegax[i] += half_dt_over_I * tx_colloid[i];
    coll_omegay[i] += half_dt_over_I * ty_colloid[i];
  }
    
}
