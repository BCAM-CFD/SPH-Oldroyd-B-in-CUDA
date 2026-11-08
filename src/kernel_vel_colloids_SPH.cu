/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>

// Colloids velocities are calculated at time t + Delta t,
// according to the collected SPH force
// V' = Vn + F_sph * Delta t
__global__ void kernel_vel_colloids_SPH(real* __restrict__ coll_vx,
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
    // i is the identity of the colloid 

  //--- Some initialization ---
  real dt_over_mass = dt / coll_mass;
  real dt_over_I    = dt / coll_I;  

  //--- Velocity at t + dt  ---
  coll_vx[i] += dt_over_mass * fx_colloid[i];
  coll_vy[i] += dt_over_mass * fy_colloid[i];
  if (dim == 3)
    coll_vz[i] += dt_over_mass * fz_colloid[i];

  //--- Angular Velocity at t + dt  ---
  coll_omegaz[i] += dt_over_I * tz_colloid[i];
  if (dim == 3) {
    coll_omegax[i] += dt_over_I * tx_colloid[i];
    coll_omegay[i] += dt_over_I * ty_colloid[i];
  }

}
