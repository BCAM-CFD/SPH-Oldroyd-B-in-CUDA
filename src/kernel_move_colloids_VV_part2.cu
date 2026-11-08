/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>
//Part 2 of the function to move colloids with a velocity Verlet method with lambda = 0.5.
// The final velocity is the velocity at t + dt
// V_n+1/2 = V_tilde + 1/2 * F_repulsion * dt / m
__global__ void kernel_move_colloids_VV_part2(real* __restrict__ coll_vx,
					      real* __restrict__ coll_vy,
					      real* __restrict__ coll_vz,
					      real* __restrict__ coll_omegax,
					      real* __restrict__ coll_omegay,
					      real* __restrict__ coll_omegaz,
					      real* __restrict__ fx_rep_colloid,
					      real* __restrict__ fy_rep_colloid,
					      real* __restrict__ fz_rep_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;

  //--- Some initializations ---
  real half_dt_over_mass = 0.5 * dt / coll_mass;
 
  // Velocity at t+dt
  coll_vx[i] += half_dt_over_mass * fx_rep_colloid[i];
  coll_vy[i] += half_dt_over_mass * fy_rep_colloid[i];
  if (dim == 3)
    coll_vz[i] += half_dt_over_mass * fz_rep_colloid[i];

  //--- Angular Velocity at t + dt was already calculated at kernel_vel_colloids_SPH.
  //    If a rotational interaction between colloids is programmed, the angular
  //    velocity change at t + dt due to such an interaction should be calculated here --  
    
}
