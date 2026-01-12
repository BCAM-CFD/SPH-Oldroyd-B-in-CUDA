#include "kernel_functions.h"
#include "config.h"

//Part 2 of the function to move particles with a velocity Verlet method with lambda = 0.5.
__global__ void kernel_move_sys_VV_part2(real* __restrict__ vx,
					 real* __restrict__ vy,
					 real* __restrict__ vz,
					 real* __restrict__ cxx,
					 real* __restrict__ cxy,
					 real* __restrict__ cyy,
					 real* __restrict__ cxz,
					 real* __restrict__ cyz,
					 real* __restrict__ czz,
					 real* __restrict__ fx,
					 real* __restrict__ fy,
					 real* __restrict__ fz,
					 real* __restrict__ dcdt_xx,
					 real* __restrict__ dcdt_xy,
					 real* __restrict__ dcdt_yy,
					 real* __restrict__ dcdt_xz,
					 real* __restrict__ dcdt_yz,
					 real* __restrict__ dcdt_zz,
					 real* __restrict__ mass,
					 int*  __restrict__ type) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  if (type[i] == 0) { //  If it is from fluid
    real half_dt           = 0.5 * dt;    
    real half_dt_over_mass = half_dt / mass[i];    
    // Velocity at t + dt/2 
    vx[i] = vx[i] + half_dt_over_mass * fx[i];
    vy[i] = vy[i] + half_dt_over_mass * fy[i];

    // Conformation tensor at t + dt/2
    cxx[i] = cxx[i] + half_dt * dcdt_xx[i]; 
    cxy[i] = cxy[i] + half_dt * dcdt_xy[i];
    cyy[i] = cyy[i] + half_dt * dcdt_yy[i];

    //-- dim == 3 variables --    
    if (dim == 3) {
      vz[i] = vz[i] + half_dt_over_mass * fz[i];
      cxz[i] = cxz[i] + half_dt * dcdt_xz[i]; 
      cyz[i] = cyz[i] + half_dt * dcdt_yz[i];
      czz[i] = czz[i] + half_dt * dcdt_zz[i];
    }
    
  }
}

