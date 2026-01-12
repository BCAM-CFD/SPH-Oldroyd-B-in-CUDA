#include "kernel_functions.h"
#include "config.h"

//Part 1 of the function to move particles with a velocity Verlet method with lambda = 0.5.
__global__ void kernel_move_sys_VV_part1(real* __restrict__ x,
					 real* __restrict__ y,
					 real* __restrict__ z,
					 real* __restrict__ vx,
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

  int type_i = type[i];
  if (type_i > 2) return; // particles of colloids are not moved here

  real xi  = x[i];
  real yi  = y[i];
  real zi;
  if (dim == 3)  
    zi  = z[i];
  real vxi;
  real vyi;
  real vzi;

  if (type_i == 0) { // i fluid 
    vxi = vx[i];
    vyi = vy[i];
    if (dim == 3)
      vzi = vz[i];

    real half_dt           = 0.5 * dt;    
    real half_dt_over_mass = half_dt / mass[i];
    // Velocity at t + dt/2 
    vx[i] = vxi + half_dt_over_mass * fx[i];
    vy[i] = vyi + half_dt_over_mass * fy[i];    

    // Conformation tensor at t + dt/2
    cxx[i] = cxx[i] + half_dt * dcdt_xx[i]; 
    cxy[i] = cxy[i] + half_dt * dcdt_xy[i];
    cyy[i] = cyy[i] + half_dt * dcdt_yy[i];

    //-- dim == 3 variables --
    if (dim == 3) {
      vz[i]  = vz[i]  + half_dt_over_mass * fz[i];
      cxz[i] = cxz[i] + half_dt * dcdt_xz[i]; 
      cyz[i] = cyz[i] + half_dt * dcdt_yz[i];
      czz[i] = czz[i] + half_dt * dcdt_zz[i];
    }
    
  }
  else if (type_i == 1) {// i bottom wall
    vxi = V_bottom;
    vyi = 0.0;
    if (dim == 3)    
      vzi = 0.0;
  }
  else if (type_i == 2) {// i top wall
    vxi = V_top;
    vyi = 0.0;
    if (dim == 3)    
      vzi = 0.0;
  }  //Particles with type_i > 2 where discarded above

  // Position at t + dt
  xi = xi + vxi * dt;
  yi = yi + vyi * dt;
  if (dim == 3)
    zi = zi + vzi * dt;
  
  // Periodic Boundary conditions
  if (xi < 0)
    xi = xi + L[0];
  if (wall == 0) //--- If no wall ---
    if (yi < 0)
      yi = yi + L[1];
  if (dim == 3)
    if (zi < 0)
      zi = zi + L[2];
  
  if (xi > L[0])
    xi = xi - L[0];
  if (wall == 0) //--- If no wall ---    
    if (yi > L[1])
      yi = yi - L[1];
  if (dim == 3)
    if (zi > L[2])
      zi = zi - L[2];

  // Data is stored  (v[i] is stored above, only when needed (when type_i == 0)
  x[i] = xi;
  y[i] = yi;
  if (dim == 3)
    z[i] = zi;    
}
