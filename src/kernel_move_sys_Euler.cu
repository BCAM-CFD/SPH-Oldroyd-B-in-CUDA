#include "kernel_functions.h"
#include "config.h"

//Function to move particles with an Euler integrator.
// If you want to use this function, you should update it.
__global__ void kernel_move_sys_Euler(real* __restrict__ x,
				real* __restrict__ y,
				real* __restrict__ z,
				real* __restrict__ vx,
				real* __restrict__ vy,
				real* __restrict__ vz,				
				real* __restrict__ fx,
				real* __restrict__ fy,
				real* __restrict__ fz,
				real* __restrict__ mass) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  x[i] = x[i] + vx[i] * dt;
  y[i] = y[i] + vy[i] * dt;
  if (dim == 3)
    z[i] = z[i] + vz[i] * dt;
    
  vx[i] = vx[i] + fx[i]/mass[i] * dt;
  vy[i] = vy[i] + fy[i]/mass[i] * dt;
  if (dim == 3)
    vz[i] = vz[i] + fz[i]/mass[i] * dt;

  // Periodic Boundary conditions
  if (x[i] < 0)
    x[i] = x[i] + L[0];
  if (wall == 0) //--- If no wall ---  
    if (y[i] < 0)
      y[i] = y[i] + L[1];
  if (dim == 3)
    if (z[i] < 0)
      z[i] = z[i] + L[2];

  if (x[i] > L[0])
    x[i] = x[i] - L[0];
  if (wall == 0) //--- If no wall ---  
    if (y[i] > L[1])
      y[i] = y[i] - L[1];
  if (dim == 3)
    if (z[i] > L[2])
      z[i] = z[i] - L[2];      
  
}
