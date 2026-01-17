/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"
#include <stdio.h>

__global__ void kernel_print_particles(real* __restrict__ x,
				       real* __restrict__ y,
				       real* __restrict__ z,
				       real* __restrict__ vx,
				       real* __restrict__ vy,
				       real* __restrict__ vz,
				       real* __restrict__ fx,
				       real* __restrict__ fy,
				       real* __restrict__ fz,
				       real* __restrict__ dens,
				       real* __restrict__ press,
				       int* __restrict__ particle_cell,
				       int* __restrict__ particle_index) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i >= N || i >= 1000) return;
  
  if (dim == 3)
    printf("%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " %d %d\n",
	   i,
	   x[i],
	   y[i],
	   z[i],
	   vx[i],
	   vy[i],
	   vz[i],
	   fx[i],
	   fy[i],
	   fz[i],	   	   
	   dens[i],
	   press[i],
	   particle_cell[i],
	   particle_index[i]);
  else // dim = 2
    printf("%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " %d %d\n",
	   i,
	   x[i],
	   y[i],
	   vx[i],
	   vy[i],
	   fx[i],
	   fy[i],	   	   
	   dens[i],
	   press[i],
	   particle_cell[i],
	   particle_index[i]);

}
