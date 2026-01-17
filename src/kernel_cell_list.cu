/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>

// Function to calculate the cell list (before sorting it)
__global__ void kernel_cell_list(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 int*  __restrict__ particle_cell,
				 int*  __restrict__ particle_index) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;
  
  // We calculate in which cell is each particle
  int cx, cy, cz;
  cx = floor(x[i] / cell_size[0]);
  if (wall == 0) //-- No wall -- 
    cy = floor(y[i] / cell_size[1]);
  else //-- With wall --
    cy = floor((y[i] + wall_width) / cell_size[1]);    
  if ( dim == 3 )
    cz = floor(z[i] / cell_size[2]);
  if (dim == 2)
    particle_cell[i] = cy * Ncells[0] + cx;
  else // dim == 3
    particle_cell[i] = cz * Ncells[0] * Ncells[1] + cy * Ncells[0] + cx;

  // Particle_index array
  particle_index[i] = i;

}

