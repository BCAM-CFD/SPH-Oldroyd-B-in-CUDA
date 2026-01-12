#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>

// Function to calculate the cell list for colloids (before sorting it)
__global__ void kernel_cell_colloids_list(real* __restrict__ coll_x,
					  real* __restrict__ coll_y,
					  real* __restrict__ coll_z,
					  int*  __restrict__ coll_cell,
					  int*  __restrict__ coll_index) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;
  
  // We calculate in which cell is each particle
  int cx, cy, cz;
  cx = floor(coll_x[i] / cell_colloids_size[0]);
  cy = floor(coll_y[i] / cell_colloids_size[1]);
  if ( dim == 3 )
    cz = floor(coll_z[i] / cell_colloids_size[2]);
  if (dim == 2)
    coll_cell[i] = cy * Ncells_colloids[0] + cx;
  else // dim == 3
    coll_cell[i] = cz * Ncells_colloids[0] * Ncells_colloids[1] +
      cy * Ncells_colloids[0] + cx;

  // Particle_index array
  coll_index[i] = i;

}
