/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include <stdio.h>

__global__ void kernel_print_cell_ranges(int* __restrict__ cell_start,
					 int* __restrict__ cell_end) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= Ntotal_cells) return;
  
  printf("%d %d %d\n",
	 i,
	 cell_start[i],
	 cell_end[i]);
}

__global__ void kernel_cell_test(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 int* __restrict__ particle_cell,
				 int* __restrict__ particle_index,
				 int* __restrict__ cell_start,
				 int* __restrict__ cell_end) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= 1) return;

  int part = 29;
  int cx_i, cy_i, cz_i;
  int cx, cy, cz;  
  cx_i = floor(x[part] / cell_size[0]);
  cy_i = floor(y[part] / cell_size[1]);
  if ( dim == 3 )
    cz_i = floor(z[part] / cell_size[2]);

  if (dim == 2)
    for (int dx = -1; dx <= 1; ++dx) 
      for (int dy = -1; dy <= 1; ++dy)
	{
	  cx = cx_i + dx;
	  cy = cy_i + dy;
	  // periodic boundary conditions
	  if (cx < 0)
	    cx = cx + Ncells[0];
	  if (cy < 0)
	    cy = cy + Ncells[1];
	  if (cx >= Ncells[0])
	    cx = cx - Ncells[0];
	  if (cy >= Ncells[1])	  
	    cy = cy - Ncells[1];      
	  int neigh_cell = cy * Ncells[0] + cx;
	  for (int k = cell_start[neigh_cell]; k <= cell_end[neigh_cell]; ++k) {
	    
	    int j = particle_index[k]; // índice de la partícula vecina
	    
	    printf("%d %d %f %f %f %f\n", j, particle_cell[k], x[j], y[j], x[part], y[part]); 
	  }
	}
  else
    for (int dx = -1; dx <= 1; ++dx) 
      for (int dy = -1; dy <= 1; ++dy)
	for (int dz = -1; dz <= 1; ++dz)	
	  {
	    cx = cx_i + dx;
	    cy = cy_i + dy;
	    cz = cz_i + dz;	    
	    // periodic boundary conditions
	    if (cx < 0)
	      cx = cx + Ncells[0];
	    if (cy < 0)
	      cy = cy + Ncells[1];
	    if (cz < 0)
	      cz = cz + Ncells[2];	    
	    if (cx >= Ncells[0])
	      cx = cx - Ncells[0];
	    if (cy >= Ncells[1])	  
	      cy = cy - Ncells[1];
	    if (cz >= Ncells[2])	  
	      cz = cz - Ncells[2];      	    
	    int neigh_cell = cz * Ncells[0] * Ncells[1] + cy * Ncells[0] + cx;
	    for (int k = cell_start[neigh_cell]; k <= cell_end[neigh_cell]; ++k) {
	      
	      int j = particle_index[k]; // índice de la partícula vecina
	      
	      printf("%d %d %f %f %f %f\n", j, particle_cell[k], x[j], y[j], x[part], y[part]); 
	    }
	  }  

}
