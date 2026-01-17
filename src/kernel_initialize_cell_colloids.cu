/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"

#include <stdio.h>

// Function to initialize cell colloids. This is important if there could be 0 or 1 particles
// in one given cell.
__global__ void kernel_initialize_cell_colloids(int* __restrict__ coll_cell_start,
						int* __restrict__ coll_cell_end) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= Ntotal_cells_colloids) return;

    //-- If these variables remain like this, that means there are no particles in the cell -
    // if coll_cell_end < coll_cell_start the loops of neighbours will not be executed
    coll_cell_start[i] = -1;
    coll_cell_end[i]   = -2;
    
}

