/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"

// Function to initialize cell colloids. This is important if there could be 0 particles
// in one given cell. This could happen in the interior of the colloid particles, where
// there are no particles

__global__ void kernel_initialize_cells(int* cell_start,
					int* cell_end) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= Ntotal_cells) return;
    
    cell_start[i] = -1;
    cell_end[i]   = -2;
}
