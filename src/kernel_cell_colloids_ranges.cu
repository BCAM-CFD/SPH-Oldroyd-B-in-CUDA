#include "kernel_functions.h"

#include <stdio.h>

// Function to calculate start and end indices of the cell_list
__global__ void kernel_cell_colloids_ranges(int* __restrict__ coll_cell,
					    int* __restrict__ coll_cell_start,
					    int* __restrict__ coll_cell_end) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= N_colloids) return;

    int current_cell = coll_cell[i];

    // It detects the beginning of a cell
    if (i == 0 || coll_cell[i - 1] != current_cell)
        coll_cell_start[current_cell] = i;

    // It detects the end of a cell
    if (i == N_colloids - 1 || coll_cell[i + 1] != current_cell)
        coll_cell_end[current_cell] = i;
    
}
