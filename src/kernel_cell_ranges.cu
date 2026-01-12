#include "kernel_functions.h"

#include <stdio.h>

// Function to calculate start and end indices of the cell_list
__global__ void kernel_cell_ranges(int* __restrict__ particle_cell,
				   int* __restrict__ cell_start,
				   int* __restrict__ cell_end) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= N) return;

    int current_cell = particle_cell[i];

    // It detects the beginning of a cell
    if (i == 0 || particle_cell[i - 1] != current_cell)
        cell_start[current_cell] = i;

    // It detects the end of a cell
    if (i == N - 1 || particle_cell[i + 1] != current_cell)
        cell_end[current_cell] = i;

}

