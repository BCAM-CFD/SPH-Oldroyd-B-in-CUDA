#include "class_system.h"
#include "config.h"
#include "kernel_functions.h"
#include <stdio.h>

// Function to calculate the density of each particle
int class_system::calculate_density(dim3 numBlocks,
				     dim3 threadsPerBlock,
				     real* k_x,
				     real* k_y,
				     real* k_z,
				     real* k_dens,
				     int*  k_particle_index,
				     int*  k_cell_start,
				     int*  k_cell_end)
{

    cudaError_t cuda_err;

    //------ Densities are set to zero -------
    kernel_density_to_zero<<<numBlocks, threadsPerBlock>>>(k_dens);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) {
      printf("Error in kernel_density to zero: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }          
    cudaDeviceSynchronize();  // We require the kernel to end to continue

    //------ Densities are calculated ----------
    kernel_densities<<<numBlocks, threadsPerBlock>>>(k_x, k_y, k_z, k_dens,
						     k_particle_index, k_cell_start, k_cell_end);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) {
      printf("Error in kernel_densities: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }      
    cudaDeviceSynchronize();  // We require the kernel to end to continue    

    return 0;

}
