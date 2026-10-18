/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

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
				    real* k_inv_dens,				    
				    int*  k_particle_index,
				    int*  k_cell_start,
				    int*  k_cell_end,
				    float4* k_posf)
{

    cudaError_t cuda_err;

    //------ Densities are calculated ----------
    kernel_densities<<<numBlocks, threadsPerBlock>>>(k_x, k_y, k_z, k_dens, k_inv_dens,
						     k_particle_index, k_cell_start,
						     k_cell_end, k_posf);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) {
      printf("Error in kernel_densities: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }      
    // cudaDeviceSynchronize();  // We require the kernel to end to continue    

    return 0;

}
