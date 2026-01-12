#include "class_system.h"
#include "config.h"
#include "kernel_functions.h"
#include <stdio.h>

// Function to calculate the pressure of each particle
int class_system::calculate_pressure(dim3 numBlocks,
				      dim3 threadsPerBlock,
				      real* k_press,
				      real* k_dens,
				      real* k_mass)
{

    cudaError_t cuda_err;
    
    kernel_pressures<<<numBlocks, threadsPerBlock>>>(k_press, k_dens, k_mass);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) {
      printf("Error in kernel_pressures: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }      
    cudaDeviceSynchronize();  // We require the kernel to end to continue

    return 0;

}
