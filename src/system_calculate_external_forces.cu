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

// Function to calculate the external force on particles
int class_system::calculate_external_forces(dim3 numBlocks,
					    dim3 threadsPerBlock,
					    real* k_x,
					    real* k_y,
					    real* k_z,
					    real* k_fx,
					    real* k_fy,
					    real* k_fz,
					    real* k_mass,
					    int*  k_type) {
  
  cudaError_t cuda_err;
  kernel_ext_force<<<numBlocks, threadsPerBlock>>>(k_x, k_y, k_z, k_fx, k_fy, k_fz,
						   k_mass, k_type);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess){
    printf("Error in kernel_ext_force: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }  	
  cudaDeviceSynchronize();  // We require the kernel to end to continue
  
  return 0;
  
}
