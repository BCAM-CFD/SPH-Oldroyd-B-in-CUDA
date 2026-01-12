#include "class_system.h"
#include "config.h"
#include "kernel_functions.h"
#include <stdio.h>

// Function to calculate forces on walls
int class_system::calculate_force_walls(dim3 numBlocks,
					dim3 threadsPerBlock,
					real* k_fx_wall,
					real* k_fy_wall,
					real* k_fz_wall,
					real* k_fx,
					real* k_fy,
					real* k_fz,
					int*  k_walls_list,
					int*  k_walls_start)  {

  cudaError_t cuda_err;
  kernel_collect_force_walls<<<numBlocks, threadsPerBlock>>>(k_fx_wall, k_fy_wall, k_fz_wall,
							     k_fx, k_fy, k_fz,
							     k_walls_list, k_walls_start);
      cuda_err = cudaGetLastError();
      if (cuda_err != cudaSuccess){
	printf("Error in kernel_force_walls: %s\n", cudaGetErrorString(cuda_err));
	return 1;
      }  	
      cudaDeviceSynchronize();  // We require the kernel to end to continue
      
      return 0;
}
