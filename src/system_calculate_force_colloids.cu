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

// Function to calculate forces and torques on colloids
int class_system::calculate_force_colloids(dim3  numBlocks,
					   dim3  threadsPerBlock,
					   real* k_coll_x,
					   real* k_coll_y,
					   real* k_coll_z,
					   real* k_fx_colloid,
					   real* k_fy_colloid,
					   real* k_fz_colloid,
					   real* k_tx_colloid,
					   real* k_ty_colloid,
					   real* k_tz_colloid,
					   real* k_fx,
					   real* k_fy,
					   real* k_fz,
					   real* k_x_center,
					   real* k_y_center,
					   real* k_z_center,
					   int*  k_colloids_list,
					   int*  k_colloids_start,
					   int*  k_coll_index,
					   int*  k_coll_cell_start,
					   int*  k_coll_cell_end)  {

  cudaError_t cuda_err;

  //---- Colloidal forces and torques are set to zero ----
  kernel_colloid_forces_to_zero<<<numBlocks, threadsPerBlock>>>(k_fx_colloid, k_fy_colloid,
								k_fz_colloid, k_tx_colloid,
								k_ty_colloid, k_tz_colloid);
  cudaDeviceSynchronize();  // We require the kernel to end to continue  

  //---- Forces and torques of the computational particles are collected ----
  kernel_collect_force_colloids<<<numBlocks, threadsPerBlock>>>(k_fx_colloid, k_fy_colloid,
								k_fz_colloid,
								k_tx_colloid, k_ty_colloid,
								k_tz_colloid,
								k_fx, k_fy, k_fz,
								k_x_center, k_y_center,
								k_z_center,
								k_colloids_list,
								k_colloids_start);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess){
    printf("Error in kernel_force_colloids: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }  	
  cudaDeviceSynchronize();  // We require the kernel to end to continue

  //---- Repulsion forces are calculated ----
  if (N_colloids >= 2)  {
    kernel_repulsion_force_colloids<<<numBlocks, threadsPerBlock>>>(k_coll_x, k_coll_y,
								    k_coll_z, k_fx_colloid,
								    k_fy_colloid,
								    k_fz_colloid,
								    k_coll_index,
								    k_coll_cell_start,
								    k_coll_cell_end);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess){
      printf("Error in kernel_repulsion_force_colloids: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }  	
    cudaDeviceSynchronize();  // We require the kernel to end to continue

    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess){
      printf("Error in kernel_magnetic_force_colloids: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }  	
    cudaDeviceSynchronize();  // We require the kernel to end to continue
  }
  
  return 0;
}
