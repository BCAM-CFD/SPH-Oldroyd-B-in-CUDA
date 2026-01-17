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
// gradv_mu_nu = dv_nu / dx_mu (checked)
int class_system::calculate_vel_gradient(dim3 numBlocks,
					 dim3 threadsPerBlock,
					 real* k_x,
					 real* k_y,
					 real* k_z,
					 real* k_vx,
					 real* k_vy,
					 real* k_vz,
					 real* k_gradvxx,
					 real* k_gradvxy,
					 real* k_gradvxz,
					 real* k_gradvyx,
					 real* k_gradvyy,
					 real* k_gradvyz,
					 real* k_gradvzx,
					 real* k_gradvzy,
					 real* k_gradvzz,
					 real* k_dens,
					 int*  k_particle_index,
					 int*  k_cell_start,
					 int*  k_cell_end,
					 int*  k_type,
					 real* k_coll_x,
					 real* k_coll_y,
					 real* k_coll_z,
					 real* k_coll_vx,
					 real* k_coll_vy,
					 real* k_coll_vz,
					 real* k_coll_omegax,
					 real* k_coll_omegay,
					 real* k_coll_omegaz) {

  cudaError_t cuda_err;

  //----  Velocity gradient tensors are set to zero ----
  kernel_vel_gradients_to_zero<<<numBlocks, threadsPerBlock>>>(k_gradvxx, k_gradvxy,
		       					       k_gradvxz, k_gradvyx,
		 					       k_gradvyy, k_gradvyz,
							       k_gradvzx, k_gradvzy,
							       k_gradvzz );
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error in kernel_vel_gradients_to_zero: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }            
  cudaDeviceSynchronize();  // We require the kernel to end to continue

  //---- Velocity gradient tensors are calculated ----
  kernel_vel_gradients<<<numBlocks, threadsPerBlock>>>(k_x, k_y, k_z,
						       k_vx, k_vy, k_vz,
						       k_gradvxx, k_gradvxy,
						       k_gradvxz, k_gradvyx,
						       k_gradvyy, k_gradvyz,
						       k_gradvzx, k_gradvzy,
						       k_gradvzz, k_dens,
						       k_particle_index, k_cell_start,
						       k_cell_end, k_type,
						       k_coll_x, k_coll_y, k_coll_z,
						       k_coll_vx, k_coll_vy, k_coll_vz,
						       k_coll_omegax, k_coll_omegay,
						       k_coll_omegaz);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error in kernel_densities: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }      
  cudaDeviceSynchronize();  // We require the kernel to end to continue      

  return 0;

}
