#include "class_system.h"
#include "config.h"
#include "kernel_functions.h"
#include <stdio.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

#include <stdio.h>

// Search for neighbour cells
int class_system::neighbours_search(dim3  numBlocks,
				    dim3  threadsPerBlock,
				    real* k_x,
				    real* k_y,
				    real* k_z,
				    int*  k_particle_cell,
				    int*  k_particle_index,
				    int*  k_cell_start,
				    int*  k_cell_end,
				    real* k_coll_x,
				    real* k_coll_y,
				    real* k_coll_z,
				    int*  k_coll_cell,
				    int*  k_coll_index,
				    int*  k_coll_cell_start,
				    int*  k_coll_cell_end) {
  cudaError_t cuda_err;

  //------------- Computational particles ----------------------
  
  //--- Neighbour cell list is built ---
  kernel_cell_list<<<numBlocks, threadsPerBlock>>>(k_x, k_y, k_z,
						   k_particle_cell, k_particle_index);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error in kernel_cell_list: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }      
  cudaDeviceSynchronize();  // We require the kernel to end to continue   

  // particle_index and particle_cell are sorted according to particle_cell
  // If we have
  //     particle_index 1 2 3 4 5
  //     particle_cell  0 2 0 1 1
  // After sorting we will have
  //     particle_index 1 3 4 5 2
  //     particle_cell  0 0 1 1 2    
  thrust::device_ptr<int> d_cell(k_particle_cell);
  thrust::device_ptr<int> d_index(k_particle_index);
  thrust::sort_by_key(d_cell, d_cell + this->N, d_index);

  //-- Cell indices for beginning and end of cells are calculated --
  kernel_cell_ranges<<<numBlocks, threadsPerBlock>>>(k_particle_cell, k_cell_start, k_cell_end);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error in kernel_cell_ranges: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }      
  cudaDeviceSynchronize();  // We require the kernel to end to continue

  //------------- Colloidal particles ----------------------
  if (N_colloids >= 2) {
      //--- Neighbour cell list for colloids is built ---
    kernel_cell_colloids_list<<<numBlocks, threadsPerBlock>>>(k_coll_x, k_coll_y, k_coll_z,
							      k_coll_cell, k_coll_index);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) {
      printf("Error in kernel_cell_colloids_list: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }      
    cudaDeviceSynchronize();  // We require the kernel to end to continue

    // coll_index and coll_cell are sorted according to coll_cell
    // If we have
    //     coll_index 0 1 2 3 4
    //     coll_cell  0 2 0 1 1
    // After sorting we will have
    //     coll_index 0 2 3 4 1
    //     coll_cell  0 0 1 1 2    
    thrust::device_ptr<int> d_cell(k_coll_cell);
    thrust::device_ptr<int> d_index(k_coll_index);
    thrust::sort_by_key(d_cell, d_cell + this->N_colloids, d_index);

    //-- Initialization of coll_cell_start and coll_cell_end --
    kernel_initialize_cell_colloids<<<numBlocks, threadsPerBlock>>>(k_coll_cell_start,
								    k_coll_cell_end);    
    
    //-- Cell indices for beginning and end of cells of colloids are calculated --
    kernel_cell_colloids_ranges<<<numBlocks, threadsPerBlock>>>(k_coll_cell,
								k_coll_cell_start,
								k_coll_cell_end);
    cuda_err = cudaGetLastError();
    if (cuda_err != cudaSuccess) {
      printf("Error in kernel_colloids_cell_ranges: %s\n", cudaGetErrorString(cuda_err));
      return 1;
    }      
    cudaDeviceSynchronize();  // We require the kernel to end to continue    
  }

  return 0;

}
