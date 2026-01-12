#include "class_system.h"
#include "config.h"
#include "kernel_functions.h"
#include <stdio.h>

// Function to calculate the initial forces of the fluid particles
int class_system::calculate_initial_forces(dim3  numBlocks,
					   dim3  threadsPerBlock,
					   real* k_x,
					   real* k_y,
					   real* k_z,
					   real* k_vx,
					   real* k_vy,
					   real* k_vz,
					   real* k_fx,
					   real* k_fy,
					   real* k_fz,
					   real* k_press,
					   real* k_dens,
					   real* k_mass,
					   int*  k_particle_index,
					   int*  k_cell_start,
					   int*  k_cell_end,
					   int*  k_type,
					   real* k_fx_wall,
					   real* k_fy_wall,
					   real* k_fz_wall,
					   int*  k_walls_list,
					   int*  k_walls_start,
					   real* k_coll_x,
					   real* k_coll_y,
					   real* k_coll_z,
					   real* k_coll_vx,
					   real* k_coll_vy,
					   real* k_coll_vz,
					   real* k_coll_omegax,
					   real* k_coll_omegay,
					   real* k_coll_omegaz,
					   real* k_fx_colloid,
					   real* k_fy_colloid,
					   real* k_fz_colloid,
					   real* k_tx_colloid,
					   real* k_ty_colloid,
					   real* k_tz_colloid,
					   real* k_x_center,
					   real* k_y_center,
					   real* k_z_center, 
					   int*  k_colloids_list,
					   int*  k_colloids_start,
					   int*  k_coll_index,
					   int*  k_coll_cell_start,
					   int*  k_coll_cell_end) {

  int error;

  //---- Densities are calculated ----
  error = calculate_density(numBlocks, threadsPerBlock, k_x, k_y, k_z, k_dens,
		    k_particle_index, k_cell_start, k_cell_end);
  if (error != 0)
    return 1;

  //---- Pressure is calculated ----
  error = calculate_pressure(numBlocks, threadsPerBlock, k_press, k_dens, k_mass);
  if (error != 0)
    return 1;  

  // Interparticle forces are set to zero
  kernel_forces_to_zero<<<numBlocks, threadsPerBlock>>>(k_fx, k_fy, k_fz);
  cudaDeviceSynchronize();  // We require the kernel to end to continue
  
  // //---- Interparticle Forces are calculated ----
  // error = calculate_forces(numBlocks, threadsPerBlock,
  // 			   k_x, k_y, k_z, k_vx, k_vy, k_vz,
  // 			   k_fx, k_fy, k_fz, k_press, k_dens,
  // 			   k_mass, k_particle_index,
  // 			   k_cell_start, k_cell_end, k_type,
  // 			   k_coll_x, k_coll_y, k_coll_z,
  // 			   k_coll_vx, k_coll_vy, k_coll_vz,
  // 			   k_coll_omegax, k_coll_omegay, k_coll_omegaz);
  // if (error != 0)
  //   return 1;  

  //---- External forces are calculated ----
  if (ext_force_type != 0) {
    error = calculate_external_forces(numBlocks, threadsPerBlock,
				      k_x, k_y, k_z,
				      k_fx, k_fy, k_fz, k_mass, k_type);
    if (error != 0)
      return 1; // End of program
  }

  //---- Forces on walls are calculated ----
  if (wall) {
    error = calculate_force_walls(numBlocks, threadsPerBlock,
				  k_fx_wall, k_fy_wall, k_fz_wall,
				  k_fx, k_fy, k_fz, k_walls_list, k_walls_start);
    if (error != 0)
      return 1; // End of program
  }
  
  //---- Forces and torques on colloids are calculated ----
  if (N_colloids > 0) {
    error = calculate_force_colloids(numBlocks, threadsPerBlock,
				     k_coll_x, k_coll_y, k_coll_z,
				     k_fx_colloid, k_fy_colloid, k_fz_colloid,
				     k_tx_colloid, k_ty_colloid, k_tz_colloid,
				     k_fx, k_fy, k_fz,
				     k_x_center, k_y_center, k_z_center,
				     k_colloids_list, k_colloids_start,
				     k_coll_index, k_coll_cell_start, k_coll_cell_end);
    if (error != 0)
      return 1; // End of program
  }

  return 0;
}
