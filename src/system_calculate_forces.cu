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

// Function to calculate force of each particle.
// and density and pressure too.
int class_system::calculate_forces(dim3 numBlocks,
				   dim3 threadsPerBlock,
				   real* k_x,
				   real* k_y,
				   real* k_z,
				   real* k_vx,
				   real* k_vy,
				   real* k_vz,
				   real* k_cxx,
				   real* k_cxy,
				   real* k_cyy,
				   real* k_cxz,
				   real* k_cyz,
				   real* k_czz,				   
				   real* k_fx,
				   real* k_fy,
				   real* k_fz,
				   real* k_dcdt_xx,
				   real* k_dcdt_xy,
				   real* k_dcdt_yy,
				   real* k_dcdt_xz,
				   real* k_dcdt_yz,
				   real* k_dcdt_zz,				   
				   real* k_gradvxx,
				   real* k_gradvxy,
				   real* k_gradvxz,
				   real* k_gradvyx,
				   real* k_gradvyy,
				   real* k_gradvyz,
				   real* k_gradvzx,
				   real* k_gradvzy,
				   real* k_gradvzz,
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

  cudaError_t cuda_err;
  int error;

  //---- Densities are calculated ----
  error = calculate_density(numBlocks, threadsPerBlock, k_x, k_y, k_z, k_dens,
		    k_particle_index, k_cell_start, k_cell_end);
  if (error != 0)
    return 1;

  //--- Velocity gradient tensor is calculated --
  error = calculate_vel_gradient(numBlocks, threadsPerBlock,
				 k_x, k_y, k_z,
				 k_vx, k_vy, k_vz,
				 k_gradvxx, k_gradvxy, k_gradvxz,
				 k_gradvyx, k_gradvyy, k_gradvyz,
				 k_gradvzx, k_gradvzy, k_gradvzz,			  
				 k_dens, k_particle_index, k_cell_start, k_cell_end,
				 k_type, k_coll_x, k_coll_y, k_coll_z,
				 k_coll_vx, k_coll_vy, k_coll_vz,
				 k_coll_omegax, k_coll_omegay, k_coll_omegaz);
  if (error != 0)
    return 1;
  
  //---- Pressure is calculated ----
  error = calculate_pressure(numBlocks, threadsPerBlock, k_press, k_dens, k_mass);
  if (error != 0)
    return 1;  
  
  // Interparticle forces are set to zero
  kernel_forces_to_zero<<<numBlocks, threadsPerBlock>>>(k_fx, k_fy, k_fz);
  cudaDeviceSynchronize();  // We require the kernel to end to continue
    
  // forces are calculated
  kernel_forces<<<numBlocks, threadsPerBlock>>>(k_x, k_y, k_z, k_vx, k_vy, k_vz,
						k_cxx, k_cxy, k_cyy, k_cxz, k_cyz, k_czz,
						k_fx, k_fy, k_fz,
						k_press, k_dens,
						k_mass, k_particle_index,
						k_cell_start, k_cell_end,
						k_type,
						k_coll_x, k_coll_y, k_coll_z,
						k_coll_vx, k_coll_vy, k_coll_vz,
						k_coll_omegax, k_coll_omegay, k_coll_omegaz);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error in kernel_forces: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }    
  cudaDeviceSynchronize();  // We require the kernel to end to continue

  //---- External forces are calculated ----
  if (ext_force_type != 0) {
    error = calculate_external_forces(numBlocks, threadsPerBlock,
				      k_x, k_y, k_z,
				      k_fx, k_fy, k_fz, k_mass, k_type);
    if (error != 0)
      return 1; // End of program
  }

  //----- Derivative of conformation tensor respect to time is calculated -----
  kernel_dcdt<<<numBlocks, threadsPerBlock>>>(k_dcdt_xx, k_dcdt_xy, k_dcdt_yy,
						k_dcdt_xz, k_dcdt_yz, k_dcdt_zz,
						k_cxx, k_cxy, k_cyy, k_cxz, k_cyz, k_czz,
						k_gradvxx, k_gradvxy, k_gradvxz,
						k_gradvyx, k_gradvyy, k_gradvyz,
						k_gradvzx, k_gradvzy, k_gradvzz, k_type);

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
