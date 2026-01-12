#include "kernel_functions.h"
#include "config.h"

#include <stdio.h>

//Part 1 of the function to move colloids with a velocity Verlet method with lambda = 0.5.
__global__ void kernel_move_colloids_VV_part1(real* __restrict__ coll_x,
					      real* __restrict__ coll_y,
					      real* __restrict__ coll_z,
					      real* __restrict__ coll_vx,
					      real* __restrict__ coll_vy,
					      real* __restrict__ coll_vz,
					      real* __restrict__ coll_omegax,
					      real* __restrict__ coll_omegay,
					      real* __restrict__ coll_omegaz,
					      real* __restrict__ coll_theta,
					      real* __restrict__ coll_quat0,
					      real* __restrict__ coll_quatx,
					      real* __restrict__ coll_quaty,
					      real* __restrict__ coll_quatz,
					      real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid,
					      real* __restrict__ x,
					      real* __restrict__ y,
					      real* __restrict__ z,
					      real* __restrict__ x_center,
					      real* __restrict__ y_center,
					      real* __restrict__ z_center,
					      int*  __restrict__ colloids_list,
					      int*  __restrict__ colloids_start) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;
    // i is the identity of the colloid 

  //--- Some initialization ---
  real half_dt_over_mass = 0.5 * dt / coll_mass;
  real half_dt_over_I    = 0.5 * dt / coll_I;  
  real coll_xi = coll_x[i];
  real coll_yi = coll_y[i];
  real coll_zi;
  if (dim == 3)
    coll_zi = coll_z[i];

  //--- Velocity at t + dt/2  ---
  coll_vx[i] += half_dt_over_mass * fx_colloid[i];
  coll_vy[i] += half_dt_over_mass * fy_colloid[i];
  if (dim == 3)
    coll_vz[i] += half_dt_over_mass * fz_colloid[i];

  //--- Angular Velocity at t + dt/2  ---
  real omega; 
  coll_omegaz[i] += half_dt_over_I * tz_colloid[i];
  if (dim == 2)
    omega = coll_omegaz[i];  
  else {  // dim == 3 
    coll_omegax[i] += half_dt_over_I * tx_colloid[i];
    coll_omegay[i] += half_dt_over_I * ty_colloid[i];
    omega = sqrt(coll_omegax[i] * coll_omegax[i] +
		 coll_omegay[i] * coll_omegay[i] +
		 coll_omegaz[i] * coll_omegaz[i]);    
  }

  //--- Position at t + dt ---
  coll_xi = coll_xi + coll_vx[i] * dt;
  coll_yi = coll_yi + coll_vy[i] * dt;
  if (dim == 3)
    coll_zi = coll_zi + coll_vz[i] * dt;

  //--- Periodic Boundary conditions ---
  if (coll_xi < 0)
    coll_xi = coll_xi + L[0];
  if (coll_xi > L[0])
    coll_xi = coll_xi - L[0];
  if (wall == 0) {   //--- If no wall ---
    if (coll_yi < 0)
      coll_yi = coll_yi + L[1];
    if (coll_yi > L[1])
      coll_yi = coll_yi - L[1];
  }    
  if (dim == 3) {
    if (coll_zi < 0)
      coll_zi = coll_zi + L[2];
    if (coll_zi > L[2])
      coll_zi = coll_zi - L[2];
  }

  // The variables are stored in the GPU
  coll_x[i] = coll_xi;
  coll_y[i] = coll_yi;
  if (dim == 3)
    coll_z[i] = coll_zi;

  //--- The computational particles are rearranged ---
  // Some definitions
  // Angular variables
  real cos_angle;
  real sin_angle;
  real coll_quat_result[4];
  if (dim == 2) {
    real coll_theta_i = coll_theta[i];
    coll_theta_i = coll_theta_i + coll_omegaz[i] * dt;
    sin_angle = sin(coll_theta_i);
    cos_angle = cos(coll_theta_i);
    coll_theta[i] = coll_theta_i;
  }
  else {   // dim == 3
    // Change of angle during this time step
    real coll_dtheta_half = 0.5 * omega * dt;
    // Rotation axis
    real axis[3] = {coll_omegax[i]/omega,
		    coll_omegay[i]/omega,
		    coll_omegaz[i]/omega};
    // We define the partial quaternion
    real sin_dtheta_half = sin(coll_dtheta_half);
    real d_coll_quat[4] = {cos(coll_dtheta_half),
			   sin_dtheta_half * axis[0],
			   sin_dtheta_half * axis[1],
			   sin_dtheta_half * axis[2]};
    // The quaternion of the colloid is updated
    real coll_quat_i[4] = {coll_quat0[i], coll_quatx[i], coll_quaty[i], coll_quatz[i]};
    kernel_quaternion_product(d_coll_quat, coll_quat_i, coll_quat_result);
    // The quaternion is updated
    coll_quat0[i] = coll_quat_result[0];
    coll_quatx[i] = coll_quat_result[1];
    coll_quaty[i] = coll_quat_result[2];
    coll_quatz[i] = coll_quat_result[3];    
  }
  
  // To define the computational particles
  int j_start    = colloids_start[i];
  int j_end;
  if (i < N_colloids - 1) 
    j_end = colloids_start[i + 1] - 1; 
  else
    j_end = Nlist_colloids - 1;
  // The rearrange is done
  for (int j = j_start; j <= j_end; ++j) {
    int part = colloids_list[j];
    // -------- Rotation -------------
    real current_x_center_j;
    real current_y_center_j;
    real current_z_center_j;            
    if (dim == 2) {
      real x0_center_j = x_center[j];
      real y0_center_j = y_center[j];
      current_x_center_j = x0_center_j * cos_angle - y0_center_j * sin_angle;
      current_y_center_j = x0_center_j * sin_angle + y0_center_j * cos_angle;
    }
    else {  //--- dim = 3 ---
      real x0_center_j = x_center[j];
      real y0_center_j = y_center[j];
      real z0_center_j = z_center[j];
      real quat_point[4] = {0, x0_center_j, y0_center_j, z0_center_j};
      real coll_quat_result_inv[4] = {coll_quat_result[0],
				      -coll_quat_result[1],
				      -coll_quat_result[2],
				      -coll_quat_result[3]};
      real quat_qp[4];
      real quat_qpq[4];      
      kernel_quaternion_product(coll_quat_result, quat_point, quat_qp);
      kernel_quaternion_product(quat_qp, coll_quat_result_inv, quat_qpq);
      current_x_center_j = quat_qpq[1];
      current_y_center_j = quat_qpq[2];
      current_z_center_j = quat_qpq[3];
      // The component quat_qpq[0] should be zero and it is not used.
    }
    //-------- Translation -----------
    x[part] = current_x_center_j + coll_x[i];
    y[part] = current_y_center_j + coll_y[i];
    if (dim == 3)
      z[part] = current_z_center_j + coll_z[i];
    // Periodic boundary conditions
    if (x[part] < 0)
      x[part] = x[part] + L[0];
    if (x[part] > L[0])
      x[part] = x[part] - L[0];
    if (wall == 0) {    //--- If no wall ---
      if (y[part] < 0)
	y[part] = y[part] + L[1];
      if (y[part] > L[1])
	y[part] = y[part] - L[1];
    }      
    if (dim == 3) {
      if (z[part] < 0)
	z[part] = z[part] + L[2];
      if (z[part] > L[2])
	z[part] = z[part] - L[2];          
    }
  }

}
