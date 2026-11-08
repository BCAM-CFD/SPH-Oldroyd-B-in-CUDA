/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include <iostream>
#include <math.h>
#include "cuda_runtime.h"
#include "class_system.h"

// Lubrication force on colloids.
// Function to calculate the velocity of the colloids implicitly, as described in
// Vázquez-Quesada, A., & Ellero, M. (2016).
// Journal of Non-Newtonian Fluid Mechanics, 233, 37-47.
//***** Careful with 2D. The equations used are for 3D. If you want 2D simulations
//      with the correct lubrication correction you should change some things. See
// Bian, X., & Ellero, M. (2014). Computer Physics Communications, 185(1), 53-62.
int class_system::colloids_vel_implicit(real* k_coll_x,
					real* k_coll_y,
					real* k_coll_z,
					real* k_coll_vx,
					real* k_coll_vy,
					real* k_coll_vz,
					real* k_coll_omegax,
					real* k_coll_omegay,
					real* k_coll_omegaz,
					int*  k_coll_cell_start,
					int*  k_coll_cell_end,
					int*  k_coll_index,
					real  eta,
					real  dt) {

  int error_out = 0;

  int  sweep_decrease;
  real error_decrease;
  real error_increase;
  int N_sweep1;
  int N_sweep2;  
  
  // Data is copied to the CPU
  cudaMemcpy(this->coll_x, k_coll_x,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->coll_y, k_coll_y,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)  
    cudaMemcpy(this->coll_z, k_coll_z,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  
  cudaMemcpy(this->coll_vx, k_coll_vx,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->coll_vy, k_coll_vy,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)  
    cudaMemcpy(this->coll_vz, k_coll_vz,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);

  int* coll_cell_start = new int[Ntotal_cells_colloids];
  int* coll_cell_end   = new int[Ntotal_cells_colloids];
  int* coll_index      = new int[N_colloids];
  cudaMemcpy(coll_cell_start, k_coll_cell_start,
  	     Ntotal_cells_colloids * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(coll_cell_end, k_coll_cell_end,
  	     Ntotal_cells_colloids * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(coll_index, k_coll_index,
  	     N_colloids * sizeof(int), cudaMemcpyDeviceToHost);    

  //-- Some needed arrays are allocated --
  real* V0x      = new real[N_colloids];
  real* V0y      = new real[N_colloids];
  real* V1x      = new real[N_colloids];
  real* V1y      = new real[N_colloids];
  real* V_diff_x = new real[N_colloids];
  real* V_diff_y = new real[N_colloids];
  real* V_best_x = new real[N_colloids];
  real* V_best_y = new real[N_colloids];    
  real* V0z;
  real* V1z;
  real* V_diff_z;
  real* V_best_z;    
  if (dim == 3) {
    V0z      = new real[N_colloids];
    V1z      = new real[N_colloids];
    V_diff_z = new real[N_colloids];
    V_best_z = new real[N_colloids];            
  }
  real sum_Vdiff;
  real sum_norm;
  int N_best;

  //-- Initial velocity is storaged --
  for (int i = 0; i < N_colloids; ++i) {
    V0x[i] = coll_vx[i];
    V0y[i] = coll_vy[i];
    if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
      V0z[i] = coll_vz[i];      
  }

  //--- Implicit velocity is computed with N_sweep ---
  error_out = colloids_vel_implicit_iteration(eta, dt, N_sweep, coll_cell_start,
					      coll_cell_end, coll_index);
  if (error_out != 0)
    return 1;
  
  //-- The calculated velocity is storaged --
  for (int i = 0; i < N_colloids; ++i) {
    V1x[i] = coll_vx[i];
    V1y[i] = coll_vy[i];
    if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
      V1z[i] = coll_vz[i];          
  }

  N_sweep1 = N_sweep / 2;
  sweep_decrease = 1; // false
  error_decrease = 0.0; //-- Necessary to start the loop --
  while(N_sweep1 >= 1 && error_decrease < sweep_tol){
    //-- The initial velocity is recovered --
    for (int i = 0; i < N_colloids; ++i) {
      coll_vx[i] =  V0x[i];
      coll_vy[i] =  V0y[i];      
      if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
	coll_vz[i] =  V0z[i];
    }

    //-- Implicit velocity is computed with N_sweep1
    error_out = colloids_vel_implicit_iteration(eta, dt, N_sweep1, coll_cell_start,
						coll_cell_end, coll_index);
    if (error_out != 0)
      return 1;

    //-- Velocity difference is calculated --
    for (int i = 0; i < N_colloids; ++i) {
      V_diff_x[i] = V1x[i] - coll_vx[i];
      V_diff_y[i] = V1y[i] - coll_vy[i];
      if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
	V_diff_z[i] = V1z[i] - coll_vz[i];
    }

    sum_Vdiff = 0.0;
    sum_norm  = 0.0;
    for (int i = 0; i < N_colloids; ++i) {
      sum_Vdiff = sum_Vdiff + V_diff_x[i] * V_diff_x[i] + V_diff_y[i] * V_diff_y[i];
      sum_norm = sum_norm + V1x[i] * V1x[i] + V1y[i] * V1y[i];      
      if (dim == 3) {
	sum_Vdiff = sum_Vdiff + V_diff_z[i] * V_diff_z[i];
	sum_norm  = sum_norm  + V1z[i] * V1z[i];
      }
    }
    
    //-- Eq. (16) from JNNFM --
    error_decrease = sqrt(sum_Vdiff / sum_norm);

    //-- If error is small enough, we have to decrease N_sweep --
    if (error_decrease < sweep_tol) {
      N_best   = N_sweep;      
      N_sweep  = N_sweep1;
      N_sweep1 = N_sweep / 2;
      for (int i = 0; i < N_colloids; ++i) {
	V_best_x[i]   = V1x[i];
	V_best_y[i]   = V1y[i];
	V1x[i]        = coll_vx[i];
	V1y[i]        = coll_vy[i];
	if (dim == 3) {
	  V_best_z[i] = V1z[i];
	  V1z[i]      = coll_vz[i];
	}
      }
      sweep_decrease = 0; // true
    }
  }

  for (int i = 0; i < N_colloids; ++i) {
    coll_vx[i] = V1x[i];
    coll_vy[i] = V1y[i];
  }
  if (dim == 3)
    for (int i = 0; i < N_colloids; ++i)
      coll_vz[i] = V1z[i];
  
  //   If number of sweeps was decreased, there is no need to perform increasing and
  //   save the current sutiable number of sweeps. --
  if (sweep_decrease == 0) {
    N_sweep = N_best;    
    for (int i = 0; i < N_colloids; ++i) {
      coll_vx[i] = V_best_x[i];
      coll_vy[i] = V_best_y[i];
      if (dim == 3)
	coll_vz[i] = V_best_z[i];      	
    }
  }
  else {  // sweep_decrease != 0
    //  -- If number of sweeps was not decreased, we need to increase the number
    //     of sweeps. There may be two reasons that we did not decrease number of sweeps:
    //       1: N_sweep1 = 1 or
    //       2: error_decrease > sweep_tolerance. --
    N_sweep2 = N_sweep * 2;

    if (N_sweep1 >= 1) // set initial error_increase to be error_decrease. --
      error_increase = error_decrease;
    else // set initial error to be an artibrarily big number --
      error_increase = 1.0e2;

    // -- perform number of sweeps until maximum
    // number or error is small enough. --
    while(N_sweep2 <= N_sweep_max && error_increase > sweep_tol) {
      // -- The initial velocity is recovered --
      for (int i = 0; i < N_colloids; ++i) {
	coll_vx[i] =  V0x[i];
	coll_vy[i] =  V0y[i];      
	if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
	  coll_vz[i] =  V0z[i];
      }

      //-- Implicit velocity is computed with N_sweep2 --
      error_out = colloids_vel_implicit_iteration(eta, dt, N_sweep2, coll_cell_start,
						  coll_cell_end, coll_index);
      if (error_out != 0)
	return 1;

      //-- Velocity difference is calculated --
      for (int i = 0; i < N_colloids; ++i) {
	V_diff_x[i] = V1x[i] - coll_vx[i];
	V_diff_y[i] = V1y[i] - coll_vy[i];
	if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
	  V_diff_z[i] = V1z[i] - coll_vz[i];
      }

      sum_Vdiff = 0.0;
      sum_norm  = 0.0;
      for (int i = 0; i < N_colloids; ++i) {
	sum_Vdiff = sum_Vdiff + V_diff_x[i] * V_diff_x[i] + V_diff_y[i] * V_diff_y[i];
	sum_norm = sum_norm + coll_vx[i] * coll_vx[i] + coll_vy[i] * coll_vy[i];      
	if (dim == 3) {
	  sum_Vdiff = sum_Vdiff + V_diff_z[i] * V_diff_z[i];
	  sum_norm  = sum_norm  + coll_vz[i] * coll_vz[i];
	}
      }
      
      //-- Eq. (16) from JNNFM --
      error_increase = sqrt(sum_Vdiff / sum_norm);
      
      if (error_increase > sweep_tol) {
	N_sweep  = N_sweep2;
	N_sweep2 = N_sweep * 2;
	for (int i = 0; i < N_colloids; ++i) {
	  V1x[i] = coll_vx[i];
	  V1y[i] = coll_vy[i];
	  if (dim == 3) // ****** La eficiencia será mejor con el if fuera del bucle ***
	    V1z[i] = coll_vz[i];          
	  }
      }
      else
	N_sweep = N_sweep2;	
      
    }
    
  }

  //-- The result is returned to the GPU --
  cudaMemcpy(k_coll_vx, coll_vx, N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(k_coll_vy, coll_vy, N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  if (dim == 3)
      cudaMemcpy(k_coll_vz, coll_vz, N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  

  //-- Arrays are liberated --
  delete[] V0x;
  delete[] V0y;
  delete[] V1x;
  delete[] V1y;
  delete[] V_diff_x;
  delete[] V_diff_y;
  delete[] V_best_x;
  delete[] V_best_y;  
  delete[] coll_cell_start;
  delete[] coll_cell_end;
  delete[] coll_index;
  if (dim == 3) {
    delete[] V0z;
    delete[] V1z;
    delete[] V_diff_z;
    delete[] V_best_z;    
  }

  return 0;

}
