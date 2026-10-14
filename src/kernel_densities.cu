/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"
#include <math.h>

#include <stdio.h>

// Function to calculate the numerical density of each fluid particle
__global__ void kernel_densities(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 real* __restrict__ dens,
				 real* __restrict__ inv_dens,				 
				 int*  __restrict__ particle_index,
				 int*  __restrict__ cell_start,
				 int*  __restrict__ cell_end,
				 float4* __restrict__ posf) {
  
  // It is faster using pid to label particles in the GPU than i
  // (not at the beginning, but after a while)
  // int i = threadIdx.x + blockIdx.x * blockDim.x;
  // if (i >= N) return;
  int pid = threadIdx.x + blockIdx.x * blockDim.x;
  if (pid >= N) return;
  int i = particle_index[pid];  

  real rijsq;
  int cx_i, cy_i, cz_i;
  int cx, cy, cz;
  real xi     = __ldg(&x[i]);
  real yi     = __ldg(&y[i]);
  real zi     = __ldg(&z[i]);
  real L0     = __ldg(&L[0]);
  real L1     = __ldg(&L[1]);
  real L2     = __ldg(&L[2]);

  //-- Position and box size in float --
  float4 ri_f    = __ldg(&posf[pid]);
  float xi_f     = ri_f.x;
  float yi_f     = ri_f.y;
  float zi_f     = ri_f.z;  
  float L0f      = Lf[0];
  float L1f      = Lf[1];
  float L2f      = Lf[2];
  float half_L0f = 0.5f * L0f;
  float half_L1f = 0.5f * L1f;
  float half_L2f = 0.5f * L2f;
  
  int Ncells0 = __ldg(&Ncells[0]);
  int Ncells1 = __ldg(&Ncells[1]);
  int Ncells2 = __ldg(&Ncells[2]);  
  real dens_i = 0.0;
  
  //The cell of the particle i is calculated
  cx_i = floor(xi / cell_size[0]);
  if (wall == 0) // No wall
    cy_i = floor(yi / cell_size[1]);
  else // With wall
    cy_i = floor((y[i] + wall_width) / cell_size[1]);
  if ( dim == 3 )
    cz_i = floor(zi / cell_size[2]);

  if (dim == 2) //------------------- dim = 2 ----------------
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
	cx = cx_i + dx;
	cy = cy_i + dy;
	
	// periodic boundary conditions
	if (cx < 0)
	  cx = cx + Ncells0;
	if (cx >= Ncells0)
	  cx = cx - Ncells0;
	if (wall == 0) {    //If there are not walls
	  if (cy < 0)
	    cy = cy + Ncells1;
	  if (cy >= Ncells1)	  
	    cy = cy - Ncells1;
	}

	//Neighbour cell
	int neigh_cell = cy * Ncells0 + cx;
	
	// Particles of the cell neigh_cell
	if (cy >= 0 && cy < Ncells1) // This can happen when there are walls
	  for (int k = cell_start[neigh_cell]; k <= cell_end[neigh_cell]; ++k) {
	    int j = particle_index[k]; // index of the neighbour particle

	    //-- Checking in float if they are neighbours --
	    //   It is checked in float because it is much faster --	      
	    float4 pj    = __ldg(&posf[k]);
	    float sx_PBC = 0.f;
	    float sy_PBC = 0.f;	    
	    float dxf    = xi_f - pj.x;
	    float dyf    = yi_f - pj.y;
	    if (dxf >  half_L0f)
	      sx_PBC = 1.f;
	    if (dxf < -half_L0f)
	      sx_PBC = -1.f;
	    if (wall == 0) {
	      if (dyf >  half_L1f)
		sy_PBC = 1.f;
	      if (dyf < -half_L1f)
		sy_PBC = - 1.f;		
	    }
	    dxf -= sx_PBC * L0f;
	    dyf -= sy_PBC * L1f;	    
	    if (dxf*dxf + dyf*dyf >= rcutsq_f)
	      continue; // It jumps to the following particle  (they are not neighbours) 

	    real xij = xi - x[j] - (real)sx_PBC * L0;
	    real yij = yi - y[j] - (real)sy_PBC * L1;	    
	      
	    rijsq = xij*xij + yij*yij;
	    if (rijsq < rcutsq) {
	      real r = sqrt(rijsq);
	      dens_i = dens_i + kernel_W(r);
	    }
	  
	  }
      }
    }
  else //---------------- dim = 3 --------------------------
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
	for (int dz = -1; dz <= 1; ++dz) {	
	  cx = cx_i + dx;
	  cy = cy_i + dy;
	  cz = cz_i + dz;	  
	
	  // periodic boundary conditions
	  if (cx < 0)
	    cx = cx + Ncells0;
	  if (cz < 0)
	    cz = cz + Ncells2;	  
	  if (cx >= Ncells0)
	    cx = cx - Ncells0;
	  if (cz >= Ncells2)	  
	    cz = cz - Ncells2;
	  if (wall == 0) {    //If there are not walls	  
	    if (cy < 0)
	      cy = cy + Ncells1;
	    if (cy >= Ncells1)	  
	      cy = cy - Ncells1;
	  }

	  //Neighbour cell
	  int neigh_cell = cz * Ncells0 * Ncells1 + cy * Ncells0 + cx;

	  // Particles of the cell neigh_cell
	  if (cy >= 0 && cy < Ncells1) // This can happen when there are walls	  
	    for (int k = cell_start[neigh_cell]; k <= cell_end[neigh_cell]; ++k) {
	      int j = particle_index[k]; // index of the neighbour particle

	      //-- Checking in float if they are neighbours --
	      //   It is checked in float because it is much faster --	      
	      float4 pj = __ldg(&posf[k]);
	      float sx_PBC = 0.f;
	      float sy_PBC = 0.f;
	      float sz_PBC = 0.f;
	      float dxf = xi_f - pj.x;
	      float dyf = yi_f - pj.y;
	      float dzf = zi_f - pj.z;
	      if (dxf >  half_L0f)
		sx_PBC = 1.f;
	      if (dxf < -half_L0f)
		sx_PBC = -1.f;
	      if (wall == 0) {
		if (dyf >  half_L1f)
		  sy_PBC = 1.f;
		if (dyf < -half_L1f)
		  sy_PBC = - 1.f;		
	      }
	      if (dzf >  half_L2f)
		sz_PBC = 1.f;
	      if (dzf < -half_L2f)
		sz_PBC = -1.f;
	      dxf -= sx_PBC * L0f;
	      dyf -= sy_PBC * L1f;
	      dzf -= sz_PBC * L2f;	    	      
	      if (dxf*dxf + dyf*dyf + dzf*dzf >= rcutsq_f)
		continue; // it jumps to the following particle (they are not neighbours)
	      
	      // If they are neighbours, the calculation continues
	      real xij = xi - x[j] - (real)sx_PBC * L0;
	      real yij = yi - y[j] - (real)sy_PBC * L1;
	      real zij = zi - z[j] - (real)sz_PBC * L2;
	    
	      rijsq = xij*xij + yij*yij + zij*zij;
	      if (rijsq < rcutsq) {
		real r = sqrt(rijsq);
		dens_i = dens_i + kernel_W(r);
	      }
	    }
	}
      }
    }
  //The calculated data is stored in the GPU  
  dens[i]     = dens_i;
  inv_dens[i] = 1.0 / dens_i;

  // -------------- Code scaling N^2 ----------------------
  // for (int j = 0; j < N; ++j) {

  //   real xij = x[j] - x[i];
  //   real yij = y[j] - y[i];
  //   //--- Periodic boundary conditions ---
  //   xij = xij - round(xij/L[0]);
  //   yij = yij - round(yij/L[1]);		    

  //   if (dim == 2)
  //     rijsq = xij*xij + yij*yij;
  //   else // dim = 3
  //     {
  // 	real zij = z[j] - z[i];
  // 	//--- Periodic boundary conditions ---
  // 	zij = zij - round(zij/L[2]);      
  // 	rijsq = xij*xij + yij*yij + zij*zij;
  //     }

  //   if (rijsq < rcutsq)
  //     {
  // 	real r = sqrt(rijsq);
  // 	dens[i] = dens[i] + kernel_W(r);
  //     }
  // }
  
}
