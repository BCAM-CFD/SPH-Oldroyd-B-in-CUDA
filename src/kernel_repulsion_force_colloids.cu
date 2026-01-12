#include "kernel_functions.h"
#include "config.h"
#include <math.h>

#include <stdio.h>
// Calculation of the repulsion force between colloidal particles
// The repulsive force is given by
//     F^rep = F0 * tau * exp(-tau * s) / (1.0 - exp(-tau * s)
// where s is the adimensional distance between colloidal particles, tau
// is related to the interaction range. Specifically, for a distance
// R * tau^{-1} the force is small but not negiglible. For a distance
// coll_rep_cutoff = 5*R*tau^{-1} we consider it negiglible.
__global__ void kernel_repulsion_force_colloids(real* __restrict__ coll_x,
						real* __restrict__ coll_y,
						real* __restrict__ coll_z,
						real* __restrict__ fx_colloid,
						real* __restrict__ fy_colloid,
						real* __restrict__ fz_colloid,
						int*  __restrict__ coll_index,
						int*  __restrict__ coll_cell_start,
						int*  __restrict__ coll_cell_end) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;

  int cx_i, cy_i, cz_i;
  int cx, cy, cz;
  real rijsq;  
  real rij[3];
  real Famp = F0_rep * tau_rep;

  real coll_xi = __ldg(&coll_x[i]);
  real coll_yi = __ldg(&coll_y[i]);
  real coll_zi = __ldg(&coll_z[i]);
  real L0     = __ldg(&L[0]);
  real L1     = __ldg(&L[1]);
  real L2     = __ldg(&L[2]);
  real half_L0 = 0.5 * L0;
  real half_L1 = 0.5 * L1;
  real half_L2;
  if (dim == 3)
    half_L2 = 0.5 * L2;  
  int Ncells_colloids0 = __ldg(&Ncells_colloids[0]);
  int Ncells_colloids1 = __ldg(&Ncells_colloids[1]);
  int Ncells_colloids2 = __ldg(&Ncells_colloids[2]);  
  
  //The cell of the particle i is calculated
  cx_i = floor(coll_xi / cell_colloids_size[0]);
  cy_i = floor(coll_yi / cell_colloids_size[1]);
  if ( dim == 3 )
    cz_i = floor(coll_zi / cell_colloids_size[2]);

  if (dim == 2) //------------------- dim = 2 ----------------
    for (int dx = -1; dx <= 1; ++dx) 
      for (int dy = -1; dy <= 1; ++dy) {
	cx = cx_i + dx;
	cy = cy_i + dy;
	
	// periodic boundary conditions
	if (cx < 0)
	  cx = cx + Ncells_colloids0;
	if (cx >= Ncells_colloids0)
	  cx = cx - Ncells_colloids0;
	if (wall == 0) {    //If there are not walls
	  if (cy < 0)
	    cy = cy + Ncells_colloids1;
	  if (cy >= Ncells_colloids1)	  
	    cy = cy - Ncells_colloids1;
	}

	//Neighbour cell
	int neigh_cell = cy * Ncells_colloids0 + cx;

	// Particles of the cell neigh_cell
	if (cy >= 0 && cy < Ncells_colloids1) // Different can happen when there are walls
	  for (int k = coll_cell_start[neigh_cell]; k <= coll_cell_end[neigh_cell]; ++k) {
	    
	    int j = coll_index[k]; // index of the neighbour particle
	    
	    if (i == j)
	      continue;
	    
	    rij[0] = coll_xi - coll_x[j];
	    rij[1] = coll_yi - coll_y[j];
	    //--- Periodic boundary conditions ---
	    // This is faster than using round
	    if (rij[0] > half_L0)
	      rij[0] -= L0;
	    if (rij[0] < -half_L0)
	      rij[0] += L0;
	    if (wall == 0) {   //If there are not walls	      	    
	      if (rij[1] > half_L1)
		rij[1] -= L1;
	      if (rij[1] < -half_L1)
		rij[1] += L1;
	    }
	      
	    rijsq = rij[0]*rij[0] + rij[1]*rij[1];

	    // If both particles are too close
	    if (rijsq < rcuton_coll_rep_sq)
	      rijsq = rcuton_coll_rep_sq;

	    if (rijsq <= rcutoff_coll_rep_sq) {
	      real eij[2];
	      real r              = sqrt(rijsq);
	      real r_inv          = 1.0/sqrt(r);
	      real tau_s_neg      = -tau_rep * (r - 2.0 * coll_R)/coll_R;
	      eij[0]              = rij[0]/r;
	      eij[1]              = rij[1]/r;
	      real exp_tau_s_neg  = exp(tau_s_neg);
	      real Frep_mod       = Famp * exp_tau_s_neg / (1.0 - exp_tau_s_neg);
	      fx_colloid[i]       = fx_colloid[i] + Frep_mod * eij[0];
	      fy_colloid[i]       = fy_colloid[i] + Frep_mod * eij[1];
	    }
	  }
      }


  else //---------------- dim = 3 --------------------------
    for (int dx = -1; dx <= 1; ++dx) 
      for (int dy = -1; dy <= 1; ++dy)
	for (int dz = -1; dz <= 1; ++dz) {	
	  cx = cx_i + dx;
	  cy = cy_i + dy;
	  cz = cz_i + dz;	
	
	  // periodic boundary conditions
	  if (cx < 0)
	    cx = cx + Ncells_colloids0;
	  if (cx >= Ncells_colloids0)
	    cx = cx - Ncells_colloids0;
	  if (wall == 0) {    //If there are not walls
	    if (cy < 0)
	      cy = cy + Ncells_colloids1;
	    if (cy >= Ncells_colloids1)	  
	      cy = cy - Ncells_colloids1;
	  }
	  if (cz < 0)
	    cz = cz + Ncells_colloids2;
	  if (cz >= Ncells_colloids2)
	    cz = cz - Ncells_colloids2;	

	  //Neighbour cell
	  int neigh_cell = cz * Ncells_colloids0 * Ncells_colloids1 +
	    cy * Ncells_colloids0 + cx;

	  // Particles of the cell neigh_cell
	  if (cy >= 0 && cy < Ncells_colloids1) // Different can happen when there are walls
	    for (int k = coll_cell_start[neigh_cell]; k <= coll_cell_end[neigh_cell]; ++k) {
	    
	      int j = coll_index[k]; // index of the neighbour particle
	    
	      if (i == j)
		continue;
	    
	      rij[0] = coll_xi - coll_x[j];
	      rij[1] = coll_yi - coll_y[j];
	      rij[2] = coll_zi - coll_z[j];	    
	      //--- Periodic boundary conditions ---
	      // This is faster than using round
	      if (rij[0] > half_L0)
		rij[0] -= L0;
	      if (rij[0] < -half_L0)
		rij[0] += L0;
	      if (wall == 0) {   //If there are not walls	      	    
		if (rij[1] > half_L1)
		  rij[1] -= L1;
		if (rij[1] < -half_L1)
		  rij[1] += L1;
	      }
	      if (rij[2] > half_L2)
		rij[2] -= L2;
	      if (rij[2] < -half_L2)
		rij[2] += L2;	    
	      
	      rijsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

	      // If both particles are too close
	      if (rijsq < rcuton_coll_rep_sq)
		rijsq = rcuton_coll_rep_sq;

	      if (rijsq <= rcutoff_coll_rep_sq) {
		real eij[3];
		real r              = sqrt(rijsq);
		real r_inv          = 1.0/sqrt(r);
		real tau_s_neg      = -tau_rep * (r - 2.0 * coll_R)/coll_R;
		eij[0]              = rij[0]/r;
		eij[1]              = rij[1]/r;
		eij[2]              = rij[2]/r;	      
		real exp_tau_s_neg  = exp(tau_s_neg);
		real Frep_mod       = Famp * exp_tau_s_neg / (1.0 - exp_tau_s_neg);
		fx_colloid[i]       = fx_colloid[i] + Frep_mod * eij[0];
		fy_colloid[i]       = fy_colloid[i] + Frep_mod * eij[1];
		fz_colloid[i]       = fz_colloid[i] + Frep_mod * eij[2];	      
	      }
	    }
	}
  
}
