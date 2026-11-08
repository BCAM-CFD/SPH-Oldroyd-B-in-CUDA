/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include <iostream>
#include <math.h>
#include "class_system.h"

// Lubrication force on colloids.
// Function to calculate one iteration of the velocity of the colloids
// implicitly, as described in
// Vázquez-Quesada, A., & Ellero, M. (2016).
// Journal of Non-Newtonian Fluid Mechanics, 233, 37-47.
// This function is called by colloids_vel_implicit
int class_system::colloids_vel_implicit_iteration(real  eta,
						  real  dt,
						  int   N_sweep_used,
						  int*  coll_cell_start,
						  int*  coll_cell_end,
						  int*  coll_index) {

  real posij[3];
  real rij;
  real eij[3];
  real eij_vij_old;
  real rijsq;
  real h;
  real vxi_old;
  real vyi_old;
  real vzi_old;
  real vxj_old;
  real vyj_old;
  real vzj_old;
  real vxij_old;
  real vyij_old;
  real vzij_old;
  real vxij_new;
  real vyij_new;
  real vzij_new;
  int cx, cy, cz;
  int cx_i, cy_i, cz_i;
  real Aij;
  real Bij;  

  //-- For the periodic boundary conditions ---
  real half_L0 = 0.5 * L[0];
  real half_L1 = 0.5 * L[1];
  real half_L2;
  if (dim == 3)
    half_L2 = 0.5 * L[2];    

  real dt_sweep = dt / (real)N_sweep_used;
  //-- inverse of the reduced mass is computed --
  real massr_inv_dt = 2.0 / coll_mass * dt_sweep;

  //-- Eqs (5) from JNNFM (for equal radius) --
  real fij1 = -6.0 * M_PI * eta * coll_R * coll_R / 4.0 * massr_inv_dt;
  real fij2 = -6.0 * M_PI * eta * coll_R * 9.0 / 40.0 * massr_inv_dt;
  real gij  = -6.0 * M_PI * eta * coll_R / 6.0 * massr_inv_dt;

  // gap limit
  real h_lub_cutoff = rcutoff_coll_lub - 2.0 * coll_R;
  real h_lub_cutoff_inv = 1.0 / h_lub_cutoff;
  real h_tan_lub_cutoff = h_lub_cutoff / 4.0;
  real h_tan_lub_cutoff_inv = 1.0 / h_tan_lub_cutoff;

  // To gain efficiency
  real log_cut     = log(coll_R * h_lub_cutoff_inv);
  real log_cut_tan = log(coll_R * h_tan_lub_cutoff_inv);

  if (dim == 2)
    for (int t = 0; t < N_sweep_used; ++t)
      for (int i = 0; i < N_colloids; ++i) {
	//The cell of the particle i is calculated
	cx_i = floor(coll_x[i] / cell_colloids_size[0]);
	cy_i = floor(coll_y[i] / cell_colloids_size[1]);
	for (int dx = -1; dx <= 1; ++dx) 
	  for (int dy = -1; dy <= 1; ++dy) {
	    cx = cx_i + dx;
	    cy = cy_i + dy;
	
	    // periodic boundary conditions
	    if (cx < 0)
	      cx = cx + Ncells_colloids[0];
	    if (cx >= Ncells_colloids[0])
	      cx = cx - Ncells_colloids[0];
	    if (wall == 0) {    //If there are not walls
	      if (cy < 0)
	    	cy = cy + Ncells_colloids[1];
	      if (cy >= Ncells_colloids[1])	  
	    	cy = cy - Ncells_colloids[1];
	    }

	    //Neighbour cell
	    int neigh_cell = cy * Ncells_colloids[0] + cx;

	    // Particles of the cell neigh_cell
	    if (cy >= 0 && cy < Ncells_colloids[1]) // Different can happen when there are walls
	      for (int k = coll_cell_start[neigh_cell]; k <= coll_cell_end[neigh_cell]; ++k) {
		
		int j = coll_index[k]; // index of the neighbour particle
		
	    	if (j <= i) // To avoid looping over i,j and j,i
	    	  continue;
		
	    	posij[0] = coll_x[i] - coll_x[j];
	    	posij[1] = coll_y[i] - coll_y[j];
	    	//--- Periodic boundary conditions ---
	    	// This is faster than using round
	    	if (posij[0] > half_L0)
	    	  posij[0] -= L[0];
	    	if (posij[0] < -half_L0)
	    	  posij[0] += L[0];
	    	if (wall == 0) {   //If there are not walls	      	    
	    	  if (posij[1] > half_L1)
	    	    posij[1] -= L[1];
	    	  if (posij[1] < -half_L1)
	    	    posij[1] += L[1];
		}
	      
		rijsq = posij[0]*posij[0] + posij[1]*posij[1];

		if (rijsq < rcutoff_coll_lub_sq) {
		  rij = sqrt(rijsq);
		  //-- The old velocity is storaged --
		  vxi_old = coll_vx[i];
		  vyi_old = coll_vy[i];
		  vxj_old = coll_vx[j];
		  vyj_old = coll_vy[j];

		  vxij_old = vxi_old - vxj_old;
		  vyij_old = vyi_old - vyj_old;

		  eij[0] = posij[0]/rij;
		  eij[1] = posij[1]/rij;

		  eij_vij_old = eij[0] * vxij_old + eij[1] * vyij_old;
		  if (rij < rcuton_coll_lub)
		    rij = rcuton_coll_lub;

		  //-- The gap between particles is calculated --
		  h = rij - 2.0 * coll_R;
		  real h_inv = 1.0 / h;
		  if (h <= 0) {
		    printf("System colloids vel implicit iterarion error: particles %d and %d are penetrating each other\n", i, j);
		    return 1;
		  }

		  //-- Eqs. (13) from JNNFM --
		  real log_val = log(coll_R * h_inv);
		  Aij = fij1 * (h_inv - h_lub_cutoff_inv) +
		    fij2 * (log_val - log_cut);
		  if (h < h_tan_lub_cutoff) 
		    Bij = gij * (log_val - log_cut_tan);
		  else
		    Bij = 0.0;

		  //-- Eqs (14) and (15) from JNNFM --
		  vxij_new = 1.0 / (1.0 - Bij) *
		    (vxij_old + (Aij - Bij) / (1.0 - Aij) *
		     eij_vij_old * eij[0]);
		  vyij_new = 1.0 / (1.0 - Bij) *
		    (vyij_old + (Aij - Bij) / (1.0 - Aij) *
		     eij_vij_old * eij[1]);

		  coll_vx[j] = 0.5 * (vxi_old - vxij_new + vxj_old);
		  coll_vy[j] = 0.5 * (vyi_old - vyij_new + vyj_old);
		  coll_vx[i] = vxij_new + coll_vx[j];
		  coll_vy[i] = vyij_new + coll_vy[j];		  
		  
		}

	      }
	  }
      }
  else // dim == 3
    for (int t = 0; t < N_sweep_used; ++t)
      for (int i = 0; i < N_colloids; ++i) {
	//The cell of the particle i is calculated
	cx_i = floor(coll_x[i] / cell_colloids_size[0]);
	cy_i = floor(coll_y[i] / cell_colloids_size[1]);
	cz_i = floor(coll_z[i] / cell_colloids_size[2]);	
	for (int dx = -1; dx <= 1; ++dx) 
	  for (int dy = -1; dy <= 1; ++dy) 
	    for (int dz = -1; dz <= 1; ++dz) {
	      cx = cx_i + dx;
	      cy = cy_i + dy;
	      cz = cz_i + dz;	      
	      
	      // periodic boundary conditions
	      if (cx < 0)
		cx = cx + Ncells_colloids[0];
	      if (cx >= Ncells_colloids[0])
		cx = cx - Ncells_colloids[0];
	      if (wall == 0) {    //If there are not walls
		if (cy < 0)
		  cy = cy + Ncells_colloids[1];
		if (cy >= Ncells_colloids[1])	  
	    	cy = cy - Ncells_colloids[1];
	      }
	      if (cz < 0)
		cz = cz + Ncells_colloids[2];
	      if (cz >= Ncells_colloids[2])
		cz = cz - Ncells_colloids[2];	      

	    //Neighbour cell
	      int neigh_cell = cz * Ncells_colloids[0] * Ncells_colloids[1] +
		cy * Ncells_colloids[0] + cx;	      

	    // Particles of the cell neigh_cell
	    if (cy >= 0 && cy < Ncells_colloids[1]) // Different can happen when there are walls
	      for (int k = coll_cell_start[neigh_cell]; k <= coll_cell_end[neigh_cell]; ++k) {
		
		int j = coll_index[k]; // index of the neighbour particle
		
	    	if (j <= i) // To avoid looping over i,j and j,i
	    	  continue;
		
	    	posij[0] = coll_x[i] - coll_x[j];
	    	posij[1] = coll_y[i] - coll_y[j];
	    	posij[2] = coll_z[i] - coll_z[j];		
	    	//--- Periodic boundary conditions ---
	    	// This is faster than using round
	    	if (posij[0] > half_L0)
	    	  posij[0] -= L[0];
	    	if (posij[0] < -half_L0)
	    	  posij[0] += L[0];
	    	if (wall == 0) {   //If there are not walls	      	    
	    	  if (posij[1] > half_L1)
	    	    posij[1] -= L[1];
	    	  if (posij[1] < -half_L1)
	    	    posij[1] += L[1];
		}
	    	if (posij[2] > half_L2)
	    	  posij[2] -= L[2];
	    	if (posij[2] < -half_L2)
	    	  posij[2] += L[2];		
	      
		rijsq = posij[0]*posij[0] + posij[1]*posij[1] + posij[2]*posij[2];

		if (rijsq < rcutoff_coll_lub_sq) {
		  rij = sqrt(rijsq);
		  //-- The old velocity is storaged --
		  vxi_old = coll_vx[i];
		  vyi_old = coll_vy[i];
		  vzi_old = coll_vz[i];		  
		  vxj_old = coll_vx[j];
		  vyj_old = coll_vy[j];
		  vzj_old = coll_vz[j];		  

		  vxij_old = vxi_old - vxj_old;
		  vyij_old = vyi_old - vyj_old;
		  vzij_old = vzi_old - vzj_old;		  

		  eij[0] = posij[0]/rij;
		  eij[1] = posij[1]/rij;
		  eij[2] = posij[2]/rij;		  

		  eij_vij_old = eij[0] * vxij_old + eij[1] * vyij_old + eij[2] * vzij_old;
		  if (rij < rcuton_coll_lub)
		    rij = rcuton_coll_lub;

		  //-- The gap between particles is calculated --
		  h = rij - 2.0 * coll_R;
		  real h_inv = 1.0 / h;
		  if (h <= 0) {
		    printf("System colloids vel implicit iterarion error: particles %d and %d are penetrating each other\n", i, j);
		    return 1;
		  }

		  //-- Eqs. (13) from JNNFM --
		  real log_val = log(coll_R * h_inv);
		  Aij = fij1 * (h_inv - h_lub_cutoff_inv) +
		    fij2 * (log_val - log_cut);
		  if (h < h_tan_lub_cutoff) 
		    Bij = gij * (log_val - log_cut_tan);
		  else
		    Bij = 0.0;

		  //-- Eqs (14) and (15) from JNNFM --
		  vxij_new = 1.0 / (1.0 - Bij) *
		    (vxij_old + (Aij - Bij) / (1.0 - Aij) *
		     eij_vij_old * eij[0]);
		  vyij_new = 1.0 / (1.0 - Bij) *
		    (vyij_old + (Aij - Bij) / (1.0 - Aij) *
		     eij_vij_old * eij[1]);
		  vzij_new = 1.0 / (1.0 - Bij) *
		    (vzij_old + (Aij - Bij) / (1.0 - Aij) *
		     eij_vij_old * eij[2]);		  

		  coll_vx[j] = 0.5 * (vxi_old - vxij_new + vxj_old);
		  coll_vy[j] = 0.5 * (vyi_old - vyij_new + vyj_old);
		  coll_vz[j] = 0.5 * (vzi_old - vzij_new + vzj_old);		  
		  coll_vx[i] = vxij_new + coll_vx[j];
		  coll_vy[i] = vyij_new + coll_vy[j];
		  coll_vz[i] = vzij_new + coll_vz[j];		  		  
		  
		}

	      }
	  }
      }

  return 0;
}
