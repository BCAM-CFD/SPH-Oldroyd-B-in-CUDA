#include "kernel_functions.h"
#include "config.h"
#include <math.h>
#include <stdio.h>

// Function to calculate the numerical density of each fluid particle
__global__ void kernel_vel_gradients(real* __restrict__ x,
				     real* __restrict__ y,
				     real* __restrict__ z,
				     real* __restrict__ vx,
				     real* __restrict__ vy,
				     real* __restrict__ vz,
				     real* __restrict__ gradvxx,
				     real* __restrict__ gradvxy,
				     real* __restrict__ gradvxz,
				     real* __restrict__ gradvyx,
				     real* __restrict__ gradvyy,
				     real* __restrict__ gradvyz,
				     real* __restrict__ gradvzx,
				     real* __restrict__ gradvzy,
				     real* __restrict__ gradvzz,
				     real* __restrict__ dens,
				     int*  __restrict__ particle_index,
				     int*  __restrict__ cell_start,
				     int*  __restrict__ cell_end,
				     int*  __restrict__ type,
				     real* __restrict__ coll_x,
				     real* __restrict__ coll_y,
				     real* __restrict__ coll_z,
				     real* __restrict__ coll_vx,
				     real* __restrict__ coll_vy,
				     real* __restrict__ coll_vz,
				     real* __restrict__ coll_omegax,
				     real* __restrict__ coll_omegay,
				     real* __restrict__ coll_omegaz) {
  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  real rijsq;
  int cx_i, cy_i, cz_i;
  int cx, cy, cz;
  real xi     = __ldg(&x[i]);
  real yi     = __ldg(&y[i]);
  real zi     = __ldg(&z[i]);
  real vxi    = __ldg(&vx[i]);
  real vyi    = __ldg(&vy[i]);
  real vzi    = __ldg(&vz[i]);
  real dens_i = __ldg(&dens[i]);
  real dens_i_inv = 1.0 / dens_i;
  int type_i = type[i];
  real L0     = __ldg(&L[0]);
  real L1     = __ldg(&L[1]);
  real L2     = __ldg(&L[2]);
  real half_L0 = 0.5 * L0;
  real half_L1 = 0.5 * L1;
  real half_L2;  
  if (dim == 3)
    half_L2 = 0.5 * L2;    
  int Ncells0 = __ldg(&Ncells[0]);
  int Ncells1 = __ldg(&Ncells[1]);
  int Ncells2 = __ldg(&Ncells[2]);  
  real gradvxx_i = 0.0;
  real gradvxy_i = 0.0;
  real gradvyx_i = 0.0;
  real gradvyy_i = 0.0;
  real gradvxz_i;
  real gradvyz_i;
  real gradvzx_i;
  real gradvzy_i;
  real gradvzz_i;
  if (dim == 3) {
    gradvxz_i = 0.0;
    gradvyz_i = 0.0;
    gradvzx_i = 0.0;
    gradvzy_i = 0.0;
    gradvzz_i = 0.0;
  }
  
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

	    if (i == j)
	      continue;
	    
	    if (type_i == type[j] && type_i != 0) // Both are not from fluid
	      continue;	  
	    
	    real xij = xi - x[j];
	    real yij = yi - y[j];
	    //--- Periodic boundary conditions ---
	    // This is faster than using round	    
	    if (xij > half_L0)
	      xij -= L0;
	    if (xij < -half_L0)
	      xij += L0;
	    if (wall == 0) {   //If there are not walls	      	    
	      if (yij > half_L1)
		yij -= L1;
	      if (yij < -half_L1)
		yij += L1;
	    }
	      
	    rijsq = xij*xij + yij*yij;
	    if (rijsq < rcutsq) {
	      real eij[2];
	      real vij[2];
	      real dist_i;
	      real dist_j;
	      real beta;
	      real r = sqrt(rijsq);
	      real r_inv = 1.0/r;
	      eij[0] = xij * r_inv;
	      eij[1] = yij * r_inv;	      
	      real gradW = kernel_gradW(r);
	      //----------- i fluid - j fluid -----------------
	      if (type_i == type[j]) { 	      
		vij[0] = vxi - vx[j];
		vij[1] = vyi - vy[j];
	      }
	      //----------------------  i fluid - j wall or colloid-------------------
	      else if (type_i == 0) {
		if (type[j] == 1) { //  -------- j bottom wall -------
		  //-- Morris boundary conditions --
		  dist_i = yi - y_bottom;
		  dist_j = y_bottom - y[j];
		  beta = 1.0 + dist_j/dist_i;
		  if (beta > beta_max)
		    beta = beta_max;		
		  if (dist_i > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (vxi - V_bottom);
		    vij[1] = beta * vyi;
		  }
		  else { //To avoid weird behaviors (Maybe a bounce back would be better?)
		    vij[0] = vxi - V_bottom;
		    vij[1] = vyi;		      
		  }
		}
		else if (type[j] == 2) { // ----------- j top wall  -----------
		  //-- Morris boundary conditions --
		  dist_i = y_top - yi;
		  dist_j = y[j] - y_top;
		  beta = 1.0 + dist_j/dist_i;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_i > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (vxi - V_top);
		    vij[1] = beta * vyi;
		  }
		  else { //To avoid weird behaviors (Maybe a bounce back would be better?)
		    vij[0] = vxi - V_top;
		    vij[1] = vyi;		      
		  }
		}
		else { // ---------- j colloid (type[j] > 2) -----------
		  int coll_part = type[j] - 3;  // Colloidal particle id
		  // Distance particle-colloid center is calculated
		  real ri_coll[2];
		  ri_coll[0] = xi - coll_x[coll_part];
		  ri_coll[1] = yi - coll_y[coll_part];
		  //-- Periodic boundary conditions --
		  if (ri_coll[0] > half_L0)
		    ri_coll[0] -= L0;
		  if (ri_coll[0] < -half_L0)
		    ri_coll[0] += L0;
		  if (wall == 0) {   //If there are not walls	      	    
		    if (ri_coll[1] > half_L1)
		      ri_coll[1] -= L1;
		    if (ri_coll[1] < -half_L1)
		      ri_coll[1] += L1;
		  }
		  real dist_center  = sqrt(ri_coll[0]*ri_coll[0] + ri_coll[1]*ri_coll[1]);
		  real dist_i_surface = dist_center - coll_R;
		  if (dist_i_surface < 0) { //-- i is inside the colloid --
		    // To avoid weird behaviors (maybe a bounce back is better)
		    vij[0] = vxi - coll_vx[coll_part];
		    vij[1] = vyi - coll_vy[coll_part];
		    //**** Rotation should be considered in this case ****
		  }
		  else { //-- Morris boundary conditions --
		    // Position of the surface respect to the colloid center
		    real x_surface = ri_coll[0]/dist_center * coll_R;
		    real y_surface = ri_coll[1]/dist_center * coll_R;
		    // Velocity on the surface
		    real omegaz        = coll_omegaz[coll_part];
		    real vrotx_surface = -omegaz * y_surface;
		    real vroty_surface =  omegaz * x_surface;
		    // Colloid surface vector inwards
		    real surf_vector[2];
		    //---- Possibility 1 ----
		    // dist_j_center should be calculated previously (with coll_x, etc)
		    // real dist_j_surface = coll_R - dist_j_center;
		    //---- Possibility 2 ----
		    surf_vector[0] = ri_coll[0] / dist_center;
		    surf_vector[1] = ri_coll[1] / dist_center;		  
		    real dist_j_surface = (xij * surf_vector[0] + yij * surf_vector[1]) -
		      dist_i_surface;
		    beta = 1.0 + dist_j_surface / dist_i_surface;
		    if (beta > beta_max)
		      beta = beta_max;				  
		    vij[0] = beta * (vxi - coll_vx[coll_part] - vrotx_surface);
		    vij[1] = beta * (vyi - coll_vy[coll_part] - vroty_surface);
		  }
		}
	      }
	      //----------------------  i wall or colloid - j fluid -------------------	 
	      else { // type[j] == 0 (the case i & j are both not fluid was discarded before)
		if (type_i == 1) {  // ---------- i bottom wall ----------------
		  dist_i = y_bottom - yi;
		  dist_j = y[j] - y_bottom;
		  beta = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_bottom - vx[j]);
		    vij[1] = -beta * vy[j];
		  }
		  else { //To avoid weird behaviors (Maybe a bounce back would be better?
		    vij[0] = V_bottom - vx[j];
		    vij[1] = -vy[j];		      
		  }
		}
		else if (type_i == 2) { // ---------- i top wall ----------------
		  dist_i = yi - y_top;
		  dist_j = y_top - y[j];
		  beta = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_top - vx[j]);
		    vij[1] = -beta * vy[j];
		  }
		  else { //To avoid weird behaviors (Maybe a bounce back would be better?
		    vij[0] = V_top - vx[j];
		    vij[1] = -vy[j];		      
		  }
		}
		else {  // ------------- i colloid (type[j] > 2 ------------------ 
		  int coll_part = type_i - 3;  // Colloidal particle id
		  // Distance particle-colloid center is calculated
		  real rj_coll[2];
		  rj_coll[0] = x[j] - coll_x[coll_part];
		  rj_coll[1] = y[j] - coll_y[coll_part];
		  //-- Periodic boundary conditions --
		  if (rj_coll[0] > half_L0)
		    rj_coll[0] -= L0;
		  if (rj_coll[0] < -half_L0)
		    rj_coll[0] += L0;
		  if (wall == 0) {   //If there are not walls	      	    
		    if (rj_coll[1] > half_L1)
		      rj_coll[1] -= L1;
		    if (rj_coll[1] < -half_L1)
		      rj_coll[1] += L1;
		  }
		  real dist_center  = sqrt(rj_coll[0]*rj_coll[0] + rj_coll[1]*rj_coll[1]);
		  real dist_j_surface = dist_center - coll_R;
		  if (dist_j_surface < 0) { //-- j is inside the colloid --
		    // To avoid weird behaviors (maybe a bounce back is better)
		    vij[0] = coll_vx[coll_part] - vx[j] ;
		    vij[1] = coll_vy[coll_part] - vy[j] ;
		    //**** Rotation should be considered in this case ****
		  }
		  else { //-- Morris boundary conditions --
		    // Position of the surface respect to the colloid center
		    real x_surface = rj_coll[0]/dist_center * coll_R;
		    real y_surface = rj_coll[1]/dist_center * coll_R;
		    // Velocity on the surface
		    real omegaz        = coll_omegaz[coll_part];
		    real vrotx_surface = -omegaz * y_surface;
		    real vroty_surface =  omegaz * x_surface;
		    // Colloid surface vector inwards
		    real surf_vector[2];
		    //---- Possibility 1 ----
		    // dist_j_center should be calculated previously (with coll_x, etc)
		    // real dist_j_surface = coll_R - dist_j_center; 
		    //---- Possibility 2 ----
		    surf_vector[0] = rj_coll[0] / dist_center;
		    surf_vector[1] = rj_coll[1] / dist_center;		  
		    real dist_i_surface = (-xij * surf_vector[0] - yij * surf_vector[1])
		      - dist_j_surface;
		    beta = 1.0 + dist_i_surface / dist_j_surface;
		    if (beta > beta_max)
		      beta = beta_max;				  
		    vij[0] = beta * (coll_vx[coll_part] + vrotx_surface - vx[j]);
		    vij[1] = beta * (coll_vy[coll_part] + vroty_surface - vy[j]);
		  }
		}
	      }
	      
	      gradvxx_i += -dens_i_inv * gradW * eij[0] * vij[0];
	      gradvxy_i += -dens_i_inv * gradW * eij[0] * vij[1];
	      gradvyx_i += -dens_i_inv * gradW * eij[1] * vij[0];
	      gradvyy_i += -dens_i_inv * gradW * eij[1] * vij[1];	      	      
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
	      real xij = xi - x[j];
	      real yij = yi - y[j];
	      real zij = zi - z[j];
	      //--- Periodic boundary conditions ---
	      // This is faster than using round	    
	      if (xij > half_L0)
		xij -= L0;
	      if (xij < -half_L0)
		xij += L0;
	      if (wall == 0) {   //If there are not walls	      	    
		if (yij > half_L1)
		  yij -= L1;
		if (yij < -half_L1)
		  yij += L1;
	      }
	      if (zij > half_L2)
		zij -= L2;
	      if (zij < -half_L2)
		zij += L2;	    	  	    	    
	    
	      rijsq = xij*xij + yij*yij + zij*zij;
	      if (rijsq < rcutsq) {
		real eij[3];
		real vij[3];
		real dist_i;
		real dist_j;
		real beta;
		real r = sqrt(rijsq);
		real r_inv = 1.0/r;
		eij[0] = xij * r_inv;
		eij[1] = yij * r_inv;
		eij[2] = zij * r_inv;	      		
		real gradW = kernel_gradW(r);
		//----------- i fluid - j fluid -----------------
		if (type_i == type[j]) { 	      
		  vij[0] = vxi - vx[j];
		  vij[1] = vyi - vy[j];
		  vij[2] = vzi - vz[j];		  
		}
		//----------------------  i fluid - j wall or colloid-------------------
		else if (type_i == 0) {
		  if (type[j] == 1) { //  -------- j bottom wall -------
		    //-- Morris boundary conditions --
		    dist_i = yi - y_bottom;
		    dist_j = y_bottom - y[j];
		    beta = 1.0 + dist_j/dist_i;
		    if (beta > beta_max)
		      beta = beta_max;		
		    if (dist_i > 0) {
		      // The wall moves in the x direction
		      vij[0] = beta * (vxi - V_bottom);
		      vij[1] = beta * vyi;
		      vij[2] = beta * vzi;		      
		    }
		    else { //To avoid weird behaviors (Maybe a bounce back would be better?)
		      vij[0] = vxi - V_bottom;
		      vij[1] = vyi;
		      vij[2] = vzi;		      		      
		    }
		  }
		  else if (type[j] == 2) { // ----------- j top wall  -----------
		    //-- Morris boundary conditions --
		    dist_i = y_top - yi;
		    dist_j = y[j] - y_top;
		    beta = 1.0 + dist_j/dist_i;
		    if (beta > beta_max)
		      beta = beta_max;				
		    if (dist_i > 0) {
		      // The wall moves in the x direction
		      vij[0] = beta * (vxi - V_top);
		      vij[1] = beta * vyi;
		      vij[2] = beta * vzi;		      
		    }
		    else { //To avoid weird behaviors (Maybe a bounce back would be better?)
		      vij[0] = vxi - V_top;
		      vij[1] = vyi;
		      vij[2] = vzi;		      		      
		    }
		  } 
		else { // ---------- j colloid (type[j] > 2) -----------
		  int coll_part = type[j] - 3;  // Colloidal particle id
		  // Distance particle-colloid center is calculated
		  real ri_coll[2];
		  ri_coll[0] = xi - coll_x[coll_part];
		  ri_coll[1] = yi - coll_y[coll_part];
		  ri_coll[2] = zi - coll_z[coll_part];		  
		  //-- Periodic boundary conditions --
		  if (ri_coll[0] > half_L0)
		    ri_coll[0] -= L0;
		  if (ri_coll[0] < -half_L0)
		    ri_coll[0] += L0;
		  if (wall == 0) {   //If there are not walls	      	    
		    if (ri_coll[1] > half_L1)
		      ri_coll[1] -= L1;
		    if (ri_coll[1] < -half_L1)
		      ri_coll[1] += L1;
		  }
		  if (ri_coll[2] > half_L2)
		    ri_coll[2] -= L2;
		  if (ri_coll[2] < -half_L2)
		    ri_coll[2] += L2;		  
		  real dist_center  = sqrt(ri_coll[0] * ri_coll[0] + 
					   ri_coll[1] * ri_coll[1] + 
					   ri_coll[2] * ri_coll[2]);
		  real dist_i_surface = dist_center - coll_R;
		  if (dist_i_surface < 0) { //-- i is inside the colloid --
		    // To avoid weird behaviors (maybe a bounce back is better)
		    vij[0] = vxi - coll_vx[coll_part];
		    vij[1] = vyi - coll_vy[coll_part];
		    vij[2] = vzi - coll_vz[coll_part];		    
		    //**** Rotation should be considered in this case ****
		  }
		  else { //-- Morris boundary conditions --
		    // Position of the surface respect to the colloid center
		    real x_surface = ri_coll[0]/dist_center * coll_R;
		    real y_surface = ri_coll[1]/dist_center * coll_R;
		    real z_surface = ri_coll[2]/dist_center * coll_R;		    
		    // Velocity on the surface
		    real omegax = coll_omegax[coll_part];	
		    real omegay = coll_omegay[coll_part];
		    real omegaz = coll_omegaz[coll_part];
		    real vrotx_surface = omegay * z_surface - omegaz * y_surface;
		    real vroty_surface = omegaz * x_surface - omegax * z_surface;
		    real vrotz_surface = omegax * y_surface - omegay * x_surface;	  
		    // Colloid surface vector inwards
		    real surf_vector[3];
		    //---- Possibility 1 ----
		    // dist_j_center should be calculated previously (with coll_x, etc)
		    // real dist_j_surface = coll_R - dist_j_center;
		    //---- Possibility 2 ----
		    surf_vector[0] = ri_coll[0] / dist_center;
		    surf_vector[1] = ri_coll[1] / dist_center;
		    surf_vector[2] = ri_coll[2] / dist_center;
		    real dist_j_surface = (xij * surf_vector[0] +
					   yij * surf_vector[1] +
					   zij * surf_vector[2]) -  dist_i_surface;
		    beta = 1.0 + dist_j_surface / dist_i_surface;
		    if (beta > beta_max)
		      beta = beta_max;				  
		    vij[0] = beta * (vxi - coll_vx[coll_part] - vrotx_surface);
		    vij[1] = beta * (vyi - coll_vy[coll_part] - vroty_surface);
		    vij[2] = beta * (vzi - coll_vz[coll_part] - vrotz_surface);
		  }
		}
	      }
	      //----------------------  i wall or colloid - j fluid -------------------	 
	      else { // type[j] == 0 (the case i & j are both not fluid was discarded before)
		if (type_i == 1) {  // ---------- i bottom wall ----------------
		  dist_i = y_bottom - yi;
		  dist_j = y[j] - y_bottom;
		  beta = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_bottom - vx[j]);
		    vij[1] = -beta * vy[j];
		    vij[2] = -beta * vz[j];		    		    
		  }
		  else { //To avoid weird behaviors (Maybe a bounce back would be better?
		    vij[0] = V_bottom - vx[j];
		    vij[1] = -vy[j];
		    vij[2] = -vz[j];		      		    
		  }
		}
		else if (type_i == 2) { // ---------- i top wall ----------------
		  dist_i = yi - y_top;
		  dist_j = y_top - y[j];
		  beta = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_top - vx[j]);
		    vij[1] = -beta * vy[j];
		    vij[2] = -beta * vz[j];		    
		  }
		  else { //To avoid weird behaviors (Maybe a bounce back would be better?
		    vij[0] = V_top - vx[j];
		    vij[1] = -vy[j];
		    vij[2] = -vz[j];		      		    
		  }
		}
		else {  // ------------- i colloid (type[j] > 2 ------------------ 
		  int coll_part = type_i - 3;  // Colloidal particle id
		  // Distance particle-colloid center is calculated
		  real rj_coll[2];
		  rj_coll[0] = x[j] - coll_x[coll_part];
		  rj_coll[1] = y[j] - coll_y[coll_part];
		  rj_coll[2] = z[j] - coll_z[coll_part];
		  //-- Periodic boundary conditions --
		  if (rj_coll[0] > half_L0)
		    rj_coll[0] -= L0;
		  if (rj_coll[0] < -half_L0)
		    rj_coll[0] += L0;
		  if (wall == 0) {   //If there are not walls	      	    
		    if (rj_coll[1] > half_L1)
		      rj_coll[1] -= L1;
		    if (rj_coll[1] < -half_L1)
		      rj_coll[1] += L1;
		  }
		  if (rj_coll[2] > half_L2)
		    rj_coll[2] -= L2;
		  if (rj_coll[2] < -half_L2)
		    rj_coll[2] += L2;
		  real dist_center  = sqrt(rj_coll[0]*rj_coll[0] +
					   rj_coll[1]*rj_coll[1] +
					   rj_coll[2]*rj_coll[2]);		  
		  real dist_j_surface = dist_center - coll_R;
		  if (dist_j_surface < 0) { //-- j is inside the colloid --
		    // To avoid weird behaviors (maybe a bounce back is better)
		    vij[0] = coll_vx[coll_part] - vx[j] ;
		    vij[1] = coll_vy[coll_part] - vy[j] ;
		    vij[2] = coll_vz[coll_part] - vz[j] ;
		    //**** Rotation should be considered in this case ****
		  }
		  else { //-- Morris boundary conditions --
		    // Position of the surface respect to the colloid center
		    real x_surface = rj_coll[0]/dist_center * coll_R;
		    real y_surface = rj_coll[1]/dist_center * coll_R;
		    real z_surface = rj_coll[2]/dist_center * coll_R;
		    // Velocity on the surface
		    real omegax        = coll_omegax[coll_part];
		    real omegay        = coll_omegay[coll_part];		    
		    real omegaz        = coll_omegaz[coll_part];
		    real vrotx_surface = omegay * z_surface - omegaz * y_surface;
		    real vroty_surface = omegaz * x_surface - omegax * z_surface;
		    real vrotz_surface = omegax * y_surface - omegay * x_surface;
		    
		    // Colloid surface vector inwards
		    real surf_vector[3];
		    //---- Possibility 1 ----
		    // dist_j_center should be calculated previously (with coll_x, etc)
		    // real dist_j_surface = coll_R - dist_j_center; 
		    //---- Possibility 2 ----
		    surf_vector[0] = rj_coll[0] / dist_center;
		    surf_vector[1] = rj_coll[1] / dist_center;
		    surf_vector[2] = rj_coll[2] / dist_center;
		    real dist_i_surface = (- xij * surf_vector[0] 
					   - yij * surf_vector[1]
					   - zij * surf_vector[2]) - dist_j_surface;	 
		    beta = 1.0 + dist_i_surface / dist_j_surface;
		    if (beta > beta_max)
		      beta = beta_max;				  
		    vij[0] = beta * (coll_vx[coll_part] + vrotx_surface - vx[j]);
		    vij[1] = beta * (coll_vy[coll_part] + vroty_surface - vy[j]);
		    vij[2] = beta * (coll_vz[coll_part] + vrotz_surface - vz[j]);
		  }
		}
	      }
	      
	      gradvxx_i += -dens_i_inv * gradW * eij[0] * vij[0];
	      gradvxy_i += -dens_i_inv * gradW * eij[0] * vij[1];
	      gradvxz_i += -dens_i_inv * gradW * eij[0] * vij[2];	      
	      gradvyx_i += -dens_i_inv * gradW * eij[1] * vij[0];
	      gradvyy_i += -dens_i_inv * gradW * eij[1] * vij[1];
	      gradvyz_i += -dens_i_inv * gradW * eij[1] * vij[2];
	      gradvzx_i += -dens_i_inv * gradW * eij[2] * vij[0];
	      gradvzy_i += -dens_i_inv * gradW * eij[2] * vij[1];
	      gradvzz_i += -dens_i_inv * gradW * eij[2] * vij[2];	      
	      
	      }
	  
	    }
	  }
      }
    }

      //The calculated data is stored in the GPU  
      gradvxx[i] = gradvxx_i;
      gradvxy[i] = gradvxy_i;
      gradvyx[i] = gradvyx_i;
      gradvyy[i] = gradvyy_i;
      if (dim == 3) {
	gradvxz[i] = gradvxz_i;
	gradvyz[i] = gradvyz_i;
	gradvzx[i] = gradvzx_i;
	gradvzy[i] = gradvzy_i;
	gradvzz[i] = gradvzz_i;
      }
  
}
  
