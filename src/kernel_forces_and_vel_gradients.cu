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
// Calculation of the forces and velocity gradients
// between computational particles (Espanol and Revenga, 2003)
__global__ void kernel_forces_and_vel_gradients(real* __restrict__ x,
						real* __restrict__ y,
						real* __restrict__ z,
						real* __restrict__ vx,
						real* __restrict__ vy,
						real* __restrict__ vz,
						real* __restrict__ fx,
						real* __restrict__ fy,
						real* __restrict__ fz,
						real* __restrict__ gradvxx,
						real* __restrict__ gradvxy,
						real* __restrict__ gradvxz,
						real* __restrict__ gradvyx,
						real* __restrict__ gradvyy,
						real* __restrict__ gradvyz,
						real* __restrict__ gradvzx,
						real* __restrict__ gradvzy,
						real* __restrict__ gradvzz,			      
						real* __restrict__ inv_dens,			      
						int*  __restrict__  particle_index,
						int*  __restrict__  cell_start,
						int*  __restrict__  cell_end,
						int*  __restrict__  type,
						real* __restrict__ Pixx_over_dens_sq,
						real* __restrict__ Pixy_over_dens_sq,
						real* __restrict__ Piyy_over_dens_sq,
						real* __restrict__ Pixz_over_dens_sq,
						real* __restrict__ Piyz_over_dens_sq,
						real* __restrict__ Pizz_over_dens_sq,
						real* __restrict__ coll_x,
						real* __restrict__ coll_y,
						real* __restrict__ coll_z,
						real* __restrict__ coll_vx,
						real* __restrict__ coll_vy,
						real* __restrict__ coll_vz,
						real* __restrict__ coll_omegax,
						real* __restrict__ coll_omegay,
						real* __restrict__ coll_omegaz,
						float4* __restrict__ posf) {

  // It is faster using pid to label particles in the GPU than i
  // (not at the beginning, but after a while)
  // int i = threadIdx.x + blockIdx.x * blockDim.x;
  // if (i >= N) return;
  int pid = threadIdx.x + blockIdx.x * blockDim.x;
  if (pid >= N) return;
  int i = particle_index[pid];  

  int cx_i, cy_i, cz_i;
  int cx, cy, cz;
  real rijsq;  
  real rij[3];
  real fxi = 0.0;
  real fyi = 0.0;
  real fzi = 0.0;  
  real inv_dens_i0 = inv_dens[i];
  real Pixx_over_dens_sq_i0 = Pixx_over_dens_sq[i];
  real Pixy_over_dens_sq_i0 = Pixy_over_dens_sq[i];
  real Piyy_over_dens_sq_i0 = Piyy_over_dens_sq[i];
  real Pixz_over_dens_sq_i0;
  real Piyz_over_dens_sq_i0;
  real Pizz_over_dens_sq_i0;
  if (dim == 3) {
    Pixz_over_dens_sq_i0 = Pixz_over_dens_sq[i];
    Piyz_over_dens_sq_i0 = Piyz_over_dens_sq[i];
    Pizz_over_dens_sq_i0 = Pizz_over_dens_sq[i];
  }
  real inv_dens_i;
  real Pixx_over_dens_sq_i;
  real Pixy_over_dens_sq_i;
  real Piyy_over_dens_sq_i;
  real Pixz_over_dens_sq_i;
  real Piyz_over_dens_sq_i;
  real Pizz_over_dens_sq_i;      

  real xi     = __ldg(&x[i]);
  real yi     = __ldg(&y[i]);
  real zi     = __ldg(&z[i]);
  real vxi    = __ldg(&vx[i]);
  real vyi    = __ldg(&vy[i]);
  real vzi    = __ldg(&vz[i]);  
  int type_i  = __ldg(&type[i]);    
  real L0     = __ldg(&L[0]);
  real L1     = __ldg(&L[1]);
  real L2     = __ldg(&L[2]);
  real half_L0 = 0.5 * L0;
  real half_L1 = 0.5 * L1;
  real half_L2;
  if (dim == 3)
    half_L2 = 0.5 * L2;

  //-- Position and box size in float -- (faster than xi_f = (float) xi; etc
  float4 ri_f = __ldg(&posf[pid]);
  float xi_f = ri_f.x;
  float yi_f = ri_f.y;
  float zi_f = ri_f.z;
  float L0f = Lf[0];
  float L1f = Lf[1];
  float L2f = Lf[2];
  float half_L0f = 0.5f * L0f;
  float half_L1f = 0.5f * L1f;
  float half_L2f = 0.5f * L2f;
  
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
	if (cy >= 0 && cy < Ncells1) // Different can happen when there are walls   
	  for (int k = cell_start[neigh_cell]; k <= cell_end[neigh_cell]; ++k) {
	    
	    int j = particle_index[k]; // index of the neighbour particle
	    
	    if (i == j)
	      continue;

	    //-- Checking in float if they are neighbours.
	    //   It is checked in float because it is much faster --
	    float4 pj  = __ldg(&posf[k]);
	    //Following is faster than type_j = type[j], because j is disordered and is read from non constant GPU memory
	    int type_j = __float_as_int(pj.w);
	    if (type_i != 0 && type_j != 0) // Both are not from fluid	      	    	   
	      continue;
	    float sx_PBC = 0.f;
	    float sy_PBC = 0.f;	    
	    float dxf = xi_f - pj.x;
	    float dyf = yi_f - pj.y;
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
	      continue;	// it jumps to the following particle (they are not neighbours)

	    rij[0] = xi - x[j] - (real)sx_PBC * L0;
	    rij[1] = yi - y[j] - (real)sy_PBC * L1;	    
	      
	    rijsq = rij[0]*rij[0] + rij[1]*rij[1];
	  
	    if (rijsq < rcutsq) {
	      real eij[2];	    
	      
	      //-- Remember that, from a line above, type[i] = type[j] only if type[i] = 0 
	      real dist_i;
	      real dist_j;		
	      real vij[2];
	      real beta;
	      real inv_dens_j;
	      real Pixx_over_dens_sq_j;
	      real Pixy_over_dens_sq_j;
	      real Piyy_over_dens_sq_j;	      
	      //----------------------  i fluid - j fluid -------------------	    
	      if (type_i == type_j) { 
		inv_dens_i              = inv_dens_i0;
		inv_dens_j              = inv_dens[j];
		Pixx_over_dens_sq_i     = Pixx_over_dens_sq_i0;
		Pixy_over_dens_sq_i     = Pixy_over_dens_sq_i0;
		Piyy_over_dens_sq_i     = Piyy_over_dens_sq_i0;
		Pixx_over_dens_sq_j     = Pixx_over_dens_sq[j];
		Pixy_over_dens_sq_j     = Pixy_over_dens_sq[j];
		Piyy_over_dens_sq_j     = Piyy_over_dens_sq[j];
		vij[0] = vxi - vx[j];
		vij[1] = vyi - vy[j];		  
	      }
	      //----------------------  i fluid - j wall or colloid-------------------
	      else if (type_i == 0) {	      
		inv_dens_i              = inv_dens_i0;
		inv_dens_j              = inv_dens_i;
		Pixx_over_dens_sq_i     = Pixx_over_dens_sq_i0;
		Pixy_over_dens_sq_i     = Pixy_over_dens_sq_i0;
		Piyy_over_dens_sq_i     = Piyy_over_dens_sq_i0;
		Pixx_over_dens_sq_j     = Pixx_over_dens_sq_i0;
		Pixy_over_dens_sq_j     = Pixy_over_dens_sq_i0;
		Piyy_over_dens_sq_j     = Piyy_over_dens_sq_i0;		
		if (type_j == 1) { //  -------- j bottom wall -------
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
		else if (type_j == 2) { // ----------- j top wall  -----------
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
		  int coll_part = type_j - 3;  // Colloidal particle id
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
		    surf_vector[0]      = ri_coll[0] / dist_center;
		    surf_vector[1]      = ri_coll[1] / dist_center;		  
		    real dist_j_surface = (rij[0] * surf_vector[0] +
					   rij[1] * surf_vector[1]) - dist_i_surface;
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
		inv_dens_j              = inv_dens[j];
		inv_dens_i              = inv_dens_j;
		Pixx_over_dens_sq_j     = Pixx_over_dens_sq[j];
		Pixy_over_dens_sq_j     = Pixy_over_dens_sq[j];
		Piyy_over_dens_sq_j     = Piyy_over_dens_sq[j];
		Pixx_over_dens_sq_i     = Pixx_over_dens_sq_j;
		Pixy_over_dens_sq_i     = Pixy_over_dens_sq_j;
		Piyy_over_dens_sq_i     = Piyy_over_dens_sq_j;			
		if (type_i == 1) {  // ---------- i bottom wall ----------------
		  dist_i = y_bottom - yi;
		  dist_j = y[j] - y_bottom;
		  beta   = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_bottom - vx[j]);
		    vij[1] = -beta * vy[j];
		  }
		  else {   //To avoid weird behaviors (Maybe a bounce back would be better?
		    vij[0] = V_bottom - vx[j];
		    vij[1] = -vy[j];		      
		  }
		}
		else if (type_i == 2) { // ---------- i top wall ----------------
		  dist_i = yi - y_top;
		  dist_j = y_top - y[j];
		  beta   = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_top - vx[j]);
		    vij[1] = -beta * vy[j];
		  }
		  else {   //To avoid weird behaviors (Maybe a bounce back would be better?
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
		    real x_surface     = rj_coll[0]/dist_center * coll_R;
		    real y_surface     = rj_coll[1]/dist_center * coll_R;
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
		    surf_vector[0]      = rj_coll[0] / dist_center;
		    surf_vector[1]      = rj_coll[1] / dist_center;		  
		    real dist_i_surface = (-rij[0] * surf_vector[0] -
					   rij[1] * surf_vector[1]) - dist_j_surface;
		    beta                = 1.0 + dist_i_surface / dist_j_surface;
		    if (beta > beta_max)
		      beta = beta_max;				  
		    vij[0] = beta * (coll_vx[coll_part] + vrotx_surface - vx[j]);
		    vij[1] = beta * (coll_vy[coll_part] + vroty_surface - vy[j]);
		  }
		}
	      }
	      // rsqrt(x) is 1/sqrt(x)
	      // This is faster than calculate
	       // real r = sqrt(rijsq);
	       // real inv_r = 1.0/r;
	      real inv_r = rsqrt(rijsq);
	      real r     = rijsq * inv_r;	      

	       real gradW = kernel_gradW(r);
	    
	      eij[0]           = rij[0] * inv_r;
	      eij[1]           = rij[1] * inv_r;	
	      real eij_dot_vij = eij[0] * vij[0] + eij[1] * vij[1];
	    
	      // Reversible force
	      real aux_xx = -gradW * (Pixx_over_dens_sq_i + Pixx_over_dens_sq_j);
	      real aux_xy = -gradW * (Pixy_over_dens_sq_i + Pixy_over_dens_sq_j);
	      real aux_yy = -gradW * (Piyy_over_dens_sq_i + Piyy_over_dens_sq_j);
	      fxi        += aux_xx * eij[0] + aux_xy * eij[1];
	      fyi        += aux_xy * eij[0] + aux_yy * eij[1];

	      // Irreversible force
	      real aux  = gradW * inv_dens_i *inv_dens_j * inv_r;	      	    
	      real aux1 = a * aux;
	      real aux2 = b * aux * eij_dot_vij;
	      fxi      += aux1 * vij[0] + aux2 * eij[0]; 
	      fyi      += aux1 * vij[1] + aux2 * eij[1];

	      // Velocity gradient tensor
	      // Note that inv_dens_i is not necessarilly equal to inv_dens_i0
	      aux = -inv_dens_i0 * gradW;
	      gradvxx_i += aux * eij[0] * vij[0];
	      gradvxy_i += aux * eij[0] * vij[1];
	      gradvyx_i += aux * eij[1] * vij[0];
	      gradvyy_i += aux * eij[1] * vij[1];	      	      	      
	    
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
	  if (cx >= Ncells0)
	    cx = cx - Ncells0;
	  if (wall == 0) {    //If there are not walls
	    if (cy < 0)
	      cy = cy + Ncells1;
	    if (cy >= Ncells1)	  
	      cy = cy - Ncells1;
	  }
	  if (cz < 0)
	    cz = cz + Ncells2;
	  if (cz >= Ncells2)
	    cz = cz - Ncells2;	  

	  //Neighbour cell
	  int neigh_cell = cz * Ncells0 * Ncells1 + cy * Ncells0 + cx;		  

	// Particles of the cell neigh_cell
	if (cy >= 0 && cy < Ncells1) // Different can happen when there are walls   
	  for (int k = cell_start[neigh_cell]; k <= cell_end[neigh_cell]; ++k) {
	    int j = particle_index[k]; // index of the neighbour particle
	    if (i == j)
	      continue;

	    //-- Checking in float if they are neighbours.
	    //   It is checked in float because it is much faster --
	    float4 pj  = __ldg(&posf[k]);
	    //Following is faster than type_j = type[j], because j is disordered and is read from non constant GPU memory
	    int type_j = __float_as_int(pj.w);
	    if (type_i != 0 && type_j != 0) // Both are not from fluid
	      continue; // it jumps to the following particle 
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
	      continue;	// it jumps to the following particle (they are not neighbours)
	  
	    rij[0] = xi - x[j] - (real)sx_PBC * L0;
	    rij[1] = yi - y[j] - (real)sy_PBC * L1;
	    rij[2] = zi - z[j] - (real)sz_PBC * L2;	    

	    rijsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];	    
	  
	    if (rijsq < rcutsq) {
	      real eij[3];	    
	      
	      //-- Remember that, from a line above, type[i] = type[j] only if type[i] = 0 
	      real dist_i;
	      real dist_j;		
	      real vij[3];
	      real beta;
	      real inv_dens_j;
	      real Pixx_over_dens_sq_j;
	      real Pixy_over_dens_sq_j;
	      real Piyy_over_dens_sq_j;
	      real Pixz_over_dens_sq_j;
	      real Piyz_over_dens_sq_j;
	      real Pizz_over_dens_sq_j;	      	      	      
	      //----------------------  i fluid - j fluid -------------------	    
	      if (type_i == type_j) { 
		inv_dens_i              = inv_dens_i0;
		inv_dens_j              = inv_dens[j];
		Pixx_over_dens_sq_i     = Pixx_over_dens_sq_i0;
		Pixy_over_dens_sq_i     = Pixy_over_dens_sq_i0;
		Piyy_over_dens_sq_i     = Piyy_over_dens_sq_i0;
		Pixz_over_dens_sq_i     = Pixz_over_dens_sq_i0;
		Piyz_over_dens_sq_i     = Piyz_over_dens_sq_i0;
		Pizz_over_dens_sq_i     = Pizz_over_dens_sq_i0;
		Pixx_over_dens_sq_j     = Pixx_over_dens_sq[j];
		Pixy_over_dens_sq_j     = Pixy_over_dens_sq[j];
		Piyy_over_dens_sq_j     = Piyy_over_dens_sq[j];
		Pixz_over_dens_sq_j     = Pixz_over_dens_sq[j];
		Piyz_over_dens_sq_j     = Piyz_over_dens_sq[j];
		Pizz_over_dens_sq_j     = Pizz_over_dens_sq[j];
		vij[0] = vxi - vx[j];
		vij[1] = vyi - vy[j];
		vij[2] = vzi - vz[j];		  		
	      }
	      //----------------------  i fluid - j wall or colloid-------------------
	      else if (type_i == 0) {	      
		inv_dens_i               = inv_dens_i0;
		inv_dens_j               = inv_dens_i;
		Pixx_over_dens_sq_i      = Pixx_over_dens_sq_i0;
		Pixy_over_dens_sq_i      = Pixy_over_dens_sq_i0;
		Piyy_over_dens_sq_i      = Piyy_over_dens_sq_i0;
		Pixz_over_dens_sq_i      = Pixz_over_dens_sq_i0;
		Piyz_over_dens_sq_i      = Piyz_over_dens_sq_i0;
		Pizz_over_dens_sq_i      = Pizz_over_dens_sq_i0;
		Pixx_over_dens_sq_j      = Pixx_over_dens_sq_i0;
		Pixy_over_dens_sq_j      = Pixy_over_dens_sq_i0;
		Piyy_over_dens_sq_j      = Piyy_over_dens_sq_i0;
		Pixz_over_dens_sq_j      = Pixz_over_dens_sq_i0;
		Piyz_over_dens_sq_j      = Piyz_over_dens_sq_i0;
		Pizz_over_dens_sq_j      = Pizz_over_dens_sq_i0;		
		if (type_j == 1) { //  -------- j bottom wall -------
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
		else if (type_j == 2) { // ----------- j top wall  -----------
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
		  int coll_part = type_j - 3;  // Colloidal particle id
		  // Distance particle-colloid center is calculated
		  real ri_coll[3];
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
		    real dist_center  = sqrt(ri_coll[0]*ri_coll[0] +
					     ri_coll[1]*ri_coll[1] +
					     ri_coll[2]*ri_coll[2]);		    
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
		    real x_surface     = ri_coll[0]/dist_center * coll_R;
		    real y_surface     = ri_coll[1]/dist_center * coll_R;
		    real z_surface     = ri_coll[2]/dist_center * coll_R;		    
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
		    surf_vector[0]      = ri_coll[0] / dist_center;
		    surf_vector[1]      = ri_coll[1] / dist_center;
		    surf_vector[2]      = ri_coll[2] / dist_center;		    
		    real dist_j_surface = (rij[0] * surf_vector[0] +
					   rij[1] * surf_vector[1] +
					   rij[2] * surf_vector[2]) - dist_i_surface;
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
		inv_dens_j              = inv_dens[j];
		inv_dens_i              = inv_dens_j;
		Pixx_over_dens_sq_j     = Pixx_over_dens_sq[j];
		Pixy_over_dens_sq_j     = Pixy_over_dens_sq[j];
		Piyy_over_dens_sq_j     = Piyy_over_dens_sq[j];
		Pixz_over_dens_sq_j     = Pixz_over_dens_sq[j];
		Piyz_over_dens_sq_j     = Piyz_over_dens_sq[j];
		Pizz_over_dens_sq_j     = Pizz_over_dens_sq[j];
		Pixx_over_dens_sq_i     = Pixx_over_dens_sq_j;
		Pixy_over_dens_sq_i     = Pixy_over_dens_sq_j;
		Piyy_over_dens_sq_i     = Piyy_over_dens_sq_j;
		Pixz_over_dens_sq_i     = Pixz_over_dens_sq_j;
		Piyz_over_dens_sq_i     = Piyz_over_dens_sq_j;
		Pizz_over_dens_sq_i     = Pizz_over_dens_sq_j;
		if (type_i == 1) {  // ---------- i bottom wall ----------------
		  dist_i = y_bottom - yi;
		  dist_j = y[j] - y_bottom;
		  beta   = 1.0 + dist_i/dist_j;
		  if (beta > beta_max)
		    beta = beta_max;				
		  if (dist_j > 0) {
		    // The wall moves in the x direction
		    vij[0] = beta * (V_bottom - vx[j]);
		    vij[1] = -beta * vy[j];
		    vij[2] = -beta * vz[j];		    
		  }
		  else {   //To avoid weird behaviors (Maybe a bounce back would be better?
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
		  else {  //To avoid weird behaviors (Maybe a bounce back would be better?
		    vij[0] = V_top - vx[j];
		    vij[1] = -vy[j];
		    vij[2] = -vz[j];		      		    
		  }
		}
		else {  // ------------- i colloid (type[j] > 2 ------------------ 
		  int coll_part = type_i - 3;  // Colloidal particle id
		  // Distance particle-colloid center is calculated
		  real rj_coll[3];
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
		    real dist_i_surface = (- rij[0] * surf_vector[0] 
					   - rij[1] * surf_vector[1]
					   - rij[2] * surf_vector[2]) - dist_j_surface;	 
		    beta = 1.0 + dist_i_surface / dist_j_surface;
		    if (beta > beta_max)
		      beta = beta_max;				  
		    vij[0] = beta * (coll_vx[coll_part] + vrotx_surface - vx[j]);
		    vij[1] = beta * (coll_vy[coll_part] + vroty_surface - vy[j]);
		    vij[2] = beta * (coll_vz[coll_part] + vrotz_surface - vz[j]);
		  }
		}
	      }
	      // rsqrt(x) is 1/sqrt(x)
	      // This is faster than calculate
	       // real r = sqrt(rijsq);
	       // real inv_r = 1.0/r;
	      real inv_r = rsqrt(rijsq);
	      real r     = rijsq * inv_r;
	      
	      real gradW = kernel_gradW(r);
	    
	      eij[0] = rij[0] * inv_r;
	      eij[1] = rij[1] * inv_r;
	      eij[2] = rij[2] * inv_r;
	      real eij_dot_vij = eij[0] * vij[0] + eij[1] * vij[1] + eij[2] * vij[2];

	      // Reversible force
	      real aux_xx = -gradW * (Pixx_over_dens_sq_i + Pixx_over_dens_sq_j);
	      real aux_xy = -gradW * (Pixy_over_dens_sq_i + Pixy_over_dens_sq_j);
	      real aux_yy = -gradW * (Piyy_over_dens_sq_i + Piyy_over_dens_sq_j);
	      real aux_xz = -gradW * (Pixz_over_dens_sq_i + Pixz_over_dens_sq_j);
	      real aux_yz = -gradW * (Piyz_over_dens_sq_i + Piyz_over_dens_sq_j);
	      real aux_zz = -gradW * (Pizz_over_dens_sq_i + Pizz_over_dens_sq_j);

	      fxi += aux_xx * eij[0] + aux_xy * eij[1] + aux_xz * eij[2];
	      fyi += aux_xy * eij[0] + aux_yy * eij[1] + aux_yz * eij[2];
	      fzi += aux_xz * eij[0] + aux_yz * eij[1] + aux_zz * eij[2];	      

	      // Irreversible force
	      real aux = gradW * inv_dens_i *inv_dens_j * inv_r;	      	    
	      real aux1 = a * aux;
	      real aux2 = b * aux * eij_dot_vij;
	      fxi += aux1 * vij[0] + aux2 * eij[0]; 
	      fyi += aux1 * vij[1] + aux2 * eij[1];
	      fzi += aux1 * vij[2] + aux2 * eij[2];	      

	      // Velocity gradient tensor
	      // Note that inv_dens_i is not necessarilly equal to inv_dens_i0
	      aux = -inv_dens_i0 * gradW;
	      gradvxx_i += aux * eij[0] * vij[0];
	      gradvxy_i += aux * eij[0] * vij[1];
	      gradvxz_i += aux * eij[0] * vij[2];	      
	      gradvyx_i += aux * eij[1] * vij[0];
	      gradvyy_i += aux * eij[1] * vij[1];
	      gradvyz_i += aux * eij[1] * vij[2];
	      gradvzx_i += aux * eij[2] * vij[0];
	      gradvzy_i += aux * eij[2] * vij[1];
	      gradvzz_i += aux * eij[2] * vij[2];	      
	    }
	  
	  }
	}
      }
    }
  
  //-- The force is summed up to the kernel variable --
  fx[i] += fxi;
  fy[i] += fyi;
  if (dim == 3)
    fz[i] += fzi;

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
