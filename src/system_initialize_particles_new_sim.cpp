/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "class_system.h"
#include "config.h"
#include <math.h>

#include <stdio.h>

// Initialization of the computational particles when they are not read from a file
int class_system::initialize_particles_new_sim(int  Nxyz[3],
					       real dx[3]) {
  // Initial positions
  int save = 0;
  int counter = 0;
  if (dim == 3) {
    for (int i = 0; i < Nxyz[0]; ++i) 
      for (int j = 0; j < Nxyz[1]; ++j) 
	for (int k = 0; k < Nxyz[2]; ++k) {
	  x[counter] = dx[0]/2.0 + i * dx[0];
	  y[counter] = dx[1]/2.0 + j * dx[1];
	  if (wall)
	    y[counter] = y[counter] - wall_width;	  
	  z[counter] = dx[2]/2.0 + k * dx[2];
	  //-- We check if the particle is inside any colloid --
	  for (int l = 0; l < N_colloids; ++l) 
	    if (save == 0) {
	      real xij = x[counter] - coll_x[l];
	      real yij = y[counter] - coll_y[l];
	      real zij = z[counter] - coll_z[l];
	      //--- Periodic boundary conditions ---
	      if (xij > 0.5 * L[0])
		xij -= L[0];
	      if (xij < -0.5 * L[0])
		xij += L[0];
	      if (wall == 0) {   //If there are not walls	      	    
		if (yij > 0.5 * L[1])
		  yij -= L[1];
		if (yij < -0.5 * L[1])
		  yij += L[1];
	      }	    	      
	      if (zij > 0.5 * L[2])
		zij -= L[2];
	      if (zij < -0.5 * L[2])
		zij += L[2];
	      real rijsq = xij*xij + yij*yij + zij*zij;
	      if (rijsq  <= coll_R*coll_R)
		save = 1;
	      else
		save = 0;
	    }
	  if (save == 0)
	    counter = counter + 1;
	  else
	    save = 0;
	}
  }
  else {  //----------- dim == 2 ----------------
    for (int i = 0; i < Nxyz[0]; ++i) 
      for (int j = 0; j < Nxyz[1]; ++j) {
	x[counter] = dx[0]/2.0 + i * dx[0];
	y[counter] = dx[1]/2.0 + j * dx[1];
	if (wall)
	  y[counter] = y[counter] - wall_width;
	//-- We check if the particle is inside any colloid --
	for (int l = 0; l < N_colloids; ++l) 
	  if (save == 0) {
	    real xij = x[counter] - coll_x[l];
	    real yij = y[counter] - coll_y[l];
	    //--- Periodic boundary conditions ---
	    if (xij > 0.5 * L[0])
	      xij -= L[0];
	    if (xij < -0.5 * L[0])
	      xij += L[0];
	    if (wall == 0) {   //If there are not walls	      	    
	      if (yij > 0.5 * L[1])
		yij -= L[1];
	      if (yij < -0.5 * L[1])
		yij += L[1];
	    }	    	      
	    real rijsq = xij*xij + yij*yij;
	    if (rijsq  <= coll_R*coll_R)
	      save = 1;
	    else
	      save = 0;
	  }
	if (save == 0)
	  counter = counter + 1;
	else
	  save = 0;		
      }
  }
  // The value of the number of particles is updated
  this->N = counter;

  //-- Types for fluid and walls particles are assigned--
  if (wall == 0)
    for (int i = 0; i < N; ++i)     
      type[i] = 0;
  else 
    for (int i = 0; i < N; ++i) 
      if (y[i] < y_bottom)
	type[i] = 1;
      else if ( y[i] > y_top)
	type[i] = 2;
      else
	type[i] = 0;  
      
  // Other initial quantities
  for (int i = 0; i < N; ++i) {
    // !--- At rest ---
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    //--- Sinusoidal ---
    // vx[i] = sin(ky * y[i]);
      
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
      
    if (dim == 2)
      mass[i] = rho * L[0] * (L[1] + 2.0 * wall_width) / (Nxyz[0] * Nxyz[1]);
    else // dim = 3
      mass[i] = rho * L[0] * (L[1] + 2.0 * wall_width) * L[2] /
	(Nxyz[0] * Nxyz[1] * Nxyz[2]);

    // Conformation tensor (valid for both 2D and 3D cases)
    cxx_tensor[i] = 1.0;
    cxy_tensor[i] = 0.0;
    cyy_tensor[i] = 1.0;
    cxz_tensor[i] = 0.0;
    cyz_tensor[i] = 0.0;
    czz_tensor[i] = 1.0;    
  }

  //--------- Initialization of colloids ----------------
  if (dim == 2) {
    int Nlayers = round(coll_R /dx[0]);
    for (int i = 0; i < N_colloids; ++i) {
      real r = 0.5 * dx[0];
      for (int j = 0; j < Nlayers; ++j) {
	if (r < coll_R && r > coll_R - wall_width) {
	  // Phi angle
	  real d_phi    = dx[0] / r;
	  int Npart_phi = round(2.0 * M_PI / d_phi);
	  d_phi         = 2.0 * M_PI / Npart_phi;
	  real phi      = 0.0;
	  for (int k = 0; k < Npart_phi; ++k) {
	    this->N = this->N + 1;
	    if (this->N > this->Nmax) {
	      printf("system initialize particles new sim error: N > Nmax. Please, increase the value of Nmax in system_constructor.\n");
	      return 1;
	    }
	      
	    x[this->N - 1]    = this->coll_x[i] + r * cos(phi);
	    y[this->N - 1]    = this->coll_y[i] + r * sin(phi);
	    type[this->N - 1] = 3 + i;
	    // Periodic boundary conditions
	    if (x[this->N - 1] < 0)
	      x[this->N - 1] = x[this->N - 1] + L[0];
	    if (wall == 0) //--- If no wall ---
	      if (y[this->N - 1] < 0)
		y[this->N - 1] = y[this->N - 1] + L[1];
	    if (x[this->N - 1] > L[0])
	      x[this->N - 1] = x[this->N - 1] - L[0];
	    if (wall == 0) //--- If no wall ---    
	      if (y[this->N - 1] > L[1])
		y[this->N - 1] = y[this->N - 1] - L[1];
	    phi = phi + d_phi;
	  }
	}
	r = r + dx[0];
      }
    }
  }
  else { // dim == 3   //Not perfect but good enough, I hope.
    int Nlayers = round(coll_R/dx[0]);
    for (int i = 0; i < N_colloids; ++i) {
      real r = 0.5 * dx[0];
      for (int j = 0; j < Nlayers; ++j) {
	if (r < coll_R && r > coll_R - wall_width) {
	  //-- First the north pole --
	  this->N        = this->N + 1;
	  if (this->N > this->Nmax) {
	    printf("system initialize particles new sim error: N > Nmax. Please, increase the value of Nmax in system_constructor.\n");
	    return 1;
	  }	  
	  x[this->N - 1] = this->coll_x[i];
	  y[this->N - 1] = this->coll_y[i];
	  z[this->N - 1] = this->coll_z[i] + r;
	  type[this->N - 1] = 3 + i;

	  //-- Theta angle --
	  real d_theta    = dx[0] / r;
	  int Npart_theta = round(M_PI / d_theta);
	  d_theta         = M_PI / Npart_theta;
	  real theta      = 0.0;
	  for (int k = 0; k < Npart_theta; ++k) {
	    //-- phi angle --
	    real d_phi    = dx[0]/(r * sin(theta));
	    int Npart_phi = round(2.0 * M_PI / d_phi);
	    d_phi         = 2.0 * M_PI / Npart_phi;
	    real phi      = 0.0;
	    for (int l = 0; l < Npart_phi; ++l) {
	      this->N = this->N + 1;
	      x[this->N - 1]    = this->coll_x[i] + r * sin(theta) * cos(phi);
	      // printf("sss %d %d %d %d\n", i, j, k, l);
	      // printf("ttt %d \n", this->N - 1);
	      y[this->N - 1]    = this->coll_y[i] + r * sin(theta) * sin(phi);
	      // printf("uuu\n");
	      // 421884
	      z[this->N - 1]    = this->coll_z[i] + r * cos(theta);
	      type[this->N - 1] = 3 + i;
	      // Periodic boundary conditions
	      if (x[this->N - 1] < 0)
		x[this->N - 1] = x[this->N - 1] + L[0];
	      if (wall == 0) //--- If no wall ---
		if (y[this->N - 1] < 0)
		  y[this->N - 1] = y[this->N - 1] + L[1];
	      if (z[this->N - 1] < 0)
		z[this->N - 1] = z[this->N - 1] + L[2];
	      if (x[this->N - 1] > L[0])
		x[this->N - 1] = x[this->N - 1] - L[0];
	      if (wall == 0) //--- If no wall ---    
		if (y[this->N - 1] > L[1])
		  y[this->N - 1] = y[this->N - 1] - L[1];
	      if (z[this->N - 1] > L[2])
		z[this->N - 1] = z[this->N - 1] - L[2];	      
	      phi = phi + d_phi;
	    }
	    theta = theta + d_theta;
	  }
	  //-- Finally, the south pole --
	  this->N = this->N + 1;
	  x[this->N - 1] = this->coll_x[i];
	  y[this->N - 1] = this->coll_y[i];
	  z[this->N - 1] = this->coll_z[i] - r;
	  type[this->N - 1] = 3 + i;
	}
	r = r + dx[0];
      }
    }
  }

  return 0;

}
