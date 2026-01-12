#include "class_system.h"
#include "config.h"
#include <math.h>
#include <stdio.h>

// Constructor of the class system
void class_system::constructor(int Nxyz[3],
			       int wall,
			       real overlap,
			       int dim,
			       int N_coll,
			       int new_sim) {

  this->dim = dim;

  //-- Walls_list and walls_start are constructed in system_initialize --  
  if (wall) {
    this->Nwalls = 2;
    fx_wall    = new real[this->Nwalls];
    fy_wall    = new real[this->Nwalls];
    fz_wall    = new real[this->Nwalls];
  }
  else
    this->Nwalls = 0;

  //-- Colloids_list and colloids_start are constructed in system_initialize --
  //   x_center, y_center and z_center are also constructed in system_initialize --
  if (N_coll > 0) {
    coll_x      = new real[N_coll];
    coll_y      = new real[N_coll];
    coll_z      = new real[N_coll];
    coll_vx     = new real[N_coll];
    coll_vy     = new real[N_coll];
    coll_vz     = new real[N_coll];
    coll_omegax = new real[N_coll];
    coll_omegay = new real[N_coll];
    coll_omegaz = new real[N_coll];
    coll_theta  = new real[N_coll];
    coll_quat0  = new real[N_coll];
    coll_quatx  = new real[N_coll];
    coll_quaty  = new real[N_coll];
    coll_quatz  = new real[N_coll];   
    fx_colloid  = new real[N_coll];
    fy_colloid  = new real[N_coll];
    fz_colloid  = new real[N_coll];
    tx_colloid  = new real[N_coll];
    ty_colloid  = new real[N_coll];
    tz_colloid  = new real[N_coll];
  }

  // If the particles are read from a file, these quantities are
  // constructed ar system_initialize_particles_old_sim
  if (new_sim == 0) {
    //-- Walls are perpendicular to the y direction --
    int Ny;
    int extra_N_wall = 2*ceil(overlap);
    if (wall == 1) // With walls
      Ny = Nxyz[1] + 2*ceil(overlap);
    else // Without walls
      Ny = Nxyz[1]; 
    
    if (dim == 2) 
      N = Nxyz[0] * Ny;
    else   //-- dim = 3 --
      N = Nxyz[0] * Ny * Nxyz[2];

    //--- Just in case, some more memory is reserved ---
    if (N_coll > 0)
      N = (int)(N * 1.01);
    Nmax = N;

    x     = new real[Nmax];
    y     = new real[Nmax];
    z     = new real[Nmax];
    
    vx    = new real[Nmax];
    vy    = new real[Nmax];
    vz    = new real[Nmax];
    
    fx    = new real[Nmax];
    fy    = new real[Nmax];
    fz    = new real[Nmax];  
    
    dens  = new real[Nmax];
    
    press = new real[Nmax];
    
    mass  = new real[Nmax];
    
    type  = new int[Nmax];

    cxx_tensor = new real[Nmax];
    cxy_tensor = new real[Nmax];
    cyy_tensor = new real[Nmax];
    cxz_tensor = new real[Nmax];
    cyz_tensor = new real[Nmax];
    czz_tensor = new real[Nmax];    
  }

}
