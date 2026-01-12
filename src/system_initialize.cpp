#include "class_system.h"
#include "config.h"
#include <math.h>

#include <stdio.h>

// Initialization of the system
int class_system::initialize(int  Nxyz[3],
			     real L[3],
			     int  dim,
			     real overlap,
			     real rho,
			     real c,
			     real P0,
			     real eta,
			     real zeta,
			     real np,
			     real kT,
			     real elast_time,			     
			     int  ext_force_type,
			     real ext_force[3],
			     int  wall,
			     real Vwall_top,
			     real Vwall_bottom,
			     int  N_coll,
			     real coll_radius,
			     real coll_density,
			     real coll_x_data[MAX_COLLOIDS][MAX_VALUES_PER_COLLOID],
			     real coll_repulsion_cuton,
			     real F0_repulsion,
			     real tau_repulsion,
			     int  new_sim) {
  int error_out;
  error_out = 0;

  this->dim  = dim;
  this->L[0] = L[0];
  this->L[1] = L[1];
  this->L[2] = L[2];
  this->wall = wall;
  this->V_top = Vwall_top;
  this->V_bottom = Vwall_bottom;
  this->nkT      = np * kT;
  this->tau      = elast_time;

  rcut = L[0] / Nxyz[0] * overlap;
  rcutsq = rcut * rcut;
  rcut_inv = 1.0 / rcut;
  
  // dx
  real dx[3];
  for (int i = 0; i < dim; ++i) 
    dx[i] = L[i]/Nxyz[i];    

  // 2 * pi / L
  this->ky = 2.0 * M_PI / L[1];

  //-- Walls are perpendicular to the y direction --
  // wall_width is calculated even if there are not walls, because is used also for colloids.
  real eps = 1.0e-10;  
  if (abs(overlap - floor(overlap)) < eps)   // If overlap = floor(overlap)
    this->wall_width = floor(overlap) * dx[1];
  else 
    this->wall_width = (floor(overlap) + 1) * dx[1];
  //The simulation box is enlarged if there are walls.
  if (wall == 1) {    //-- With wall --
    Nxyz[1]  = Nxyz[1] + 2*ceil(overlap);
    y_bottom = 0.0;
    y_top    = L[1];
  }

    if (dim == 2)
      cw = 63.0 / (478.0 * M_PI * rcut * rcut);
  else
    cw = 27.0 / (120.0 * M_PI * rcut * rcut * rcut);
  c_gradw = -15.0 *cw / rcut;

  // Related to the equation of state
  this->rho = rho;
  this->c   = c;
  csq       = c*c;
  this->P0  = P0;
  rho0      = P0 * rho;
  a         = 5.0 * eta/3.0  - zeta;
  b         = 5.0 * (zeta + eta/3.0);

  //About Morris boundary conditions
  this->beta_max = 1.5;

  // Neighbour cells
  Ncells[0]    = (int)(L[0]/rcut);
  cell_size[0] = L[0] / Ncells[0];  
  if (wall == 1) {
    Ncells[1]    = (int)((L[1] + 2 * wall_width)/rcut);
    cell_size[1] = (L[1] + 2 * wall_width) / Ncells[1];
  }
  else { // No walls
    Ncells[1]    = (int)(L[1]/rcut);
    cell_size[1] = L[1] / Ncells[1];      
  }
  if (dim == 3) {
    Ncells[2]    = (int)(L[2]/rcut);
    cell_size[2] = L[2] / Ncells[2];          
  }
  else { // dim == 2
    Ncells[2]    = 1;
    cell_size[2] = 1;
  }
  Ntotal_cells = Ncells[0] * Ncells[1] * Ncells[2];
  if (Ncells[0] < 3) {
    printf("System initialize error: neighbour cells are too big\n");
    error_out = 1;
    return 1;
  }
  if (Ncells[1] < 3) {
    printf("System initialize error: neighbour cells are too big\n");
    error_out = 1;
    return 1;
  }
  if (dim == 3)
    if (Ncells[2] < 3) {
      printf("System initialize error: neighbour cells are too big\n");
      error_out = 1;
      return 1;
    }

  // External ext_force
  this->ext_force_type = ext_force_type;
  this->ext_force[0]   = ext_force[0];
  this->ext_force[1]   = ext_force[1];
  if (dim == 3)
    this->ext_force[2] = ext_force[2];

  this->kin_energy = 0.0;
  
  //----- Initialization of colloids ------
  this->N_colloids      = N_coll;
  this->coll_R          = coll_radius;
  this->coll_rho        = coll_density;
  this->F0_rep          = F0_repulsion;
  this->tau_rep         = tau_repulsion;
  this->coll_rep_cuton  = coll_repulsion_cuton;
  this->coll_rep_cutoff = 5.0 * coll_R / tau_rep;
  // Colloid mass
  if (dim == 2) {
    this->coll_mass = M_PI * coll_R * coll_R * coll_rho;
    this->coll_I    = 0.5 * coll_mass * coll_R * coll_R;
  }
  else { // dim == 3
    this->coll_mass = 4.0/3.0 * M_PI * coll_R * coll_R * coll_R * coll_rho;
    this->coll_I    = 2.0/5.0 * coll_mass * coll_R * coll_R;
  }
  
  //--- Initial colloids positions and velocities ---
  if (this->N_colloids > 0)
    for (int i = 0; i < this->N_colloids; ++i) {
      //--- Positions ---
      this->coll_x[i] = coll_x_data[i][0];
      this->coll_y[i] = coll_x_data[i][1];
      if (dim == 3)
	this->coll_z[i] = coll_x_data[i][2];

      //--- Velocities 
      this->coll_vx[i] = 0.0;
      this->coll_vy[i] = 0.0;
      if (dim == 3)
	this->coll_vz[i] = 0.0;

      //--- Angular velocities 
      this->coll_omegax[i] = 0.0;
      this->coll_omegay[i] = 0.0;
      if (dim == 3)
	this->coll_omegaz[i] = 0.0;
      this->coll_theta[i]    = 0.0;
      this->coll_quat0[i]    = 1.0;
      this->coll_quatx[i]    = 0.0;
      this->coll_quaty[i]    = 0.0;
      this->coll_quatz[i]    = 0.0;      
      
      //We check if colloid i overlaps the walls
      if (wall) {
	if (this->coll_y[i] - coll_R < y_bottom) {
	  printf("System initialize error: colloid %d is overlapping "
		 "with the bottom wall\n", i);
	  return 1;
	}
	if (this->coll_y[i] + coll_R > y_top) {
	  printf("System initialize error: colloid %d is overlapping "
		 "with the top wall\n", i);
	  return 1;
	}
      }	
    }

  //--- Neighbour cells for colloids ----
  if (N_colloids > 0) {  
    rcuton_coll_rep       = 2.0* coll_R + coll_rep_cuton;
    rcuton_coll_rep_sq    = rcuton_coll_rep * rcuton_coll_rep;    
    rcutoff_coll_rep      = 2.0 * coll_R + coll_rep_cutoff;
    rcutoff_coll_rep_sq   = rcutoff_coll_rep * rcutoff_coll_rep;
    rcutoff_coll          = rcutoff_coll_rep;
    rcutoff_coll_sq       = rcutoff_coll * rcutoff_coll;
    Ncells_colloids[0]    = (int)(L[0]/rcutoff_coll);
    cell_colloids_size[0] = L[0] / Ncells_colloids[0];  
    Ncells_colloids[1]    = (int)(L[1]/rcutoff_coll);
    cell_colloids_size[1] = L[1] / Ncells_colloids[1];      
    if (dim == 3) {
      Ncells_colloids[2]    = (int)(L[2]/rcutoff_coll);
      cell_colloids_size[2] = L[2] / Ncells_colloids[2];          
    }
    else {
      Ncells_colloids[2]    = 1;
      cell_colloids_size[2] = 1;
    }
    Ntotal_cells_colloids = Ncells_colloids[0] * Ncells_colloids[1] * Ncells_colloids[2];
    if (Ncells_colloids[0] < 3) {
      printf("System initialize error: neighbour colloids cells are too big\n");
      printf("Ncells_colloids[0] = %d\n", Ncells_colloids[0]);
      error_out = 1;
      return 1;
    }
    if (Ncells_colloids[1] < 3) {
      printf("System initialize error: neighbour colloids cells are too big\n");
      printf("Ncells_colloids[1] = %d\n", Ncells_colloids[1]);    
      error_out = 1;
      return 1;
    }
    if (dim == 3)
      if (Ncells_colloids[2] < 3) {
	printf("System initialize error: neighbour colloids cells are too big\n");
	printf("Ncells_colloids[2] = %d\n", Ncells_colloids[2]);      
	error_out = 1;
      return 1;
      }
  }

  //--- Computational particles are initialized ---
  if (new_sim == 0) {
    error_out = this->initialize_particles_new_sim(Nxyz, dx);
    if (error_out != 0)
      return 1;
  }
  else  { // If the computational particles are read from a file
    error_out = this->initialize_particles_old_sim();
    if (error_out != 0)
      return 1;
  }
  
  // List of particles of walls is built
  // Construction and initialization of walls_list
  if (wall) {   //Nwalls > 0
    walls_start = new int[Nwalls];  
    Nlist_walls = 0;
    for (int i = 0; i < this->N; ++i) 
      if (type[i] == 1 || type[i] == 2)
	Nlist_walls = Nlist_walls + 1;
    // The array is allocated
    walls_list = new int[Nlist_walls];
    // walls_list is initialized
    int counter = 0;
    // Walls types are i+1
    for (int i = 0; i < Nwalls; ++i) {
      int counter_i = 0;
      for (int j = 0; j < this->N; ++j)	
	if (type[j] == i+1) {
	  walls_list[counter] = j;
	  if (counter_i == 0)
	    walls_start[i] = counter;
	  counter_i = counter_i + 1;
	  counter   = counter + 1;
	}
    }
  }

  // List of particles of colloids is built
  // Construction and initialization of colloids_list
  if (N_colloids > 0) {
    colloids_start = new int[N_colloids];    
    Nlist_colloids = 0;
    for (int i = 0; i < this->N; ++i) 
      if (type[i] > 2)
	Nlist_colloids = Nlist_colloids + 1;
    // The arrays are allocated
    colloids_list = new int[Nlist_colloids];
    // The array colloids_list is initialized
    // Colloids types are type + 3
    int counter = 0;
    for (int i = 0; i < N_colloids; ++i) {
      int counter_i = 0;
      for (int j = 0; j < this->N; ++j)
	if (type[j] == i + 3) {
	  colloids_list[counter] = j;
	  if (counter_i == 0) 
	    colloids_start[i] = counter;
	  counter_i = counter_i + 1;
	  counter = counter + 1;
	}
    }
    // Construction and initialization of x_center, y_center and z_center
    x_center = new real[Nlist_colloids];
    y_center = new real[Nlist_colloids];
    z_center = new real[Nlist_colloids];
    for (int i = 0; i < Nlist_colloids; ++i) {
      int part    = colloids_list[i];
      int colloid = type[part] - 3;
      x_center[i] = x[part] - coll_x[colloid];
      y_center[i] = y[part] - coll_y[colloid];
      if (dim == 3)
	z_center[i] = z[part] - coll_z[colloid];
      //--- Periodic boundary conditions ---
      if (x_center[i] > 0.5 * L[0])
	x_center[i] -= L[0];
      if (x_center[i] < -0.5 * L[0])
	x_center[i] += L[0];
      if (wall == 0) {   //If there are not walls	      	    
	if (y_center[i] > 0.5 * L[1])
	  y_center[i] -= L[1];
	if (y_center[i] < -0.5 * L[1])
	  y_center[i] += L[1];
      }
      if (dim == 3) {
	if (z_center[i] > 0.5 * L[2])
	  z_center[i] -= L[2];
	if (z_center[i] < -0.5 * L[2])
	  z_center[i] += L[2];
      }
    }
  }

  return 0;
  
}
