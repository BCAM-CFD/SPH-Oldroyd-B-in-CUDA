/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include <iostream>
#include "class_system.h"

// function to print system info
void class_system::print_info() {
  std::cout << " \n";
  std::cout << "------------------ System info -----------------------------"  << " \n";
  std::cout << "Number of particles, N                 = " << N << " \n";
  std::cout << "Number of particles reserved, Nmax     = " << Nmax << " \n";  
  if (dim == 3) 
    std::cout << "Box size, L                            = " << L[0] << " " << L[1] << " " << L[2] << " \n";
  else
    std::cout << "Box size, L                            = " << L[0] << " " << L[1]  << " \n";
  std::cout << "Number of dimensions, dim              = " << dim << " \n";
  std::cout << "Cutoff radius, rcut                    = " << rcut << " \n";
  std::cout << "Normalization const of the kernel, cw  = " << cw << " \n"; 
  std::cout << "density, rho                           = " << rho << " \n";
  std::cout << "Speed of sound, c                      = " << c << " \n";
  std::cout << "Speed of sound square, csq             = " << csq << " \n";
  std::cout << "P0                                     = " << P0 << " \n";
  std::cout << "Base density, rho0                     = " << rho0 << " \n";
  std::cout << "Transport coefficient a, a             = " << a << " \n";
  std::cout << "Transport coefficient b, b             = " << b << " \n";
  std::cout << "np * kT                                = " << nkT << "\n";
  std::cout << "tau                                    = " << tau << "\n";    
  std::cout << "Cell size                              = " << cell_size[0] << " " <<
    cell_size[1] << " " << cell_size[2] << " \n";
  std::cout << "Number of cells                        = " << Ncells[0] << " " << Ncells[1] <<
    " " << Ncells[2] << "\n";
  std::cout << "Total number of cells                  = " << Ntotal_cells << "\n";
  std::cout << "Type of external force, ext_force_type = " <<  ext_force_type << " \n";  
  if (dim == 3)
    std::cout << "External force per unit mass, f        = " << ext_force[0] << " " << ext_force[1] << " " << ext_force[2] << " \n";
  else
    std::cout << "External force per unit mass, f        = " << ext_force[0] << " " << ext_force[1]  << " \n";
  std::cout << "2 * pi * Ly = ky                       = " << ky << " \n";
  if (wall == 0)
    std::cout << "Wall                                   = No\n";
  else {//-- Wall = 1
    std::cout << "Wall                                   = Yes\n";
    std::cout << "Number of walls                        = " << Nwalls << "\n";    
    std::cout << "y coordinate of the bottom wall        = " << y_bottom << "\n";
    std::cout << "y coordinate of the top wall           = " << y_top << "\n";
    std::cout << "Velocity of the bottom wall            = " << V_bottom << "\n";
    std::cout << "Velocity of the top wall               = " << V_top << "\n";
    std::cout << "Number of particles in walls list, Nlist_walls = " << Nlist_walls << "\n";
    for (int i = 0; i < Nwalls; ++i)
      std::cout << "Start walls list " << i << "                     = " << walls_start[i] << "\n";      
  }
  std::cout << "Number of colloids, N_colloids         = " << N_colloids << "\n";
  if (N_colloids > 0) {
    std::cout << "Colloids radius                        = " << coll_R << "\n";
    std::cout << "Colloids density                       = " << coll_rho << "\n";
    std::cout << "Colloids mass                          = " << coll_mass << "\n";
    std::cout << "Colloids moment of inertia             = " << coll_I << "\n";        
    if (dim == 3)
      for (int i = 0; i < N_colloids; ++i)
	std::cout << "Initial Position of colloid " << i << "          = " 
		  << coll_x[i] << " "  	<< coll_y[i] << " " << coll_z[i] << "\n";
    else
      for (int i = 0; i < N_colloids; ++i)      
        std::cout << "Initial Position of colloid " << i << "          = " 
		  << coll_x[i] << " " << coll_y[i] << "\n";
    std::cout << "Number of particles in colloids list, Nlist_colloids = " << Nlist_colloids << "\n";
  }
  if (N_colloids >= 2) {
    std::cout << "colloids repulsion cutoff (surface-surface) = " << coll_rep_cutoff << "\n";
    std::cout << "colloids repulsion cuton (surface-surface)  = " << coll_rep_cuton << "\n";
    std::cout << "colloids repulsion cutoff (center-center)   = " << rcutoff_coll_rep << "\n";
    std::cout << "colloids repulsion cuton (center-center)    = " << rcuton_coll_rep << "\n";
    std::cout << "F0_repulsion                                = " << F0_rep << "\n";
    std::cout << "tau_repulsion                               = " << tau_rep << "\n";
    std::cout << "colloids cutoff                             = " << rcutoff_coll << "\n";
    std::cout << "Colloids cell size                          = " << cell_colloids_size[0]
	      << " " << cell_colloids_size[1] << " " << cell_colloids_size[2] << " \n";
    std::cout << "Number of colloids cells                    = " << Ncells_colloids[0] << " "
	      << Ncells_colloids[1] << " " << Ncells_colloids[2] << "\n";
    std::cout << "Total number of colloids cells              = "
	      << Ntotal_cells_colloids << "\n";   
  }

  if (wall == 1 || N_colloids > 0)
    std::cout << "Wall and/or colloids width                   = " << wall_width << "\n";

      
  std::cout << "------------------------------------------------------------"  << " \n";
  std::cout << " \n";  
}
