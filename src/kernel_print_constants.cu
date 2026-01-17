/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"
#include <stdio.h>

//--- Constant variables storaged in the GPU are written ---
__global__ void kernel_print_constants() {

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i == 0)    {
    printf("----------------- Info storaged in the GPU -------------------------\n");
    printf("Number of particles, N            = %d\n", N);
    printf("Number of dimensions, dim         = %d\n", dim);
    printf("Box size, L                       = " REAL_FMT " " REAL_FMT " " REAL_FMT " \n", L[0], L[1], L[2]);
    printf("Kernel constant, cw               = " REAL_FMT " \n", cw);
    printf("Gradient kernel constant, c_gradw = " REAL_FMT " \n", c_gradw);
    printf("Cutoff radius, rcut               = " REAL_FMT " \n", rcut);
    printf("Square cutoff radius, rcutsq      = " REAL_FMT " \n", rcutsq);
    printf("Square speed of sound, csq        = " REAL_FMT " \n", csq);
    printf("Reference density, rho0           = " REAL_FMT " \n", rho0);
    printf("Transport coefficient a, a        = " REAL_FMT " \n", a);
    printf("Transport coefficient b, b        = " REAL_FMT " \n", b);
    printf("np * k * T (Oldroyd-B model)      = " REAL_FMT " \n", nkT);
    printf("relaxation time tau (Oldroyd-B model) = " REAL_FMT " \n", tau);        
    printf("Time step, dt                     = " REAL_FMT " \n", dt);
    printf("Cell size                         = " REAL_FMT " " REAL_FMT " " REAL_FMT " \n", cell_size[0], cell_size[1], cell_size[2]);
    printf("Number of cells                   = %d %d %d \n", Ncells[0], Ncells[1], Ncells[2]);
    printf("Total number of cells             = %d \n", Ntotal_cells);
    printf("External force type               = %d \n", ext_force_type);
    printf("External_force                    = " REAL_FMT " " REAL_FMT " " REAL_FMT "\n", ext_force[0], ext_force[1], ext_force[2]);
    printf("2*pi/Ly = ky                      = " REAL_FMT "\n", ky);
    if (wall == 1){
      printf("Wall                              = yes\n");
      printf("y coord. of the bottom wall       = " REAL_FMT "\n", y_bottom);
      printf("y coord. of the top wall          = " REAL_FMT "\n", y_top);
      printf("Vel of the bottom wall            = " REAL_FMT "\n", V_bottom);
      printf("Vel of the top wall               = " REAL_FMT "\n", V_top);
      printf("Number of walls, Nwalls           = %d\n",Nwalls);
      printf("Number of particles of walls      = %d\n",Nlist_walls);
    }
    else
      printf("Wall                              = no\n");
    printf("Number of colloids                = %d \n", N_colloids);
    if (N_colloids > 0) {
      printf("Colloids radius                   = " REAL_FMT "\n", coll_R);
      printf("Colloids density                  = " REAL_FMT "\n", coll_rho);
      printf("Colloids mass                     = " REAL_FMT "\n", coll_mass);
      printf("Colloids moment of inertia        = " REAL_FMT "\n", coll_I);      
      printf("Number of particles of colloids   = %d\n",Nlist_colloids);
    }
    if (N_colloids >= 2) {
      printf("Colloids repulsion cutoff (surface-surface) = " REAL_FMT "\n", coll_rep_cutoff);
      printf("Colloids repulsion cuton (surface-surface)  = " REAL_FMT "\n", coll_rep_cuton);
      printf("Colloids repulsion cutoff (center-center)   = " REAL_FMT "\n", rcutoff_coll_rep);
      printf("Colloids repulsion cuton (center-center)    = " REAL_FMT "\n", rcuton_coll_rep);      
      printf("F0 repulsion coefficient          = " REAL_FMT "\n", F0_rep);
      printf("tau repulsion coefficient         = " REAL_FMT "\n", tau_rep);
      printf("Colloids cutoff                   = " REAL_FMT "\n", rcutoff_coll);
      printf("Colloids cell size                = " REAL_FMT " " REAL_FMT " " REAL_FMT " \n",
	     cell_colloids_size[0], cell_colloids_size[1], cell_colloids_size[2]);
      printf("Number of colloids cells          = %d %d %d \n", Ncells_colloids[0],
	     Ncells_colloids[1], Ncells_colloids[2]);
      printf("Total number of colloids cells    = %d \n", Ntotal_cells_colloids);      
    }
    if (wall == 1 || N_colloids > 0)	
      printf("Wall and/or colloids width        = " REAL_FMT "\n", wall_width);
    
    printf("--------------------------------------------------------------------\n");      
  }

}
