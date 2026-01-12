#include <iostream>
#include "cuda_runtime.h"
#include "class_system.h"

// function to print system info
void class_system::print_output(dim3 numBlocks,
				dim3 threadsPerBlock,
				real* k_x,
				real* k_y,
				real* k_z,
				real* k_vx,
				real* k_vy,
				real* k_vz,
				real* k_fx,
				real* k_fy,
				real* k_fz,
				real* k_mass,				   
				real* k_dens,
				real* k_press,
				real* k_cxx,
				real* k_cxy,
				real* k_cyy,
				real* k_cxz,
				real* k_cyz,
				real* k_czz,
				real* k_fx_wall,
				real* k_fy_wall,
				real* k_fz_wall,
				real* k_coll_x,
				real* k_coll_y,
				real* k_coll_z,
				real* k_coll_vx,
				real* k_coll_vy,
				real* k_coll_vz,
				real* k_coll_omegax,
				real* k_coll_omegay,
				real* k_coll_omegaz,
				real* k_fx_colloid,
				real* k_fy_colloid,
				real* k_fz_colloid,
				real* k_tx_colloid,
				real* k_ty_colloid,
				real* k_tz_colloid,
				real* k_kin_energy,
				int   freq_micro,
				int   freq_walls,
				int   freq_colloids,
				int   freq_macro,				
				int   step) {

  //---- Particles are written ----
  if (step % freq_micro == 0)
    print_particles(k_x, k_y, k_z, k_vx, k_vy, k_vz, k_mass, k_dens, k_press,
		    k_cxx, k_cxy, k_cyy, k_cxz, k_cyz, k_czz, step);
  //---- Walls are written ----
  if (wall == 1)
    if (step % freq_walls== 0)
      print_walls(k_fx_wall, k_fy_wall, k_fz_wall, step);
  //---- Colloids are written ----
  if (N_colloids > 0)
    if (step % freq_colloids == 0)
      print_colloids(k_coll_x, k_coll_y, k_coll_z,
		     k_coll_vx, k_coll_vy, k_coll_vz,
		     k_coll_omegax, k_coll_omegay, k_coll_omegaz,			 
		     k_fx_colloid, k_fy_colloid, k_fz_colloid,
		     k_tx_colloid, k_ty_colloid, k_tz_colloid, step);    
  //---- Macro variables are written ----
  if (step % freq_macro == 0) 
    print_macro_vars(numBlocks, threadsPerBlock, k_mass, k_vx, k_vy, k_vz,
		     k_kin_energy, step);

}
