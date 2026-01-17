/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#ifndef KERNEL_FUNCTIONS_H
#define KERNEL_FUNCTIONS_H

#include "config.h"

// Constants
//-- The kernel constants are declared at kernel_declare_constants.cu --
extern __constant__ int N;
extern __constant__ int dim;
extern __constant__ real L[3];
extern __constant__ real cw;
extern __constant__ real c_gradw;
extern __constant__ real rcut;
extern __constant__ real rcutsq;
extern __constant__ real rcut_inv;
extern __constant__ real csq;
extern __constant__ real rho0;
extern __constant__ real a;
extern __constant__ real b;
extern __constant__ real nkT;
extern __constant__ real tau;
extern __constant__ real dt;
extern __constant__ real cell_size[3];
extern __constant__ int  Ncells[3];
extern __constant__ int  Ntotal_cells;
extern __constant__ int  ext_force_type;
extern __constant__ real ext_force[3];
extern __constant__ real ky;
extern __constant__ int wall;
extern __constant__ real y_bottom;
extern __constant__ real y_top;
extern __constant__ real V_bottom;
extern __constant__ real V_top;
extern __constant__ int  Nwalls;
extern __constant__ int  Nlist_walls;
extern __constant__ int  N_colloids;
extern __constant__ int  Nlist_colloids;
extern __constant__ real coll_R;
extern __constant__ real coll_rho;
extern __constant__ real coll_mass;
extern __constant__ real coll_I;
extern __constant__ real beta_max;
extern __constant__ real wall_width;
extern __constant__ real coll_rep_cutoff;
extern __constant__ real coll_rep_cuton;
extern __constant__ real cell_colloids_size[3];
extern __constant__ int  Ncells_colloids[3];
extern __constant__ int  Ntotal_cells_colloids;
extern __constant__ real rcutoff_coll_rep;
extern __constant__ real rcuton_coll_rep;
extern __constant__ real rcutoff_coll_rep_sq;
extern __constant__ real rcuton_coll_rep_sq;
extern __constant__ real rcutoff_coll;
extern __constant__ real rcutoff_coll_sq;
extern __constant__ real F0_rep;
extern __constant__ real tau_rep;

// Functions executed by the GPU
__global__ void kernel_print_particles(real* __restrict__ x,
				       real* __restrict__ y,
				       real* __restrict__ z,
				       real* __restrict__ vx,
				       real* __restrict__ vy,
				       real* __restrict__ vz,
				       real* __restrict__ fx,
				       real* __restrict__ fy,
				       real* __restrict__ fz,
				       real* __restrict__ dens,
				       real* __restrict__ press,
				       int*  __restrict__ particle_cell,
				       int*  __restrict__ particle_index);
__global__ void kernel_print_constants();
__global__ void kernel_print_cell_list(int* __restrict__ particle_cell,
				       int* __restrict__ particle_index);
__global__ void kernel_print_cell_ranges(int* __restrict__ cell_start,
					 int* __restrict__ cell_end);
__global__ void kernel_cell_test(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 int*  __restrict__ particle_cell,
				 int*  __restrict__ particle_index,
				 int*  __restrict__ cell_start,
				 int*  __restrict__ cell_end);
__global__ void kernel_cell_list(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 int*  __restrict__ particle_cell,
				 int*  __restrict__ particle_index);
__global__ void kernel_order_cell_list(int* __restrict__ particle_cell,
				       int* __restrict__ particle_index);
__global__ void kernel_cell_ranges(int* __restrict__ particle_cell,
				   int* __restrict__ cell_start,
				   int* __restrict__ cell_end);
__global__ void kernel_density_to_zero(real* __restrict__ dens);
__global__ void kernel_densities(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 real* __restrict__ dens,
				 int*  __restrict__ particle_index,
				 int*  __restrict__ cell_start,
				 int*  __restrict__ cell_end);
__global__ void kernel_pressures(real* __restrict__ press,
				 real* __restrict__ dens,
				 real* __restrict__ mass);
__global__ void kernel_vel_gradients_to_zero(real* __restrict__ gradvxx,
					     real* __restrict__ gradvxy,
					     real* __restrict__ gradvxz,
					     real* __restrict__ gradvyx,
					     real* __restrict__ gradvyy,
					     real* __restrict__ gradvyz,
					     real* __restrict__ gradvzx,
					     real* __restrict__ gradvzy,
					     real* __restrict__ gradvzz);
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
				     real* __restrict__ coll_omegaz);
__global__ void kernel_forces_to_zero(real* __restrict__ fx,
				      real* __restrict__ fy,
				      real* __restrict__ fz);
__global__ void kernel_forces(real* __restrict__ x,
			      real* __restrict__ y,
			      real* __restrict__ z,
			      real* __restrict__ vx,
			      real* __restrict__ vy,
			      real* __restrict__ vz,
			      real* __restrict__ cxx,
			      real* __restrict__ cxy,
			      real* __restrict__ cyy,
			      real* __restrict__ cxz,
			      real* __restrict__ cyz,
			      real* __restrict__ czz,
			      real* __restrict__ fx,
			      real* __restrict__ fy,
			      real* __restrict__ fz,
			      real* __restrict__ press,
			      real* __restrict__ dens,
			      real* __restrict__ mass,			      
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
			      real* __restrict__ coll_omegaz);
__global__ void kernel_collect_force_walls(real* __restrict__ fx_wall,
					   real* __restrict__ fy_wall,
					   real* __restrict__ fz_wall,
					   real* __restrict__ fx,
					   real* __restrict__ fy,
					   real* __restrict__ fz,
					   int*  __restrict__ walls_list,
					   int*  __restrict__ walls_start);
__global__ void kernel_collect_force_colloids(real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid,
					      real* __restrict__ fx,
					      real* __restrict__ fy,
					      real* __restrict__ fz,
					      real* __restrict__ x_center,
					      real* __restrict__ y_center,
					      real* __restrict__ z_center,
					      int*  __restrict__ colloids_list,
					      int*  __restrict__ colloids_start);
__global__ void kernel_ext_force(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 real* __restrict__ fx,
				 real* __restrict__ fy,
				 real* __restrict__ fz,
				 real* __restrict__ mass,
				 int*  __restrict__ type);
__global__ void kernel_move_sys_VV_part1(real* __restrict__ x,
					 real* __restrict__ y,
					 real* __restrict__ z,
					 real* __restrict__ vx,
					 real* __restrict__ vy,
					 real* __restrict__ vz,
					 real* __restrict__ cxx,
					 real* __restrict__ cxy,
					 real* __restrict__ cyy,
					 real* __restrict__ cxz,
					 real* __restrict__ cyz,
					 real* __restrict__ czz,
					 real* __restrict__ fx,
					 real* __restrict__ fy,
					 real* __restrict__ fz,
					 real* __restrict__ dcdt_xx,
					 real* __restrict__ dcdt_xy,
					 real* __restrict__ dcdt_yy,
					 real* __restrict__ dcdt_xz,
					 real* __restrict__ dcdt_yz,
					 real* __restrict__ dcdt_zz,
					 real* __restrict__ mass,
					 int*  __restrict__ type);
__global__ void kernel_move_sys_VV_part2(real* __restrict__ vx,
					 real* __restrict__ vy,
					 real* __restrict__ vz,
					 real* __restrict__ cxx,
					 real* __restrict__ cxy,
					 real* __restrict__ cyy,
					 real* __restrict__ cxz,
					 real* __restrict__ cyz,
					 real* __restrict__ czz,
					 real* __restrict__ fx,
					 real* __restrict__ fy,
					 real* __restrict__ fz,
					 real* __restrict__ dcdt_xx,
					 real* __restrict__ dcdt_xy,
					 real* __restrict__ dcdt_yy,
					 real* __restrict__ dcdt_xz,
					 real* __restrict__ dcdt_yz,
					 real* __restrict__ dcdt_zz,
					 real* __restrict__ mass,
					 int*  __restrict__ type);
__global__ void kernel_move_sys_Euler(real* __restrict__ x,
				      real* __restrict__ y,
				      real* __restrict__ z,
				      real* __restrict__ vx,
				      real* __restrict__ vy,
				      real* __restrict__ vz,
				      real* __restrict__ fx,
				      real* __restrict__ fy,
				      real* __restrict__ fz,
				      real* __restrict__ mass);
__global__ void kernel_move_colloids_VV_part1(real* __restrict__ coll_x,
					      real* __restrict__ coll_y,
					      real* __restrict__ coll_z,
					      real* __restrict__ coll_vx,
					      real* __restrict__ coll_vy,
					      real* __restrict__ coll_vz,
					      real* __restrict__ coll_omegax,
					      real* __restrict__ coll_omegay,
					      real* __restrict__ coll_omegaz,
					      real* __restrict__ coll_theta,
					      real* __restrict__ coll_quat0,
					      real* __restrict__ coll_quatx,
					      real* __restrict__ coll_quaty,
					      real* __restrict__ coll_quatz,
					      real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid,
					      real* __restrict__ x,
					      real* __restrict__ y,
					      real* __restrict__ z,
					      real* __restrict__ x_center,
					      real* __restrict__ y_center,
					      real* __restrict__ z_center,
					      int*  __restrict__ colloids_list,
					      int*  __restrict__ colloids_start);
__global__ void kernel_move_colloids_VV_part2(real* __restrict__ coll_vx,
					      real* __restrict__ coll_vy,
					      real* __restrict__ coll_vz,
					      real* __restrict__ coll_omegax,
					      real* __restrict__ coll_omegay,
					      real* __restrict__ coll_omegaz,
					      real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid);
__global__ void kernel_calculate_macro_vars(real* __restrict__ mass,
					    real* __restrict__ vx,
					    real* __restrict__ vy,
					    real* __restrict__ vz,
					    real* __restrict__ kin_energy);
__global__ void kernel_macro_vars_to_zero(real* __restrict__ kin_energy);
__global__ void kernel_cell_colloids_list(real* __restrict__ coll_x,
					  real* __restrict__ coll_y,
					  real* __restrict__ coll_z,
					  int*  __restrict__ coll_cell,
					  int*  __restrict__ coll_index);
__global__ void kernel_initialize_cell_colloids(int* __restrict__ coll_cell_start,
						int* __restrict__ coll_cell_end);
__global__ void kernel_cell_colloids_ranges(int* __restrict__ coll_cell,
					    int* __restrict__ coll_cell_start,
					    int* __restrict__ coll_cell_end);
__global__ void kernel_colloid_forces_to_zero(real* __restrict__ fx_colloid,
					      real* __restrict__ fy_colloid,
					      real* __restrict__ fz_colloid,
					      real* __restrict__ tx_colloid,
					      real* __restrict__ ty_colloid,
					      real* __restrict__ tz_colloid);
__global__ void kernel_repulsion_force_colloids(real* __restrict__ coll_x,
						real* __restrict__ coll_y,
						real* __restrict__ coll_z,
						real* __restrict__ fx_colloid,
						real* __restrict__ fy_colloid,
						real* __restrict__ fz_colloid,
						int*  __restrict__ coll_index,
						int*  __restrict__ coll_cell_start,
						int*  __restrict__ coll_cell_end);
__global__ void kernel_dcdt(real* __restrict__ dcdt_xx,
			    real* __restrict__ dcdt_xy,
			    real* __restrict__ dcdt_yy,
			    real* __restrict__ dcdt_xz,
			    real* __restrict__ dcdt_yz,
			    real* __restrict__ dcdt_zz,			    
			    real* __restrict__ cxx,
			    real* __restrict__ cxy,
			    real* __restrict__ cyy,
			    real* __restrict__ cxz,
			    real* __restrict__ cyz,
			    real* __restrict__ czz,
			    real* __restrict__ gradvxx,
			    real* __restrict__ gradvxy,
			    real* __restrict__ gradvxz,
			    real* __restrict__ gradvyx,
			    real* __restrict__ gradvyy,
			    real* __restrict__ gradvyz,
			    real* __restrict__ gradvzx,
			    real* __restrict__ gradvzy,
			    real* __restrict__ gradvzz,
			    int*  __restrict__ type);
__device__ real kernel_W(real r);
__device__ real kernel_gradW(real r);
__device__ void kernel_quaternion_product(real* c,
					  real* d,
					  real* result);

#endif // KERNEL_FUNCTIONS_H
