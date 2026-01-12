#include "config.h"
#include <cuda_runtime.h>

struct class_system {
  int N;            // Number of particles
  int Nmax;         // reserved number of particles

  real* x;          // position
  real* y;
  real* z;

  real* vx;         // Velocity
  real* vy;
  real* vz;

  real* fx;         // force per mass unit
  real* fy;
  real* fz;

  real* dens;       // Numerical density
  
  real* press;      // pressure

  real* mass;       // mass

  int* type;        // type of particle

  real* cxx_tensor;   // Conformation tensor of the Oldroyd-B model
  real* cxy_tensor;
  real* cyy_tensor;
  real* cxz_tensor;
  real* cyz_tensor;
  real* czz_tensor;
  real nkT;  //np * k * T  (Oldroyd-B model)
  real tau;  // relaxation time of the Oldroyd-B model    

  real L[3];        // Box size

  int dim;          // Number of dimensions

  int   wall;       // wall = 0: No walls.  Wall = 1: walls perpendicular to y-axis.
  real  y_bottom;   // Coordinate y of the bottom wall
  real  y_top;      // Coordinate y of the top wall
  real  V_bottom;   // Velocity of the bottom wall
  real  V_top;      // Velocity of the top wall
  int   Nwalls;     // Number of walls
  real* fx_wall;    // force on walls
  real* fy_wall;
  real* fz_wall;
  int   Nlist_walls;  // Number of particles in the list of the walls
  int*  walls_list;   // List of computational particles of the walls
  int*  walls_start;  // Start of each wall in the walls_list array
  real  wall_width;   // Wall width

  int   N_colloids;     // Number of colloids
  real  coll_R;         // radius of colloids
  real  coll_rho;       // Density of colloids
  int   Nlist_colloids; // Number of particles in the list of colloids
  int*  colloids_list;  // List of computational particles of the colloids
  int*  colloids_start; // Start of each colloid in the colloids_list array  
  real* coll_x;         // positions of the colloids
  real* coll_y;
  real* coll_z;
  real* coll_vx;        // Velocities of colloids
  real* coll_vy;   
  real* coll_vz;
  real* coll_omegax;    // Angular velocities of colloids
  real* coll_omegay;   
  real* coll_omegaz;
  real* coll_theta;   // polar angle of the colloid
  real* coll_quat0;   // Component of quaternions of colloids
  real* coll_quatx;   // Component of quaternions of colloids
  real* coll_quaty;   // Component of quaternions of colloids
  real* coll_quatz;   // Component of quaternions of colloids  
  real  coll_mass;    // Colloids mass
  real  coll_I;       // Colloids moment of inertia  
  real* x_center;     // For colloids. x_center[i] is the x-distance from solid_particles[i]
  real* y_center;     // to the center of the corresponding colloid. Same with
  real* z_center;     // y_center and z_center
  // For 3D rotations of the colloids, we use quaternions q = (quat0, quatx, quaty, quatz) 
  real* quat0;
  real* quatx;
  real* quaty;
  real* quatz;            
  real* fx_colloid;   // force on colloids
  real* fy_colloid; 
  real* fz_colloid;
  real* tx_colloid;  // torque on colloids
  real* ty_colloid; 
  real* tz_colloid;
  real  coll_rep_cutoff;   // Colloid surface-colloid surface repulsion cutoff 
  real  coll_rep_cuton;    // Colloid surface-colloid surface repulsion cuton
  real  rcutoff_coll_rep;  // Colloid-colloid repulsion cutoff
  real  rcuton_coll_rep;   // Colloid-colloid repulsion cutoff
  real  rcutoff_coll_rep_sq;  // (rcutoff_coll_rep)^2    
  real  rcuton_coll_rep_sq;  // (rcuton_coll_rep)^2  
  real  rcutoff_coll;      // Colloid-colloid cutoff  (maximum of all cutoffs)
  real  rcutoff_coll_sq;   // (rcutoff_coll)^2
  real  F0_rep;          // F0 * tau_rep is the magnitude of the force
  real  tau_rep;

  real rcut;        // cutoff radius
  real rcutsq;      // Square cutoff radius
  real rcut_inv;    // Inverse cutoff radius  
  real cw;          // Normalization constant of W
  real c_gradw;     // Constant of gradW
  real c;           // speed of sound
  real csq;         // speed of sound squared  

  real rho;         // density  
  real P0;          // related to base density  
  real rho0;        // base density

  real a, b;        // transport coefficients (related to viscosity)

  real cell_size[3];    // Size of the cells
  int  Ncells[3];       // Number of cells
  int  Ntotal_cells;    // Total number of cells:
  real cell_colloids_size[3];    // Size of the cells for colloids
  int  Ncells_colloids[3];       // Number of cells for colloids
  int  Ntotal_cells_colloids;    // Total number of cells of colloids
  // Type of external force.
  // 0 -> No external force
  // 1 -> Constant external force
  // 2 -> Sinusoidal external force
  int  ext_force_type;  
  real ext_force[3];        // External force per SPH particle
  real ky;                  // 2 * pi / Ly
  real kin_energy;          // Kinetic energy

  real beta_max;   // Related to Morris boundary conditions

  /********* Subroutines **********/
  int read_input(int   N[3],
		 real  L[3],
		 int&  dim,
		 real& delta_t,
		 int&  Nsteps,
		 real& overlap,
		 real& rho,
		 real& c,
		 real& P0,
		 real& eta,
		 real& zeta,
		 real& np,
		 real& kT,
		 real& elast_time,		 
		 int&  ext_force_type,
		 real  ext_force[3],
		 int&  freq_micro,
		 int&  freq_macro,
		 int&  freq_walls,
		 int&  freq_colloid,
		 int&  wall,
		 real& Vwall_top,
		 real& Vwall_bottom,
		 int&  N_coll,
		 real& coll_radius,
		 real& col_density,
		 real  coll_x_data[MAX_COLLOIDS][MAX_VALUES_PER_COLLOID],
		 int&  coll_move,
		 real& coll_repulsion_cuton,
		 real& F0_repulsion,
		 real& tau_repulsion,
		 int&  new_sim);
  void initialize_pointers();
  void constructor(int  Nxyz[3],
		   int  wall,
		   real overlap,
		   int  dim,
		   int  N_coll,
		   int  new_sim);
  int initialize(int  Nxyz[3],
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
		 int  new_sim);
  int initialize_particles_new_sim(int  Nxyz[3],
				   real dx[3]);
  int initialize_particles_old_sim();    
  void print_info();
  void print_particles(real* k_x,
		       real* k_y,
		       real* k_z,
		       real* k_vx,
		       real* k_vy,
		       real* k_vz,
		       real* k_mass,		       
		       real* k_dens,
		       real* k_press,
		       real* k_cxx,
		       real* k_cxy,
		       real* k_cyy,
		       real* k_cxz,
		       real* k_cyz,
		       real* k_czz,		       
		       int step);
  void print_walls(real* k_fx_wall,
		   real* k_fy_wall,
		   real* k_fz_wall,
		   int step);
  void print_colloids(real* k_coll_x,
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
		      int step);
  void print_particles_gradv(real* k_x,
			     real* k_y,
			     real* k_z,
			     real* k_gradvxx,
			     real* k_gradvxy,
			     real* k_gradvxz,
			     real* k_gradvyx,
			     real* k_gradvyy,
			     real* k_gradvyz,
			     real* k_gradvzx,
			     real* k_gradvzy,
			     real* k_gradvzz,
			     int   step);
  int copy_pointers_to_device(real** k_x,
			      real** k_y,
			      real** k_z,
			      real** k_vx,
			      real** k_vy,
			      real** k_vz,
			      real** k_mass,
			      real** k_fx,
			      real** k_fy,
			      real** k_fz,
			      real** k_dens,
			      real** k_press,
			      int**  k_particle_cell,
			      int**  k_particle_index,
			      int**  k_cell_start,
			      int**  k_cell_end,
			      real** k_kin_energy,
			      int**  k_type,
			      real** k_cxx,
			      real** k_cxy,
			      real** k_cyy,
			      real** k_cxz,
			      real** k_cyz,
			      real** k_czz,
			      real** k_dcdt_xx,
			      real** k_dcdt_xy,
			      real** k_dcdt_yy,
			      real** k_dcdt_xz,
			      real** k_dcdt_yz,
			      real** k_dcdt_zz,
			      real** k_gradvxx,
			      real** k_gradvxy,
			      real** k_gradvxz,
			      real** k_gradvyx,
			      real** k_gradvyy,
			      real** k_gradvyz,
			      real** k_gradvzx,
			      real** k_gradvzy,
			      real** k_gradvzz,			       
			      real** k_fx_wall,
			      real** k_fy_wall,
			      real** k_fz_wall,
			      int**  k_walls_list,
			      int**  k_walls_start,
			      int**  k_colloids_list,
			      int**  k_colloids_start,		      
			      real** k_coll_x,
			      real** k_coll_y,
			      real** k_coll_z,
			      real** k_coll_vx,
			      real** k_coll_vy,
			      real** k_coll_vz,
			      real** k_coll_omegax,
			      real** k_coll_omegay,
			      real** k_coll_omegaz,
			      real** k_coll_theta,
			      real** k_coll_quat0,
			      real** k_coll_quatx,
			      real** k_coll_quaty,
			      real** k_coll_quatz,
			      real** k_x_center,
			      real** k_y_center,
			      real** k_z_center,
			      real** k_fx_colloid,
			      real** k_fy_colloid,
			      real** k_fz_colloid,
			      real** k_tx_colloid,
			      real** k_ty_colloid,
			      real** k_tz_colloid,
			      int**  k_coll_cell,
			      int**  k_coll_index,
			      int**  k_coll_cell_start,
			      int**  k_coll_cell_end);
  int neighbours_search(dim3  numBlocks,
			dim3  threadsPerBlock,
			real* k_x,
			real* k_y,
			real* k_z,
			int*  k_particle_cell,
			int*  k_particle_index,
			int*  k_cell_start,
			int*  k_cell_end,
			real* k_coll_x,
			real* k_coll_y,
			real* k_coll_z,
			int*  k_coll_cell,
			int*  k_coll_index,
			int*  k_coll_cell_start,
			int*  k_coll_cell_end);
  int calculate_density(dim3  numBlocks,
			dim3  threadsPerBlock,
			real* k_x,
			real* k_y,
			real* k_z,
			real* k_dens,
			int*  k_particle_index,
			int*  k_cell_start,
			int*  k_cell_end);
  int calculate_pressure(dim3 numBlocks,
			 dim3  threadsPerBlock,
			 real* k_press,
			 real* k_dens,
			 real* k_mass);
  int calculate_forces(dim3  numBlocks,
		       dim3  threadsPerBlock,
		       real* k_x,
		       real* k_y,
		       real* k_z,
		       real* k_vx,
		       real* k_vy,
		       real* k_vz,
		       real* k_cxx,
		       real* k_cxy,
		       real* k_cyy,
		       real* k_cxz,
		       real* k_cyz,
		       real* k_czz,		       
		       real* k_fx,
		       real* k_fy,
		       real* k_fz,
		       real* k_dcdt_xx,
		       real* k_dcdt_xy,
		       real* k_dcdt_yy,
		       real* k_dcdt_xz,
		       real* k_dcdt_yz,
		       real* k_dcdt_zz,				   
		       real* k_gradvxx,
		       real* k_gradvxy,
		       real* k_gradvxz,
		       real* k_gradvyx,
		       real* k_gradvyy,
		       real* k_gradvyz,
		       real* k_gradvzx,
		       real* k_gradvzy,
		       real* k_gradvzz,				   		       
		       real* k_press,
		       real* k_dens,
		       real* k_mass,				    
		       int*  k_particle_index,
		       int*  k_cell_start,
		       int*  k_cell_end,
		       int*  k_type,
		       real* k_fx_wall,
		       real* k_fy_wall,
		       real* k_fz_wall,
		       int*  k_walls_list,
		       int*  k_walls_start,		       
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
		       real* k_x_center,
		       real* k_y_center,
		       real* k_z_center, 
		       int*  k_colloids_list,
		       int*  k_colloids_start,
		       int*  k_coll_index,
		       int*  k_coll_cell_start,
		       int*  k_coll_cell_end);
  int calculate_external_forces(dim3  numBlocks,
				dim3  threadsPerBlock,
				real* k_x,
				real* k_y,
				real* k_z,
				real* k_fx,
				real* k_fy,
				real* k_fz,
				real* k_mass,
				int*  type);
  int calculate_force_walls(dim3  numBlocks,
			    dim3  threadsPerBlock,
			    real* k_fx_wall,
			    real* k_fy_wall,
			    real* k_fz_wall,
			    real* k_fx,
			    real* k_fy,
			    real* k_fz,
			    int*  k_walls_list,
			    int*  k_walls_start);
  int calculate_force_colloids(dim3  numBlocks,
			       dim3  threadsPerBlock,
			       real* k_coll_x,
			       real* k_coll_y,
			       real* k_coll_z,			       
			       real* k_fx_colloid,
			       real* k_fy_colloid,
			       real* k_fz_colloid,
			       real* k_tx_colloid,
			       real* k_ty_colloid,
			       real* k_tz_colloid,			       
			       real* k_fx,
			       real* k_fy,
			       real* k_fz,
			       real* k_x_center,
			       real* k_y_center,
			       real* k_z_center,
			       int*  k_colloids_list,
			       int*  k_colloids_start,
			       int*  k_coll_index,
			       int*  k_coll_cell_start,
			       int*  k_coll_cell_end);  
  int V_Verlet_step1(dim3 numBlocks,
		     dim3 threadsPerBlock,
		     real* k_x,
		     real* k_y,
		     real* k_z,
		     real* k_vx,
		     real* k_vy,
		     real* k_vz,
		     real* k_cxx,
		     real* k_cxy,
		     real* k_cyy,
		     real* k_cxz,
		     real* k_cyz,
		     real* k_czz,		     
		     real* k_fx,
		     real* k_fy,
		     real* k_fz,
		     real* k_dcdt_xx,
		     real* k_dcdt_xy,
		     real* k_dcdt_yy,
		     real* k_dcdt_xz,
		     real* k_dcdt_yz,
		     real* k_dcdt_zz,				 		     
		     real* k_mass,
		     int*  k_type,
		     real* k_coll_x,
		     real* k_coll_y,
		     real* k_coll_z,
		     real* k_coll_vx,
		     real* k_coll_vy,
		     real* k_coll_vz,
		     real* k_coll_omegax,
		     real* k_coll_omegay,
		     real* k_coll_omegaz,
		     real* k_coll_theta,
		     real* k_coll_quat0,
		     real* k_coll_quatx,
		     real* k_coll_quaty,
		     real* k_coll_quatz,
		     real* k_fx_colloid,
		     real* k_fy_colloid,
		     real* k_fz_colloid,
		     real* k_tx_colloid,
		     real* k_ty_colloid,
		     real* k_tz_colloid,		     
		     real* k_x_center,
		     real* k_y_center,
		     real* k_z_center,
		     int*  k_colloids_list,
		     int*  k_colloids_start,
		     int   coll_move);
  int V_Verlet_step2(dim3 numBlocks,
		     dim3 threadsPerBlock,
		     real* k_vx,
		     real* k_vy,
		     real* k_vz,
		     real* k_cxx,
		     real* k_cxy,
		     real* k_cyy,
		     real* k_cxz,
		     real* k_cyz,
		     real* k_czz,				 		     
		     real* k_fx,
		     real* k_fy,
		     real* k_fz,
		     real* k_dcdt_xx,
		     real* k_dcdt_xy,
		     real* k_dcdt_yy,
		     real* k_dcdt_xz,
		     real* k_dcdt_yz,
		     real* k_dcdt_zz,		     
		     real* k_mass,
		     int*  k_type,
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
		     int   coll_move);
  void print_macro_vars(dim3 numBlocks,
			dim3 threadsPerBlock,
			real* k_mass,
			real* k_vx,
			real* k_vy,
			real* k_vz,
			real* k_energy,
			int step);
  void print_output(dim3 numBlocks,
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
		    int   step);
  int calculate_vel_gradient(dim3 numBlocks,
			     dim3 threadsPerBlock,
			     real* k_x,
			     real* k_y,
			     real* k_z,
			     real* k_vx,
			     real* k_vy,
			     real* k_vz,
			     real* k_gradvxx,
			     real* k_gradvxy,
			     real* k_gradvxz,
			     real* k_gradvyx,
			     real* k_gradvyy,
			     real* k_gradvyz,
			     real* k_gradvzx,
			     real* k_gradvzy,
			     real* k_gradvzz,
			     real* k_dens,
			     int*  k_particle_index,
			     int*  k_cell_start,
			     int*  k_cell_end,
			     int*  k_type,
			     real* k_coll_x,
			     real* k_coll_y,
			     real* k_coll_z,
			     real* k_coll_vx,
			     real* k_coll_vy,
			     real* k_coll_vz,
			     real* k_coll_omegax,
			     real* k_coll_omegay,
			     real* k_coll_omegaz);  
  void destructor();
};
