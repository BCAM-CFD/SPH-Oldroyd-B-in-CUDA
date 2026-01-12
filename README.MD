# SPH in CUDA
SPH for the Giesekus model with wall boundaries and suspensions of spherical particles, implemented in CUDA.

Developed by Adolfo Vázquez-Quesada. 

mail: a.vazquez-quesada@fisfun.uned.es
--------------------------------------------------

- Host variables and functions are declared in class_system.h.
  Each function is defined in its own file, specifically in system_* files.

- Device variables and functions are declared in kernel_functions.h.
  Each function is defined in its own file, specifically in kernel_* files.

- In config.h we can define if the real variables are float or double.

- main.cu is the main file of the code.

*******************************
******** input variables ******
*******************************
N       -> Number of particles in each direction ( in a rectangular grid)

L       -> Box length in each direction.

dim     -> Number of dimensions

dt      -> Time steps

Nsteps  -> Number of steps

overlap -> rcut / dx   (where dx = L/N)

rho     -> fluid density

c       -> Speed of sound

P0      -> parameter of the equation of state

eta     -> shear viscosity

zeta    -> bulk viscosity

np      -> number of polymers per particle

kT      -> Boltzmann constant multiplied by temperature

tau     -> Relaxation time of the Oldroyd-B model

alpha   -> alpha parameter of the Giesekus model

ext_force_type -> Type of external force:

	       0: No external force

		   1: Constant force

		   2: sine force

ext_force -> External force vector

freq_micro    -> Frequency to write micro file.

freq_macro    -> Frequency to write macro file.

freq_walls    -> Frequency to write walls file.

freq_colloids -> Frequency to write colloids file.

wall -> walls are normal to the y-direction.

		0: no walls

     	1: walls

Vwall_bottom -> x-velocity of the bottom wall.

Vwall_top    -> x-velocity of the top wall.

N_colloids -> Number of colloidal particles.

coll_R     -> Colloids radius.

coll_rho   -> Colloids density.

coll_move  -> 
			  0: colloids are moving.

	   		  1: colloids are not moving.

coll_x     -> Position of one colloid. Put as many coll_x as needed.

coll_repulsion_cuton -> Cut on of the repulsion force between colloids.

F0_repulsion         -> Magnitude of the repulsion force.

tau_repulsion        -> Tau parameter of the repulsion force (related to the cut off radius of the repulsion force).

cutoff_magnetic -> Cut off radius of the magnetic force.

F0_magnetic     -> F0 parameter of the magnetic force.

omega_magnetic  -> Angular velocity of the rotation of the magnetic field.

new_sim -> 
		   0: new simulation.

	       1: existing simulation (the orientation of the magnetic field is not currently read in this option)

************************************
Key of the output files
************************************
micro: 

		-- 2D --

       1  -> Particle id

       2  -> x

       3  -> y

       4  -> vel x

       5  -> vel y

       6  -> mass

       7  -> density

       8  -> pressure

       9  -> type of particle

       10 -> cxx  

       11 -> cxy  

       12 -> cyy  	

       -- 3D --

       1  -> Particle id

       2  -> x

       3  -> y

       4  -> z   

       5  -> vel x

       6  -> vel y

       7  -> vel z    

       8  -> mass

       9  -> density

       10 -> pressure

       11 -> type of particle

       12 -> cxx

       13 -> cxy

       14 -> cyy

       15 -> cxz

       16 -> cyz

       17 -> czz       

walls: 

		-- 2D --

       1 -> wall identity

       2 -> force x

       3 -> force y

       -- 3D --

       1 -> wall identity

       2 -> force x

       3 -> force y

       4 -> force z

colloids: 

	  -- 2D --

	  1 -> x

	  2 -> y

	  3 -> vel x

	  4 -> vel y

	  5 -> ang vel

	  6 -> force x

	  7 -> force y

	  8 -> torque

	  -- 3D --

	  1  -> x

	  2  -> y

	  3  -> z	 

	  4  -> vel x

	  5  -> vel y

	  6  -> vel z	

	  7  -> ang vel x

	  8  -> ang vel y

	  9  -> ang vel z	

	  10 -> force x

	  11 -> force y

	  12 -> force z	  

	  13 -> torque x

	  14 -> torque y

	  15 -> torque z

macro: 

	   1 -> Step

       2 -> Total kinetic energy

*****************************
To restart a simulation
*****************************
1.- Use the same input file of the initial simulation (or some variation of it).

2.- In the input file: new_sim 1

3.- Create a file 'input_particles.in' with the following columns:

       --- 2D --

	   1.- x

	   2.- y

	   3.- Vel x

	   4.- Vel y

	   5.- mass

	   6.- type of particle

	   7.- cxx

	   8.- cxy

	   9.- cyy	   

       --- 3D --

	   1.-  x

	   2.-  y

	   3.-  z	 

	   4.-  Vel x

	   5.-  Vel y

	   6.-  Vel z	

	   7.-  mass

	   8.-  type of particle

	   9.-  cxx

	   10.- cxy

	   11.- cyy

	   12.- cxz

	   13.- cyz

	   14.- czz	  

4.- Run the simulation
