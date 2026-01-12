#include "kernel_functions.h"
#include "config.h"

// Calculation of external forces
__global__ void kernel_ext_force(real* __restrict__ x,
				 real* __restrict__ y,
				 real* __restrict__ z,
				 real* __restrict__ fx,
				 real* __restrict__ fy,
				 real* __restrict__ fz,
				 real* __restrict__ mass,
				 int*  __restrict__ type) {
  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  if (type[i] == 0) { // Only fluids  
    if (ext_force_type == 1) { // constant force
      real fi[3] = {0.0, 0.0, 0.0};        
      real mass_i = mass[i];
      fi[0] = mass_i * ext_force[0];
      fi[1] = mass_i * ext_force[1];
      if (dim == 3)
	fi[2] = mass_i * ext_force[2];
      
      //-- The force is summed up to the kernel variable --
      fx[i] += fi[0];
      fy[i] += fi[1];
      if (dim == 3)
	fz[i] += fi[2];
    }
    else if (ext_force_type == 2) { //Sine force along x-direction
      real fi[3] = {0.0, 0.0, 0.0};            
      real mass_i = mass[i];    
      fi[0] = mass_i * ext_force[0] * sin(ky * y[i]);
      fi[1] = 0.0;
      if (dim == 3)
	fi[2] = 0.0;
      
      //-- The force is summed up to the kernel variable --
      fx[i] += fi[0];
      fy[i] += fi[1];
      if (dim == 3)
	fz[i] += fi[2];
    }
  }
}
