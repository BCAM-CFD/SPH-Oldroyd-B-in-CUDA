/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// Forces are set to zero
__global__ void kernel_vel_gradients_to_zero(real* __restrict__ gradvxx,
					     real* __restrict__ gradvxy,
					     real* __restrict__ gradvxz,
					     real* __restrict__ gradvyx,
					     real* __restrict__ gradvyy,
					     real* __restrict__ gradvyz,
					     real* __restrict__ gradvzx,
					     real* __restrict__ gradvzy,
					     real* __restrict__ gradvzz)  {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  //-------- Velocity gradient tensors are set to zero ---------
  gradvxx[i] = 0.0;
  gradvxy[i] = 0.0;
  gradvyx[i] = 0.0;
  gradvyy[i] = 0.0;  
  if (dim == 3) {
    gradvxz[i] = 0.0;
    gradvyz[i] = 0.0;
    gradvzx[i] = 0.0;
    gradvzy[i] = 0.0;
    gradvzz[i] = 0.0;
  }

}
