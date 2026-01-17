/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"
#include <math.h>

#include <stdio.h>
// Calculation of the time derivative of the conformation tensor
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
			    int*  __restrict__ type) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;
  if (type[i] != 0) return; //If the particle i is not of the fluid.

  real cxxi     = cxx[i];
  real cxyi     = cxy[i];
  real cyyi     = cyy[i];
  real cxzi;
  real cyzi;
  real czzi;  
  real gradvxxi = gradvxx[i];
  real gradvxyi = gradvxy[i];
  real gradvxzi;
  real gradvyxi = gradvyx[i];
  real gradvyyi = gradvyy[i];
  real gradvyzi;
  real gradvzxi;
  real gradvzyi;
  real gradvzzi;
  if (dim == 3) {
  cxzi     = cxz[i];
  cyzi     = cyz[i];
  czzi     = czz[i];
  gradvxzi = gradvxz[i];
  gradvyzi = gradvyz[i];
  gradvzxi = gradvzx[i];
  gradvzyi = gradvzy[i];
  gradvzzi = gradvzz[i];  
  }

  real cgradvxxi;
  real cgradvxyi;
  real cgradvxzi;
  real cgradvyxi;
  real cgradvyyi;
  real cgradvyzi;
  real cgradvzxi;
  real cgradvzyi;
  real cgradvzzi;    
  if (dim == 2) { // Remember that c is symmetric
    cgradvxxi = cxxi * gradvxxi + cxyi * gradvyxi;
    cgradvxyi = cxxi * gradvxyi + cxyi * gradvyyi;
    cgradvyxi = cxyi * gradvxxi + cyyi * gradvyxi;
    cgradvyyi = cxyi * gradvxyi + cyyi * gradvyyi;
  }
  else {   // dim == 3
    cgradvxxi = cxxi * gradvxxi + cxyi * gradvyxi + cxzi * gradvzxi;
    cgradvxyi = cxxi * gradvxyi + cxyi * gradvyyi + cxzi * gradvzyi;
    cgradvxzi = cxxi * gradvxzi + cxyi * gradvyzi + cxzi * gradvzzi;
    cgradvyxi = cxyi * gradvxxi + cyyi * gradvyxi + cyzi * gradvzxi;
    cgradvyyi = cxyi * gradvxyi + cyyi * gradvyyi + cyzi * gradvzyi;
    cgradvyzi = cxyi * gradvxzi + cyyi * gradvyzi + cyzi * gradvzzi;
    cgradvzxi = cxzi * gradvxxi + cyzi * gradvyxi + czzi * gradvzxi;
    cgradvzyi = cxzi * gradvxyi + cyzi * gradvyyi + czzi * gradvzyi;
    cgradvzzi = cxzi * gradvxzi + cyzi * gradvyzi + czzi * gradvzzi;    
  }

  real tau_inv = 1.0/tau;
  dcdt_xx[i] = 2.0 * cgradvxxi + tau_inv * (1.0 - cxxi);
  dcdt_xy[i] = cgradvxyi + cgradvyxi - tau_inv * cxyi;
  dcdt_yy[i] = 2.0 * cgradvyyi + tau_inv * (1.0 - cyyi);
  if (dim == 3) {
  dcdt_xz[i] = cgradvxzi + cgradvzxi - tau_inv * cxzi;
  dcdt_yz[i] = cgradvyzi + cgradvzyi - tau_inv * cyzi;
  dcdt_zz[i] = 2.0 * cgradvzzi + tau_inv * (1.0 - czzi);  
  }
  
}
