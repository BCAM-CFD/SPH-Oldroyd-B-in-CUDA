/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

// Calculation of the pressure with an equation of state,
// and some auxiliar variables
__global__ void kernel_pressures(real* __restrict__ press,
				 real* __restrict__ dens,
				 real* __restrict__ mass,
				 real* __restrict__ inv_dens,
				 real* __restrict__ cxx,
				 real* __restrict__ cxy,
				 real* __restrict__ cyy,
				 real* __restrict__ cxz,
				 real* __restrict__ cyz,
				 real* __restrict__ czz,				 
				 real* __restrict__ Pixx_over_dens_sq,
				 real* __restrict__ Pixy_over_dens_sq,
				 real* __restrict__ Piyy_over_dens_sq,
				 real* __restrict__ Pixz_over_dens_sq,
				 real* __restrict__ Piyz_over_dens_sq,
				 real* __restrict__ Pizz_over_dens_sq) {
  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;
  
  real press_i = csq * (mass[i] * dens[i] - rho0);
  press[i] = press_i;

  real inv_dens_i = inv_dens[i];
  real inv_dens_i_sq = inv_dens_i * inv_dens_i;

  // We consider c tensor symmetric, so Pi tensor is also symmetric
  //-- Pressure tensor Pi = p - nkT * c  
  Pixx_over_dens_sq[i] = (press_i - nkT * cxx[i]) * inv_dens_i_sq;
  Pixy_over_dens_sq[i] =          - nkT * cxy[i]  * inv_dens_i_sq;
  Piyy_over_dens_sq[i] = (press_i - nkT * cyy[i]) * inv_dens_i_sq;
  if (dim == 3) {
    Pixz_over_dens_sq[i] =         - nkT * cxz[i]   * inv_dens_i_sq;
    Piyz_over_dens_sq[i] =         - nkT * cyz[i]   * inv_dens_i_sq;
    Pizz_over_dens_sq[i] = (press_i - nkT * czz[i]) * inv_dens_i_sq;  
  }
}
