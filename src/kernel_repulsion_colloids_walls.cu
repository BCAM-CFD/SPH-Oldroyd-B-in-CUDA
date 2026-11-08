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
// Calculation of the repulsion force between wall and colloidal particles
// The repulsive force is given by
//     F^rep = F0 * tau * exp(-tau * s) / (1.0 - exp(-tau * s)
// where s is the adimensional distance between wall and colloid particle, tau
// is related to the interaction range. Specifically, for a distance
// R * tau^{-1} the force is small but not negiglible. For a distance
// coll_rep_cutoff = 5*R*tau^{-1} we consider it negiglible.
// *** In this subroutine we are supposing that there are walls at
//         y = 0 and y = Ly
__global__ void kernel_repulsion_colloids_walls(real* __restrict__ coll_y,
						real* __restrict__ fy_rep_colloid) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N_colloids) return;

  real coll_yi = __ldg(&coll_y[i]);
  real L1      = __ldg(&L[1]);
  real Famp    = F0_rep * tau_rep;
  real r;

  if (coll_yi - coll_R < coll_rep_cutoff) {   // Bottom wall
    if (coll_yi - coll_R < coll_rep_cuton) //the particle is too close to the bottom wall
      r = coll_rep_cuton;
    else
      r = coll_yi - coll_R;
    real tau_s_neg      = -tau_rep * r / coll_R;
    real exp_tau_s_neg  = exp(tau_s_neg);
    real Frep_mod       = Famp * exp_tau_s_neg / (1.0 - exp_tau_s_neg);
    fy_rep_colloid[i]   = fy_rep_colloid[i] + Frep_mod;    
  }
  else 
    if (coll_yi + coll_R > L1 - coll_rep_cutoff) {   // Top wall
      if (coll_yi + coll_R > L1 - coll_rep_cuton) //it is too close to the top wall
	r = coll_rep_cuton;
      else
	r = L1 - coll_yi - coll_R;
      real tau_s_neg      = -tau_rep * r / coll_R;
      real exp_tau_s_neg  = exp(tau_s_neg);
      real Frep_mod       = Famp * exp_tau_s_neg / (1.0 - exp_tau_s_neg);
      fy_rep_colloid[i]       = fy_rep_colloid[i] - Frep_mod;
    }
  
}
