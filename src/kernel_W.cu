#include "kernel_functions.h"
#include "config.h"

// Calculation of the SPH kernel: quintic spline
__device__ real kernel_W(real r) {
  
  real W;
  real r_rcut = r * rcut_inv;
  real r_rcut_3 = r_rcut * 3.0;    
  if (r_rcut_3 < 1.0) {
    real W1 = 3.0 - r_rcut_3;
    real W1_2 = W1 * W1;
    real W1_4 = W1_2 * W1_2;
    real W1_5 = W1_4 * W1;
    real W2 = 2.0 - r_rcut_3;
    real W2_2 = W2 * W2;
    real W2_4 = W2_2 * W2_2;      
    real W2_5 = W2_4 * W2;
    W2_5 = 6.0 * W2_5;
    real W3 = 1.0 - r_rcut_3;
    real W3_2 = W3 * W3;
    real W3_4 = W3_2 * W3_2;      
    real W3_5 = W3_4 * W3;
    W3_5 = 15.0 * W3_5;
    W = W1_5 - W2_5;
    W = W + W3_5;
    W = cw * W;
  }
  else if (r_rcut_3 < 2.0) {
    real W1 = 3.0 - r_rcut_3;
    real W1_2 = W1 * W1;
    real W1_4 = W1_2 * W1_2;
    real W1_5 = W1_4 * W1;
    real W2 = 2.0 - r_rcut_3;
    real W2_2 = W2 * W2;
    real W2_4 = W2_2 * W2_2;      
    real W2_5 = W2_4 * W2;
    W2_5 = 6.0 * W2_5;
    W = W1_5 - W2_5;
    W = cw * W;	
  }	
  else  {
    real W1 = 3.0 - r_rcut_3;
    real W1_2 = W1 * W1;
    real W1_4 = W1_2 * W1_2;
    real W1_5 = W1_4 * W1;	
    W = cw * W1_5;
  }

  return W;
}
