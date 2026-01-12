#include "kernel_functions.h"
#include "config.h"

// Calculation of the gradient of the SPH kernel: quintic spline
__device__ real kernel_gradW(real r) {

  real gradW;
  real r_rcut = r * rcut_inv;
  real r_rcut_3 = 3.0 * r_rcut;
  if (r_rcut_3 < 1.0) {
    real gW1 = 3.0 - r_rcut_3;
    real gW1_2 = gW1 * gW1;
    real gW1_4 = gW1_2 * gW1_2;
    real gW2 = 2.0 - r_rcut_3;
    real gW2_2 = gW2 * gW2;
    real gW2_4 = gW2_2 * gW2_2;
    gW2_4 = 6.0 * gW2_4;
    real gW3 = 1.0 - r_rcut_3;
    real gW3_2 = gW3 * gW3;
    real gW3_4 = gW3_2 * gW3_2;
    gW3_4 = 15.0 * gW3_4;
    gradW = gW1_4 - gW2_4;
    gradW = gradW + gW3_4;
    gradW = c_gradw * gradW;
    }
  else if (r_rcut_3 < 2.0) {
    real gW1 = 3.0 - r_rcut_3;
    real gW1_2 = gW1 * gW1;
    real gW1_4 = gW1_2 * gW1_2;
    real gW2 = 2.0 - r_rcut_3;
    real gW2_2 = gW2 * gW2;
    real gW2_4 = gW2_2 * gW2_2;
    gW2_4 = 6.0 * gW2_4;
    gradW = gW1_4 - gW2_4;
    gradW = c_gradw * gradW;	
  }
  else {
    real gW1 = 3.0 - r_rcut_3;
    real gW1_2 = gW1 * gW1;
    real gW1_4 = gW1_2 * gW1_2;
    gradW = c_gradw * gW1_4;
  }
  
  return gradW;
}
