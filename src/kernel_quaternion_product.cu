#include "kernel_functions.h"
#include "config.h"

// Calculation of the product of the quaternions
//   c = c0 + cx i + cy j + cz k
//   d = d0 + dx i + dy j + dz k
__device__ void kernel_quaternion_product(real* c,
					  real* d,
					  real* result) {

  real c0 = c[0];
  real cx = c[1];
  real cy = c[2];
  real cz = c[3];
  real d0 = d[0];
  real dx = d[1];
  real dy = d[2];
  real dz = d[3];  
  
  result[0] = c0 * d0 - cx * dx - cy * dy - cz * dz;
  result[1] = c0 * dx + cx * d0 + cy * dz - cz * dy;
  result[2] = c0 * dy - cx * dz + cy * d0 + cz * dx;
  result[3] = c0 * dz + cx * dy - cy * dx + cz * d0;
}
