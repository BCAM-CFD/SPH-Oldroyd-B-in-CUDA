#include "kernel_functions.h"
#include "config.h"
#include <stdio.h>

//--- The table W is printed in the shell ---
__global__ void kernel_print_W() {

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i == 0){
    printf("----------------- W table storaged in the GPU -------------------------\n");
    for (int j = 0; j < Ntable; ++j)
      printf("%d %f\n", j, W[j]);
    printf("--------------------------------------------------------------------\n");      
  }
}
