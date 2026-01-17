/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include <stdio.h>

__global__ void kernel_print_cell_list(int* __restrict__ particle_cell,
				       int* __restrict__ particle_index) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i >= N) return;
  
  printf("%d %d %d\n",
	 i,
	 particle_cell[i],
	 particle_index[i]);
}
