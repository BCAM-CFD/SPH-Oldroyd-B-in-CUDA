#include <iostream>
#include "cuda_runtime.h"
#include "class_system.h"
#include "kernel_functions.h"

// function to print system info
void class_system::print_macro_vars(dim3 numBlocks,
				    dim3 threadsPerBlock,
				    real* k_mass,
				    real* k_vx,
				    real* k_vy,
				    real* k_vz,
				    real* kin_energy,
				    int step) {
  //--  Macro variables are set to zero --
  kernel_macro_vars_to_zero<<<1, 1>>>(kin_energy);
  cudaDeviceSynchronize();  // We require the kernel to end to continue

  //-- Macro variables are calculated --
  kernel_calculate_macro_vars<<<numBlocks, threadsPerBlock>>>(k_mass, k_vx, k_vy, k_vz, kin_energy);

  //-- Macro variables are copied into the host --
  cudaMemcpy(&this->kin_energy, kin_energy, sizeof(real), cudaMemcpyDeviceToHost);  

  //-- The file is opened --
  char filename[50];
  sprintf(filename, "macro.dat");
  FILE* file;
  if (step == 0)
    file = fopen(filename, "w");
  else
    file = fopen(filename, "a");
  if (!file)
    {
      printf("System print macro vars error: Error opening the file %s\n", filename);
      return;
    }
  
  fprintf(file, "%d %f\n", step, this->kin_energy);

  fclose(file);

  //  printf("Macro file written. Time step %d\n", step);
}
