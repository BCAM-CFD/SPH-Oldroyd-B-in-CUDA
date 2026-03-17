/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include <iostream>
#include "cuda_runtime.h"
#include "class_system.h"

// function to print walls info
void class_system::print_walls(real* k_fx_wall,
			       real* k_fy_wall,
			       real* k_fz_wall,
			       int step) {

  cudaMemcpy(this->fx_wall, k_fx_wall, Nwalls * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->fy_wall, k_fy_wall, Nwalls * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)
    cudaMemcpy(this->fz_wall, k_fz_wall, Nwalls * sizeof(real), cudaMemcpyDeviceToHost);  

  //-- The file is opened --
  char filename[50];
  sprintf(filename, "walls.dat");
  FILE* file;
  if (step == 0)
    file = fopen(filename, "w");
  else
    file = fopen(filename, "a");
  if (!file)
    {
      printf("System print walls error: Error opening the file %s\n", filename);
      return;
    }

  if (dim == 2)
    fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT  " " REAL_FMT "\n",
	    step,
	    fx_wall[0],
	    fy_wall[0],
	    fx_wall[1],
	    fy_wall[1]);
  else //--- dim == 3  ---
    fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT  " " REAL_FMT " " REAL_FMT " " REAL_FMT "\n",
	    step,
	    fx_wall[0],
	    fy_wall[0],
	    fz_wall[0],
	    fx_wall[1],
	    fy_wall[1],
	    fz_wall[1]);    
  fclose(file);

  //  printf("Walls file written. Time step %d\n", step);
}
