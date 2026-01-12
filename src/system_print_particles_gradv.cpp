#include <iostream>
#include "cuda_runtime.h"
#include "class_system.h"

// function to print system info
void class_system::print_particles_gradv(real* k_x,
					 real* k_y,
					 real* k_z,
					 real* k_gradvxx,
					 real* k_gradvxy,
					 real* k_gradvxz,
					 real* k_gradvyx,
					 real* k_gradvyy,
					 real* k_gradvyz,
					 real* k_gradvzx,
					 real* k_gradvzy,
					 real* k_gradvzz,
					 int step) {
  int Npart = this->N;

  char filename[50];
  sprintf(filename, "gradv-%d.dat", step);
  FILE* file = fopen(filename, "w");
  if (!file) {
    printf("System print particle gradv error: Error opening the file %s\n", filename);
    return;
  }

  real* gradvxx = new real[Npart];
  real* gradvxy = new real[Npart];
  real* gradvxz = new real[Npart];
  real* gradvyx = new real[Npart];
  real* gradvyy = new real[Npart];
  real* gradvyz = new real[Npart];
  real* gradvzx = new real[Npart];
  real* gradvzy = new real[Npart];
  real* gradvzz = new real[Npart];  
  cudaMemcpy(this->x, k_x,       Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->y, k_y,       Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->z, k_z,       Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvxx, k_gradvxx, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvxy, k_gradvxy, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvxz, k_gradvxz, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvyx, k_gradvyx, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvyy, k_gradvyy, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvyz, k_gradvyz, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvzx, k_gradvzx, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvzy, k_gradvzy, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(gradvzz, k_gradvzz, Npart * sizeof(real), cudaMemcpyDeviceToHost);    
  if (dim == 2)
    for (int i = 0; i < Npart; ++i)
      fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT "\n",
	      i,
	      this->x[i],
	      this->y[i],
	      gradvxx[i],
	      gradvxy[i],
	      gradvyx[i],
	      gradvyy[i]);
  else //--- dim == 3  ---
    for (int i = 0; i < Npart; ++i)
      fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT "\n",
	      i,
	      this->x[i],
	      this->y[i],
	      this->z[i],
	      gradvxx[i],
	      gradvxy[i],
	      gradvxz[i],
	      gradvyx[i],
	      gradvyy[i],
	      gradvyz[i],
	      gradvzx[i],
	      gradvzy[i],
	      gradvzz[i]);    
  fclose(file);

  printf("gradv file written. Time step %d\n", step);
}
