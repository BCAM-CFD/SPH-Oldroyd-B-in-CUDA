/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include <iostream>
#include "cuda_runtime.h"
#include "class_system.h"

// function to print colloids info
void class_system::print_colloids(real* k_coll_x,
				  real* k_coll_y,
				  real* k_coll_z,
				  real* k_coll_vx,
				  real* k_coll_vy,
				  real* k_coll_vz,
				  real* k_coll_omegax,
				  real* k_coll_omegay,
				  real* k_coll_omegaz,				  
				  real* k_fx_colloid,
				  real* k_fy_colloid,
				  real* k_fz_colloid,
				  real* k_tx_colloid,
				  real* k_ty_colloid,
				  real* k_tz_colloid,
				  int step) {

  char filename[50];
  sprintf(filename, "colloids-%d.dat", step);
  FILE* file = fopen(filename, "w");
  if (!file) {
    printf("System print colloids error: Error opening the file %s\n", filename);
    return;
  }

  cudaMemcpy(this->coll_x, k_coll_x,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->coll_y, k_coll_y,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)  
    cudaMemcpy(this->coll_z, k_coll_z,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  
  cudaMemcpy(this->coll_vx, k_coll_vx,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->coll_vy, k_coll_vy,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)  
    cudaMemcpy(this->coll_vz, k_coll_vz,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);

  cudaMemcpy(this->coll_omegaz, k_coll_omegaz,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);  
  if (dim == 3)  {
    cudaMemcpy(this->coll_omegax, k_coll_omegax,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(this->coll_omegay, k_coll_omegay,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  }
  
  cudaMemcpy(this->fx_colloid, k_fx_colloid,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->fy_colloid, k_fy_colloid,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)
    cudaMemcpy(this->fz_colloid, k_fz_colloid,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);

  cudaMemcpy(this->tz_colloid, k_tz_colloid,
	     N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 3)    {
    cudaMemcpy(this->tx_colloid, k_tx_colloid,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(this->ty_colloid, k_ty_colloid,
	       N_colloids * sizeof(real), cudaMemcpyDeviceToHost);
  }  

  if (dim == 2)
    for (int i = 0; i < N_colloids; ++i)
      fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT "\n",
	      i,
	      this->coll_x[i],
	      this->coll_y[i],
	      this->coll_vx[i],
	      this->coll_vy[i],
	      this->coll_omegaz[i],	      
	      this->fx_colloid[i],
	      this->fy_colloid[i],
	      this->tz_colloid[i]);
  else //--- dim == 3  ---
    for (int i = 0; i < N_colloids; ++i)
      fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT "\n",
	      i,
	      this->coll_x[i],
	      this->coll_y[i],
	      this->coll_z[i],	      
	      this->coll_vx[i],
	      this->coll_vy[i],
	      this->coll_vz[i],
	      this->coll_omegax[i],
	      this->coll_omegay[i],
	      this->coll_omegaz[i],	      	      
	      this->fx_colloid[i],
	      this->fy_colloid[i],
	      this->fz_colloid[i],
	      this->tx_colloid[i],
	      this->ty_colloid[i],
	      this->tz_colloid[i]);
  fclose(file);

  printf("Colloids file written. Time step %d\n", step);
}
