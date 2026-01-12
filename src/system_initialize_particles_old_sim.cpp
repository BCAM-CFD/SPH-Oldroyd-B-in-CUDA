#include "class_system.h"
#include "config.h"

#include <iostream>
#include <fstream>
#include <vector>

// Initialization of the computational particles when they are read from a file
int class_system::initialize_particles_old_sim() {


  std::ifstream file;  // Declaration
  file.open("input_particles.in");  // Opening file
  
  if (!file.is_open()) {
    std::cerr << "Error al abrir el archivo." << std::endl;
    return 1;
  }

  //------ Variables to be read --------
  std::vector<real> x_r, y_r, z_r, vx_r, vy_r, vz_r, mass_r, type_r,
    cxx_r, cxy_r, cyy_r, cxz_r, cyz_r, czz_r;
  float xi, yi, zi, vxi, vyi, vzi, massi, typei,
    cxxi, cxyi, cyyi, cxzi, cyzi, czzi;

  //------ The file is read -------
  if (dim == 2)
    while (file >> xi >> yi >> vxi >> vyi >> massi >> typei >> cxxi >> cxyi >> cyyi) {
      x_r.push_back(xi);
      y_r.push_back(yi);
      vx_r.push_back(vxi);
      vy_r.push_back(vyi);
      mass_r.push_back(massi);
      type_r.push_back(typei);
      cxx_r.push_back(cxxi);
      cxy_r.push_back(cxyi);
      cyy_r.push_back(cyyi);      
    }
  else // dim == 3
    while (file >> xi >> yi >> zi >> vxi >> vyi >> vzi >> massi >> typei >> cxxi >> cxyi >> cyyi >> cxzi >> cyzi >> czzi) {
      x_r.push_back(xi);
      y_r.push_back(yi);
      z_r.push_back(zi);      
      vx_r.push_back(vxi);
      vy_r.push_back(vyi);
      vz_r.push_back(vzi);      
      mass_r.push_back(massi);
      type_r.push_back(typei);
      cxx_r.push_back(cxxi);
      cxy_r.push_back(cxyi);
      cyy_r.push_back(cyyi);      
      cxz_r.push_back(cxzi);
      cyz_r.push_back(cyzi);
      czz_r.push_back(czzi);
    }    

  file.close();

  //----- The vector is passed to the variables of the system object ------
  // N is defined
  this->N = x_r.size();
  // The arrays are constructed.
  x     = new real[N];
  y     = new real[N];
  z     = new real[N];
  vx    = new real[N];
  vy    = new real[N];
  vz    = new real[N];
  fx    = new real[N];
  fy    = new real[N];
  fz    = new real[N];  
  dens  = new real[N];
  press = new real[N];
  mass  = new real[N];
  type  = new int[N];
  cxx_tensor = new real[N];
  cxy_tensor = new real[N];
  cyy_tensor = new real[N];
  cxz_tensor = new real[N];
  cyz_tensor = new real[N];
  czz_tensor = new real[N];
  // The vectors are copied
  for (size_t i = 0; i < N; ++i) {
    x[i]    = x_r[i];
    y[i]    = y_r[i];
    vx[i]   = vx_r[i];
    vy[i]   = vy_r[i];
    mass[i] = mass_r[i];    
    type[i] = type_r[i];
    cxx_tensor[i] = cxx_r[i];
    cxy_tensor[i] = cxy_r[i];
    cyy_tensor[i] = cyy_r[i];    
    if (dim == 3) {
      z[i]  = z_r[i];
      vz[i] = vz_r[i];
      cxz_tensor[i] = cxz_r[i];
      cyz_tensor[i] = cyz_r[i];
      czz_tensor[i] = czz_r[i];
    }
  }

  return 0;

}
