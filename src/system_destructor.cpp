#include "class_system.h"

// Destructor of the class system
void class_system::destructor() {
  if (x)
    delete[] x;
  if (y)
    delete[] y;
  if (z)
    delete[] z;

  if (vx)
    delete[] vx;
  if (vy)  
  delete[] vy;
  if (vz)  
  delete[] vz;

  if (fx)
    delete[] fx;
  if (fy)  
  delete[] fy;
  if (fz)  
  delete[] fz;

  if (cxx_tensor)
    delete[] cxx_tensor;
  if (cxy_tensor)
    delete[] cxy_tensor;
  if (cyy_tensor)
    delete[] cyy_tensor;
  if (cxz_tensor)
    delete[] cxz_tensor;
  if (cyz_tensor)
    delete[] cyz_tensor;
  if (czz_tensor)
    delete[] czz_tensor;
  
  if (mass)
    delete[] mass;

  if (dens)
    delete[] dens;

  if (press)
    delete[] press;

  if (type)
    delete[] type;

  if (fx_wall)
    delete[] fx_wall;
  if (fy_wall)
    delete[] fy_wall;
  if (fz_wall)
    delete[] fz_wall;
  if (walls_list)
    delete[] walls_list;

  if (coll_x)
    delete[] coll_x;
  if (coll_y)
    delete[] coll_y;
  if (coll_z)
    delete[] coll_z;
  if (coll_vx)
    delete[] coll_vx;
  if (coll_vy)
    delete[] coll_vy;
  if (coll_vz)
    delete[] coll_vz;
  if (coll_omegax)
    delete[] coll_omegax;
  if (coll_omegay)
    delete[] coll_omegay;
  if (coll_omegaz)
    delete[] coll_omegaz;
  if (coll_theta)
    delete[] coll_theta;
  if (coll_quat0)
    delete[] coll_quat0;
  if (coll_quatx)  
    delete[] coll_quatx;
  if (coll_quaty)  
    delete[] coll_quaty;
  if (coll_quatz)  
    delete[] coll_quatz;
  if (colloids_list)
    delete[] colloids_list;
  if (x_center)
    delete[] x_center;
  if (y_center)
    delete[] y_center;
  if (z_center)
    delete[] z_center;
  if (fx_colloid)
    delete[] fx_colloid;
  if (fy_colloid)
    delete[] fy_colloid;
  if (fz_colloid)
    delete[] fz_colloid;
  if (tx_colloid)
    delete[] tx_colloid;
  if (ty_colloid)
    delete[] ty_colloid;
  if (tz_colloid)
    delete[] tz_colloid;      
  
}
