#include "class_system.h"
#include "config.h"


// Constructor of the class system
void class_system::initialize_pointers() {

  fx_wall     = nullptr;
  fy_wall     = nullptr;
  fz_wall     = nullptr;
  coll_x      = nullptr;
  coll_y      = nullptr;
  coll_z      = nullptr;
  coll_vx     = nullptr;
  coll_vy     = nullptr;
  coll_vz     = nullptr;
  coll_omegax = nullptr;
  coll_omegay = nullptr;
  coll_omegaz = nullptr;
  coll_theta  = nullptr;
  coll_quat0  = nullptr;
  coll_quatx  = nullptr;
  coll_quaty  = nullptr;
  coll_quatz  = nullptr;   
  fx_colloid  = nullptr;
  fy_colloid  = nullptr;
  fz_colloid  = nullptr;
  tx_colloid  = nullptr;
  ty_colloid  = nullptr;
  tz_colloid  = nullptr;  
  x           = nullptr;
  y           = nullptr;
  z           = nullptr;
  vx          = nullptr;
  vy          = nullptr;
  vz          = nullptr;
  fx          = nullptr;
  fy          = nullptr;
  fz          = nullptr;  
  dens        = nullptr;
  press       = nullptr;
  mass        = nullptr;
  type        = nullptr;
  cxx_tensor  = nullptr;
  cxy_tensor  = nullptr;
  cyy_tensor  = nullptr;
  cxz_tensor  = nullptr;
  cyz_tensor  = nullptr;
  czz_tensor  = nullptr;  

  walls_start    = nullptr;
  walls_list     = nullptr;
  colloids_start = nullptr;
  colloids_list  = nullptr;
  x_center       = nullptr;
  y_center       = nullptr;
  z_center       = nullptr;  
}
