/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "class_system.h"
#include "config.h"
#include <math.h>

#define MAX_LINE 256

int class_system::read_input(int   N[3],
			     real  L[3],
			     int&  dim,
			     real& dt,
			     int&  Nsteps,
			     real& overlap,
			     real& rho,
			     real& c,
			     real& P0,
			     real& eta,
			     real& zeta,
			     int&  ext_force_type,
			     real  ext_force[3],
			     int&  freq_micro,
			     int&  freq_macro,
			     int&  freq_walls,
			     int&  freq_colloids,
			     int&  wall,
			     real& Vwall_top,
			     real& Vwall_bottom,
			     int&  N_coll,
			     real& coll_radius,
			     real& coll_density,
			     real  coll_x_data[MAX_COLLOIDS][MAX_VALUES_PER_COLLOID],
			     int&  coll_move,
			     real& coll_repulsion_cuton,
			     real& F0_repulsion,
			     real& tau_repulsion,
			     int&  new_sim) {

  // Array to store initial positions of the colloids	   
  int coll_x_count = 0;  
  
  FILE *file = fopen("input", "r");
  if (!file) {
    perror("Error al abrir el archivo");
    return 0;
  }

  char line[MAX_LINE];
  int dim_N, dim_L, dim_ext_force;
  int N_read               = 1;
  int L_read               = 1;
  int dim_read             = 1;
  int dt_read              = 1;
  int Nsteps_read          = 1;
  int overlap_read         = 1;
  int rho_read             = 1;
  int c_read               = 1;
  int P0_read              = 1;
  int eta_read             = 1;
  int zeta_read            = 1;
  int ext_force_type_read  = 1;    
  int ext_force_read       = 1;
  int freq_micro_read      = 1;
  int freq_walls_read      = 1;    
  int freq_macro_read      = 1;
  int freq_colloids_read   = 1;
  int wall_read            = 1;
  int Vwall_top_read       = 1;
  int Vwall_bottom_read    = 1;
  int N_coll_read          = 1;
  int coll_radius_read     = 1;
  int coll_density_read    = 1;
  int coll_move_read       = 1;
  int coll_repulsion_cuton_read  = 1;
  int F0_repulsion_read    = 1;
  int tau_repulsion_read   = 1;
  int new_sim_read         = 1;    
    
  while (fgets(line, sizeof(line), file)) {
    char var[50];
    double val1, val2, val3;

    // Elimina espacios al principio y al final
    char *start = line;
    while (*start == ' ' || *start == '\t') start++;
    char *end = start + strlen(start) - 1;
    while (end > start && (*end == ' ' || *end == '\t' || *end == '\n')) *end-- = '\0';

    // Extrae la variable y hasta 3 valores
    int count = sscanf(start, "%s %lf %lf %lf", var, &val1, &val2, &val3);
	
    if (strcmp(var, "N") == 0 && count == 4) {
      N[0]                 = (int)val1;
      N[1]                 = (int)val2;
      N[2]                 = (int)val3;
      dim_N                = 3;
      N_read               = 0;
    }
    else if (strcmp(var, "N") == 0 && count == 3) {
	N[0]               = (int)val1;
	N[1]               = (int)val2;
	dim_N              = 2;
	N_read             = 0;
      }
    else if (strcmp(var, "L") == 0 && count == 4) {
      L[0]                 = val1;
      L[1]                 = val2;
      L[2]                 = val3;
      dim_L                = 3;
      L_read               = 0;
    }
    else if (strcmp(var, "L") == 0 && count == 3) {
      L[0]                 = val1;
      L[1]                 = val2;
      dim_L                = 2;
      L_read               = 0;
      }
    else if (strcmp(var, "dim") == 0 && count >= 2) {
      dim                  = (int)val1;
      dim_read             = 0;
    }
    else if (strcmp(var, "dt") == 0 && count >= 2) {
      dt                   = val1;
      dt_read              = 0;
    } else if (strcmp(var, "Nsteps") == 0 && count >= 2) {
      Nsteps               = (int)val1;
      Nsteps_read          = 0;
    } else if (strcmp(var, "overlap") == 0 && count >= 2) {
      overlap              = val1;
      overlap_read         = 0;  
    } else if (strcmp(var, "rho") == 0 && count >= 2) {
      rho                  = val1;
      rho_read             = 0;
    } else if (strcmp(var, "c") == 0 && count >= 2) {
      c                    = val1;
      c_read               = 0;
    } else if (strcmp(var, "P0") == 0 && count >= 2) {
      P0                   = val1;
      P0_read              = 0;  
    } else if (strcmp(var, "eta") == 0 && count >= 2) {
      eta                  = val1;
      eta_read             = 0;
    } else if (strcmp(var, "zeta") == 0 && count >= 2) {
      zeta                 = val1;
      zeta_read            = 0;
    } else if (strcmp(var, "ext_force_type") == 0 && count >= 2) {
      ext_force_type       = (int)val1;
      ext_force_type_read  = 0;	  
    } else if (strcmp(var, "ext_force") == 0 && count == 4) {
      ext_force[0]         = val1;
      ext_force[1]         = val2;
      ext_force[2]         = val3;
      dim_ext_force        = 3;
      ext_force_read       = 0;
    } else if (strcmp(var, "ext_force") == 0 && count == 3) {
      ext_force[0]         = val1;
      ext_force[1]         = val2;
      dim_ext_force        = 2;
      ext_force_read       = 0;
    } else if (strcmp(var, "freq_micro") == 0 && count >= 2) {
      freq_micro           = val1;
      freq_micro_read      = 0;
    } else if (strcmp(var, "freq_macro") == 0 && count >= 2) {
      freq_macro           = val1;
      freq_macro_read      = 0;
    } else if (strcmp(var, "freq_walls") == 0 && count >= 2) {
      freq_walls           = val1;
      freq_walls_read      = 0;
    } else if (strcmp(var, "freq_colloids") == 0 && count >= 2) {
      freq_colloids        = val1;
      freq_colloids_read   = 0;      
    } else if (strcmp(var, "wall") == 0 && count >= 2) {
      wall                 = val1;
      wall_read            = 0;
    } else if (strcmp(var, "Vwall_top") == 0 && count >= 2) {
      Vwall_top            = val1;
      Vwall_top_read       = 0;
    } else if (strcmp(var, "Vwall_bottom") == 0 && count >= 2) {
      Vwall_bottom         = val1;
      Vwall_bottom_read    = 0;
    } else if (strcmp(var, "N_colloids") == 0 && count >= 2) {
      N_coll               = val1;
      N_coll_read          = 0;
    } else if (strcmp(var, "coll_R") == 0 && count >= 2) {
      coll_radius          = val1;
      coll_radius_read     = 0;
    } else if (strcmp(var, "coll_rho") == 0 && count >= 2) {
      coll_density         = val1;
      coll_density_read    = 0;
    } else if (strcmp(var, "coll_move") == 0 && count >= 2) {
      coll_move            = val1;
      coll_move_read       = 0;
    } else if (strcmp(var, "coll_repulsion_cuton") == 0 && count >= 2) {
      coll_repulsion_cuton      = val1;
      coll_repulsion_cuton_read = 0;
    } else if (strcmp(var, "F0_repulsion") == 0 && count >= 2) {
      F0_repulsion         = val1;
      F0_repulsion_read    = 0;
    } else if (strcmp(var, "tau_repulsion") == 0 && count >= 2) {
      tau_repulsion        = val1;
      tau_repulsion_read   = 0;
    } else if (strcmp(var, "new_sim") == 0 && count >= 2) {
      new_sim              = val1;
      new_sim_read         = 0;            
     }
    else if (strcmp(var, "coll_x") == 0){
      if (count >= 3 && N_coll > 0) {
	if (coll_x_count >= N_coll) {
	  printf("system read input error: different coll_x entries than N_colloids\n");
	  return 1;
	}
	char* token = strtok(start, " "); // The line is divided in tokens
	int value_index = 0;
	while (token != nullptr) {
	  if (strcmp(token, "coll_x") != 0) { // If the token is not equal to "coll_x"
	    coll_x_data[coll_x_count][value_index] = atof(token); //atof converts the string into a number
	    value_index++;
	  }
	  token = strtok(nullptr, " "); // The next token is obtained.
	}
	if (value_index != dim) {
	  printf("system_read_input error: the number of dimensions of coll_x "
		 "is not compatible with the number of dimensions of the system\n");
	  printf("%d\n", value_index);
	  return 1;
	}
	coll_x_count++;      
      }      
    }  
  }

  // Checking that every variable was read
  if (N_read == 1)  {
    printf("system_read_input error: N was not read\n");
    return 1;
  }
  if (L_read == 1)  {
    printf("system_read_input error: L was not read\n");
    return 1;
  }
  if (dim_read == 1)  {
    printf("system_read_input error: dim was not read\n");
    return 1;
  }
  if (dt_read == 1)  {
    printf("system_read_input error: dt was not read\n");
    return 1;
  }
  if (Nsteps_read == 1)  {
    printf("system_read_input error: Nsteps was not read\n");
    return 1;
  }
  if (overlap_read == 1)  {
    printf("system_read_input error: overlap was not read\n");
    return 1;
  }
  if (rho_read == 1)  {
    printf("system_read_input error: rho was not read\n");
    return 1;
  }
  if (c_read == 1)  {
    printf("system_read_input error: c was not read\n");
    return 1;
  }
  if (P0_read == 1)  {
    printf("system_read_input error: P0 was not read\n");
    return 1;
  }
  if (eta_read == 1)  {
    printf("system_read_input error: eta was not read\n");
    return 1;
  }
  if (zeta_read == 1)  {
    printf("system_read_input error: zeta was not read\n");
    return 1;
  }
  if (ext_force_type_read == 1)  {
    printf("system_read_input error: ext_force_type was not read\n");
    return 1;
  }       
  if (ext_force_read == 1)  {
    printf("system_read_input error: ext_force was not read\n");
    return 1;
  }
  if (freq_micro_read == 1)  {
    printf("system_read_input error: freq_micro was not read\n");
    return 1;
  }
  if (freq_macro_read == 1)  {
    printf("system_read_input error: freq_macro was not read\n");
    return 1;
  } 
  if (freq_walls_read == 1)  {
    printf("system_read_input error: freq_walls was not read\n");
    return 1;
  }
  if (freq_colloids_read == 1)  {
    printf("system_read_input error: freq_colloids was not read\n");
    return 1;
  }     
  if (wall_read == 1)  {
    printf("system_read_input error: wall was not read\n");
    return 1;
  }
  if (Vwall_top_read == 1)  {
    printf("system_read_input error: Vwall_top was not read\n");
    return 1;
  }
  if (Vwall_bottom_read == 1)  {
    printf("system_read_input error: Vwall_bottom was not read\n");
    return 1;
  }
  if (N_coll_read == 1)  {
    printf("system_read_input error: N_colloids was not read\n");
    return 1;
  }
  if (coll_radius_read == 1)  {
    printf("system_read_input error: coll_R was not read\n");
    return 1;
  }
  if (coll_density_read == 1)  {
    printf("system_read_input error: coll_rho was not read\n");
    return 1;
  }
  if (coll_move_read == 1)  {
    printf("system_read_input error: coll_move was not read\n");
    return 1;
  }
  if (coll_repulsion_cuton_read == 1)  {
    printf("system_read_input error: coll_repulsion_cuton was not read\n");
    return 1;
  }
  if (F0_repulsion_read == 1)  {
    printf("system_read_input error: F0_repulsion was not read\n");
    return 1;
  }
  if (tau_repulsion_read == 1)  {
    printf("system_read_input error: tau_repulsion was not read\n");
    return 1;
  }
  if (new_sim_read == 1)  {
    printf("system_read_input error: new_sim was not read\n");
    return 1;
  }    
  if (coll_x_count != N_coll) {
    printf("system read input error: different number of coll_x variables than N_colloids\n");
    printf("coll_x_count = %d, N_colloids = %d\n", coll_x_count, N_coll);
    return 1;
  }  

  if (dim != dim_N)  {
    printf("system_read_input error: the number of dimensions of the system is not compatible with the number of dimensions of N\n");
    return 1;
  }
  if (dim != dim_L)  {
    printf("system_read_input error: the number of dimensions of the system is not compatible with the number of dimensions of L\n");
    return 1;
  }
  if (ext_force_type != 0) 
    if (dim != dim_ext_force)  {
      printf("system_read_input error: the number of dimensions of the system is not compatible with the number of dimensions of ext_force\n");
      return 1;
    }
  if (ext_force_type < 0 || ext_force_type > 3)  {
    printf("system_read_input error: the option ext_force_type = %d is not available.\n", ext_force_type);
    return 1;
  }
  if (N_coll > MAX_COLLOIDS) {
    printf("System read input error: N_colloids is too high. Consider to change the value of MAX_COLLOIDS\n");
    return 1;
  }

  //---- Checking if colloids are inside the simulation box ----
  if (N_coll >= 2)
    for (int i = 0; i < N_coll; ++i) {
      if (coll_x_data[i][0] < 0 || coll_x_data[i][0] > L[0]) {
	printf("System read input error: Colloid is out of the simulation box\n");
	return 1;	
      }
      if (coll_x_data[i][1] < 0 || coll_x_data[i][1] > L[1]) {
	printf("System read input error: Colloid is out of the simulation box\n");
	return 1;	
      }
      if (dim == 3)
	if (coll_x_data[i][1] < 0 || coll_x_data[i][1] > L[2]) {
	  printf("System read input error: Colloid is out of the simulation box\n");
	  return 1;	
	}	
    }

  //---- Checking if colloids are overlapping -------
  if (N_coll >= 2)
    for (int i = 0; i < N_coll; ++i)
      for (int j = i+1; j < N_coll; ++j) 
	if (dim == 2) {
	  real xij = coll_x_data[i][0] - coll_x_data[j][0];
	  real yij = coll_x_data[i][1] - coll_x_data[j][1];
	  //--- Periodic boundary conditions ---
	  if (xij > 0.5 * L[0])
	    xij -= L[0];
	  if (xij < -0.5 * L[0])
	    xij += L[0];
	  if (wall == 0) {   //If there are not walls	      	    
	    if (yij > 0.5 * L[1])
	      yij -= L[1];
	    if (yij < -0.5 * L[1])
	      yij += L[1];
	  }
	  real rij = sqrt(xij*xij + yij*yij);
	  if (rij < 2 * coll_radius) {
	    printf("System read input error: Overlapping colloids\n");
	    return 1;
	  }
	}
	else {  // dim == 3
	  real xij = coll_x_data[i][0] - coll_x_data[j][0];
	  real yij = coll_x_data[i][1] - coll_x_data[j][1];
	  real zij = coll_x_data[i][1] - coll_x_data[j][1];	
	  //--- Periodic boundary conditions ---
	  if (xij > 0.5 * L[0])
	    xij -= L[0];
	  if (xij > -0.5 * L[0])
	    xij += L[0];
	  if (wall == 0) {   //If there are not walls	      	    
	    if (yij > 0.5 * L[1])
	      yij -= L[1];
	    if (yij < -0.5 * L[1])
	      yij += L[1];
	  }
	  if (zij > 0.5 * L[2])
	    zij -= L[2];
	  if (zij > -0.5 * L[2])
	    zij += L[2];	  
	
	  real rij = sqrt(xij*xij + yij*yij + zij*zij);
	  if (rij < 2 * coll_radius) {
	    printf("System read input error: Overlapping colloids\n");
	    return 1;
	  }
	}      

  fclose(file);

  return 0;
}
