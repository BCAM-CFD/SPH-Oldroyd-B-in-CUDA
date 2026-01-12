NVCC = nvcc
CFLAGS = -O3 --std=c++14 -arch=sm_61 --relocatable-device-code=true # sm_86 depends on your GPU
# -g -G
#beta: c-1-6 -> sm_86
#beta: c-1-2 -> sm_61
#dseta: c-4-0 -> sm_61
#dseta: c-4-2 -> sm_86
#portatil: sm_75

# Source files
HEADERS = config.h class_system.h kernel_functions.h
CUDA_SOURCES = main.cu kernel_print_particles.cu kernel_print_constants.cu kernel_declare_constants.cu kernel_forces_to_zero.cu kernel_density_to_zero.cu kernel_forces.cu kernel_ext_force.cu kernel_move_sys_VV_part1.cu kernel_move_sys_VV_part2.cu kernel_move_sys_Euler.cu kernel_W.cu kernel_gradW.cu kernel_densities.cu kernel_pressures.cu kernel_print_cell_list.cu kernel_cell_list.cu kernel_cell_ranges.cu kernel_print_cell_ranges.cu kernel_calculate_macro_vars.cu kernel_macro_vars_to_zero.cu kernel_collect_force_walls.cu kernel_collect_force_colloids.cu kernel_move_colloids_VV_part1.cu kernel_move_colloids_VV_part2.cu kernel_cell_colloids_list.cu kernel_initialize_cell_colloids.cu kernel_cell_colloids_ranges.cu kernel_colloid_forces_to_zero.cu kernel_repulsion_force_colloids.cu kernel_quaternion_product.cu kernel_vel_gradients_to_zero.cu kernel_vel_gradients.cu
CPP_SYSTEM_SOURCES  = system_read_input.cpp system_initialize_pointers.cpp system_constructor.cpp system_initialize.cpp system_print_info.cpp system_print_particles.cpp system_copy_pointers_to_device.cpp system_destructor.cpp system_neighbours_search.cu system_calculate_density.cu system_calculate_pressure.cu system_calculate_forces.cu system_calculate_external_forces.cu system_V_Verlet_step1.cu system_V_Verlet_step2.cu system_print_macro_vars.cu system_print_output.cpp system_calculate_force_walls.cu system_print_walls.cpp system_calculate_force_colloids.cu system_print_colloids.cpp system_initialize_particles_new_sim.cpp system_initialize_particles_old_sim.cpp system_calculate_vel_gradient.cu kernel_dcdt.cu system_print_particles_gradv.cpp
CPP_CONFIG_SOURCES  = config.h

# Objects are created depparately for CUDA and C++
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=.o)
CPP_SYSTEM_OBJECTS  = $(CPP_SYSTEM_SOURCES:.cpp=.o)

EXECUTABLE = SDPD_program

all: $(EXECUTABLE)

$(EXECUTABLE): $(CUDA_OBJECTS) $(CPP_SYSTEM_OBJECTS)
	$(NVCC) $(CFLAGS) -o $@ $^

%.o: %.cu $(HEADERS)
	$(NVCC) $(CFLAGS) -c $< -o $@

%.o: %.cpp $(HEADERS)
	$(NVCC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)
