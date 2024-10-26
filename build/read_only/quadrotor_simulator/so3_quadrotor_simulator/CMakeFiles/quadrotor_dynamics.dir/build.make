# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tyh/DB_plan_Project/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tyh/DB_plan_Project/build

# Include any dependencies generated for this target.
include read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/depend.make

# Include the progress variables for this target.
include read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/progress.make

# Include the compile flags for this target's objects.
include read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/flags.make

read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.o: read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/flags.make
read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.o: /home/tyh/DB_plan_Project/src/read_only/quadrotor_simulator/so3_quadrotor_simulator/src/dynamics/Quadrotor.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tyh/DB_plan_Project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.o"
	cd /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.o -c /home/tyh/DB_plan_Project/src/read_only/quadrotor_simulator/so3_quadrotor_simulator/src/dynamics/Quadrotor.cpp

read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.i"
	cd /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tyh/DB_plan_Project/src/read_only/quadrotor_simulator/so3_quadrotor_simulator/src/dynamics/Quadrotor.cpp > CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.i

read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.s"
	cd /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tyh/DB_plan_Project/src/read_only/quadrotor_simulator/so3_quadrotor_simulator/src/dynamics/Quadrotor.cpp -o CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.s

# Object files for target quadrotor_dynamics
quadrotor_dynamics_OBJECTS = \
"CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.o"

# External object files for target quadrotor_dynamics
quadrotor_dynamics_EXTERNAL_OBJECTS =

/home/tyh/DB_plan_Project/devel/lib/libquadrotor_dynamics.so: read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/src/dynamics/Quadrotor.cpp.o
/home/tyh/DB_plan_Project/devel/lib/libquadrotor_dynamics.so: read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/build.make
/home/tyh/DB_plan_Project/devel/lib/libquadrotor_dynamics.so: read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/tyh/DB_plan_Project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library /home/tyh/DB_plan_Project/devel/lib/libquadrotor_dynamics.so"
	cd /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quadrotor_dynamics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/build: /home/tyh/DB_plan_Project/devel/lib/libquadrotor_dynamics.so

.PHONY : read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/build

read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/clean:
	cd /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator && $(CMAKE_COMMAND) -P CMakeFiles/quadrotor_dynamics.dir/cmake_clean.cmake
.PHONY : read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/clean

read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/depend:
	cd /home/tyh/DB_plan_Project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tyh/DB_plan_Project/src /home/tyh/DB_plan_Project/src/read_only/quadrotor_simulator/so3_quadrotor_simulator /home/tyh/DB_plan_Project/build /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator /home/tyh/DB_plan_Project/build/read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : read_only/quadrotor_simulator/so3_quadrotor_simulator/CMakeFiles/quadrotor_dynamics.dir/depend
