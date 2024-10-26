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
CMAKE_BINARY_DIR = /home/tyh/DB_plan_Project/src/armadillo-9.870.2

# Utility rule file for quadrotor_msgs_generate_messages_py.

# Include the progress variables for this target.
include read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/progress.make

read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_AuxCommand.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Corrections.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Gains.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PPROutputData.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Serial.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_StatusData.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PolynomialTrajectory.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_LQRTrajectory.py
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py


devel/lib/python3/dist-packages/quadrotor_msgs/msg/_AuxCommand.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_AuxCommand.py: ../read_only/Utils/quadrotor_msgs/msg/AuxCommand.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating Python from MSG quadrotor_msgs/AuxCommand"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/AuxCommand.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Corrections.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Corrections.py: ../read_only/Utils/quadrotor_msgs/msg/Corrections.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating Python from MSG quadrotor_msgs/Corrections"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/Corrections.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Gains.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Gains.py: ../read_only/Utils/quadrotor_msgs/msg/Gains.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating Python from MSG quadrotor_msgs/Gains"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/Gains.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py: ../read_only/Utils/quadrotor_msgs/msg/OutputData.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Vector3.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Quaternion.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Generating Python from MSG quadrotor_msgs/OutputData"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/OutputData.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py: ../read_only/Utils/quadrotor_msgs/msg/PositionCommand.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Vector3.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Point.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Generating Python from MSG quadrotor_msgs/PositionCommand"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/PositionCommand.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PPROutputData.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PPROutputData.py: ../read_only/Utils/quadrotor_msgs/msg/PPROutputData.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PPROutputData.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Generating Python from MSG quadrotor_msgs/PPROutputData"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/PPROutputData.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Serial.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Serial.py: ../read_only/Utils/quadrotor_msgs/msg/Serial.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Serial.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Generating Python from MSG quadrotor_msgs/Serial"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/Serial.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py: ../read_only/Utils/quadrotor_msgs/msg/SO3Command.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py: ../read_only/Utils/quadrotor_msgs/msg/AuxCommand.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Vector3.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Quaternion.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Generating Python from MSG quadrotor_msgs/SO3Command"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/SO3Command.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_StatusData.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_StatusData.py: ../read_only/Utils/quadrotor_msgs/msg/StatusData.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_StatusData.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Generating Python from MSG quadrotor_msgs/StatusData"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/StatusData.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py: ../read_only/Utils/quadrotor_msgs/msg/TRPYCommand.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py: ../read_only/Utils/quadrotor_msgs/msg/AuxCommand.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Generating Python from MSG quadrotor_msgs/TRPYCommand"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/TRPYCommand.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: ../read_only/Utils/quadrotor_msgs/msg/Odometry.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Vector3.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Quaternion.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/PoseWithCovariance.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/TwistWithCovariance.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Twist.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /opt/ros/noetic/share/nav_msgs/msg/Odometry.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Point.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py: /home/tyh/me5400a_ws2/src/geometry_msgs/msg/Pose.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Generating Python from MSG quadrotor_msgs/Odometry"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/Odometry.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PolynomialTrajectory.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PolynomialTrajectory.py: ../read_only/Utils/quadrotor_msgs/msg/PolynomialTrajectory.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PolynomialTrajectory.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Generating Python from MSG quadrotor_msgs/PolynomialTrajectory"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/PolynomialTrajectory.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/_LQRTrajectory.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_LQRTrajectory.py: ../read_only/Utils/quadrotor_msgs/msg/LQRTrajectory.msg
devel/lib/python3/dist-packages/quadrotor_msgs/msg/_LQRTrajectory.py: /opt/ros/noetic/share/std_msgs/msg/Header.msg
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Generating Python from MSG quadrotor_msgs/LQRTrajectory"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg/LQRTrajectory.msg -Iquadrotor_msgs:/home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs/msg -Igeometry_msgs:/home/tyh/me5400a_ws2/src/geometry_msgs/msg -Inav_msgs:/opt/ros/noetic/share/nav_msgs/cmake/../msg -Istd_msgs:/opt/ros/noetic/share/std_msgs/cmake/../msg -Iactionlib_msgs:/opt/ros/noetic/share/actionlib_msgs/cmake/../msg -p quadrotor_msgs -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg

devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: /opt/ros/noetic/lib/genpy/genmsg_py.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_AuxCommand.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Corrections.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Gains.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PPROutputData.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Serial.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_StatusData.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PolynomialTrajectory.py
devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_LQRTrajectory.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tyh/DB_plan_Project/src/armadillo-9.870.2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Generating Python msg __init__.py for quadrotor_msgs"
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && ../../../catkin_generated/env_cached.sh /usr/bin/python3 /opt/ros/noetic/share/genpy/cmake/../../../lib/genpy/genmsg_py.py -o /home/tyh/DB_plan_Project/src/armadillo-9.870.2/devel/lib/python3/dist-packages/quadrotor_msgs/msg --initpy

quadrotor_msgs_generate_messages_py: read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_AuxCommand.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Corrections.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Gains.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_OutputData.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PositionCommand.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PPROutputData.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Serial.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_SO3Command.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_StatusData.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_TRPYCommand.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_Odometry.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_PolynomialTrajectory.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/_LQRTrajectory.py
quadrotor_msgs_generate_messages_py: devel/lib/python3/dist-packages/quadrotor_msgs/msg/__init__.py
quadrotor_msgs_generate_messages_py: read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/build.make

.PHONY : quadrotor_msgs_generate_messages_py

# Rule to build all files generated by this target.
read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/build: quadrotor_msgs_generate_messages_py

.PHONY : read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/build

read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/clean:
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs && $(CMAKE_COMMAND) -P CMakeFiles/quadrotor_msgs_generate_messages_py.dir/cmake_clean.cmake
.PHONY : read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/clean

read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/depend:
	cd /home/tyh/DB_plan_Project/src/armadillo-9.870.2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tyh/DB_plan_Project/src /home/tyh/DB_plan_Project/src/read_only/Utils/quadrotor_msgs /home/tyh/DB_plan_Project/src/armadillo-9.870.2 /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs /home/tyh/DB_plan_Project/src/armadillo-9.870.2/read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : read_only/Utils/quadrotor_msgs/CMakeFiles/quadrotor_msgs_generate_messages_py.dir/depend

