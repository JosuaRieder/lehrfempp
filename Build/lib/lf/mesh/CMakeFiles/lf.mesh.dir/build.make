# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_SOURCE_DIR = /u/magina/Documents/lehrfempp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/magina/Documents/lehrfempp/Build

# Include any dependencies generated for this target.
include lib/lf/mesh/CMakeFiles/lf.mesh.dir/depend.make

# Include the progress variables for this target.
include lib/lf/mesh/CMakeFiles/lf.mesh.dir/progress.make

# Include the compile flags for this target's objects.
include lib/lf/mesh/CMakeFiles/lf.mesh.dir/flags.make

lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh.cc.o: lib/lf/mesh/CMakeFiles/lf.mesh.dir/flags.make
lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh.cc.o: ../lib/lf/mesh/mesh.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/magina/Documents/lehrfempp/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh.cc.o"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.dir/mesh.cc.o -c /u/magina/Documents/lehrfempp/lib/lf/mesh/mesh.cc

lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.dir/mesh.cc.i"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/magina/Documents/lehrfempp/lib/lf/mesh/mesh.cc > CMakeFiles/lf.mesh.dir/mesh.cc.i

lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.dir/mesh.cc.s"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/magina/Documents/lehrfempp/lib/lf/mesh/mesh.cc -o CMakeFiles/lf.mesh.dir/mesh.cc.s

lib/lf/mesh/CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.o: lib/lf/mesh/CMakeFiles/lf.mesh.dir/flags.make
lib/lf/mesh/CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.o: ../lib/lf/mesh/tp_triag_mesh_builder.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/magina/Documents/lehrfempp/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object lib/lf/mesh/CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.o"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.o -c /u/magina/Documents/lehrfempp/lib/lf/mesh/tp_triag_mesh_builder.cc

lib/lf/mesh/CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.i"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/magina/Documents/lehrfempp/lib/lf/mesh/tp_triag_mesh_builder.cc > CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.i

lib/lf/mesh/CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.s"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/magina/Documents/lehrfempp/lib/lf/mesh/tp_triag_mesh_builder.cc -o CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.s

lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh_interface.cc.o: lib/lf/mesh/CMakeFiles/lf.mesh.dir/flags.make
lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh_interface.cc.o: ../lib/lf/mesh/mesh_interface.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/magina/Documents/lehrfempp/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh_interface.cc.o"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.dir/mesh_interface.cc.o -c /u/magina/Documents/lehrfempp/lib/lf/mesh/mesh_interface.cc

lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh_interface.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.dir/mesh_interface.cc.i"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/magina/Documents/lehrfempp/lib/lf/mesh/mesh_interface.cc > CMakeFiles/lf.mesh.dir/mesh_interface.cc.i

lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh_interface.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.dir/mesh_interface.cc.s"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/magina/Documents/lehrfempp/lib/lf/mesh/mesh_interface.cc -o CMakeFiles/lf.mesh.dir/mesh_interface.cc.s

# Object files for target lf.mesh
lf_mesh_OBJECTS = \
"CMakeFiles/lf.mesh.dir/mesh.cc.o" \
"CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.o" \
"CMakeFiles/lf.mesh.dir/mesh_interface.cc.o"

# External object files for target lf.mesh
lf_mesh_EXTERNAL_OBJECTS =

lib/lf/mesh/liblf.mesh.a: lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh.cc.o
lib/lf/mesh/liblf.mesh.a: lib/lf/mesh/CMakeFiles/lf.mesh.dir/tp_triag_mesh_builder.cc.o
lib/lf/mesh/liblf.mesh.a: lib/lf/mesh/CMakeFiles/lf.mesh.dir/mesh_interface.cc.o
lib/lf/mesh/liblf.mesh.a: lib/lf/mesh/CMakeFiles/lf.mesh.dir/build.make
lib/lf/mesh/liblf.mesh.a: lib/lf/mesh/CMakeFiles/lf.mesh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/u/magina/Documents/lehrfempp/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library liblf.mesh.a"
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && $(CMAKE_COMMAND) -P CMakeFiles/lf.mesh.dir/cmake_clean_target.cmake
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lf.mesh.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/lf/mesh/CMakeFiles/lf.mesh.dir/build: lib/lf/mesh/liblf.mesh.a

.PHONY : lib/lf/mesh/CMakeFiles/lf.mesh.dir/build

lib/lf/mesh/CMakeFiles/lf.mesh.dir/clean:
	cd /u/magina/Documents/lehrfempp/Build/lib/lf/mesh && $(CMAKE_COMMAND) -P CMakeFiles/lf.mesh.dir/cmake_clean.cmake
.PHONY : lib/lf/mesh/CMakeFiles/lf.mesh.dir/clean

lib/lf/mesh/CMakeFiles/lf.mesh.dir/depend:
	cd /u/magina/Documents/lehrfempp/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/magina/Documents/lehrfempp /u/magina/Documents/lehrfempp/lib/lf/mesh /u/magina/Documents/lehrfempp/Build /u/magina/Documents/lehrfempp/Build/lib/lf/mesh /u/magina/Documents/lehrfempp/Build/lib/lf/mesh/CMakeFiles/lf.mesh.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/lf/mesh/CMakeFiles/lf.mesh.dir/depend

