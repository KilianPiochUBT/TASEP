# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/bt308081/build/eigen-3.4.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bt308081/build/eigen-3.4.0/build

# Include any dependencies generated for this target.
include unsupported/test/CMakeFiles/autodiff_scalar_2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include unsupported/test/CMakeFiles/autodiff_scalar_2.dir/compiler_depend.make

# Include the progress variables for this target.
include unsupported/test/CMakeFiles/autodiff_scalar_2.dir/progress.make

# Include the compile flags for this target's objects.
include unsupported/test/CMakeFiles/autodiff_scalar_2.dir/flags.make

unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o: unsupported/test/CMakeFiles/autodiff_scalar_2.dir/flags.make
unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o: ../unsupported/test/autodiff_scalar.cpp
unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o: unsupported/test/CMakeFiles/autodiff_scalar_2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o"
	cd /home/bt308081/build/eigen-3.4.0/build/unsupported/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o -MF CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o.d -o CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o -c /home/bt308081/build/eigen-3.4.0/unsupported/test/autodiff_scalar.cpp

unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.i"
	cd /home/bt308081/build/eigen-3.4.0/build/unsupported/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bt308081/build/eigen-3.4.0/unsupported/test/autodiff_scalar.cpp > CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.i

unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.s"
	cd /home/bt308081/build/eigen-3.4.0/build/unsupported/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bt308081/build/eigen-3.4.0/unsupported/test/autodiff_scalar.cpp -o CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.s

# Object files for target autodiff_scalar_2
autodiff_scalar_2_OBJECTS = \
"CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o"

# External object files for target autodiff_scalar_2
autodiff_scalar_2_EXTERNAL_OBJECTS =

unsupported/test/autodiff_scalar_2: unsupported/test/CMakeFiles/autodiff_scalar_2.dir/autodiff_scalar.cpp.o
unsupported/test/autodiff_scalar_2: unsupported/test/CMakeFiles/autodiff_scalar_2.dir/build.make
unsupported/test/autodiff_scalar_2: unsupported/test/CMakeFiles/autodiff_scalar_2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable autodiff_scalar_2"
	cd /home/bt308081/build/eigen-3.4.0/build/unsupported/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/autodiff_scalar_2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unsupported/test/CMakeFiles/autodiff_scalar_2.dir/build: unsupported/test/autodiff_scalar_2
.PHONY : unsupported/test/CMakeFiles/autodiff_scalar_2.dir/build

unsupported/test/CMakeFiles/autodiff_scalar_2.dir/clean:
	cd /home/bt308081/build/eigen-3.4.0/build/unsupported/test && $(CMAKE_COMMAND) -P CMakeFiles/autodiff_scalar_2.dir/cmake_clean.cmake
.PHONY : unsupported/test/CMakeFiles/autodiff_scalar_2.dir/clean

unsupported/test/CMakeFiles/autodiff_scalar_2.dir/depend:
	cd /home/bt308081/build/eigen-3.4.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bt308081/build/eigen-3.4.0 /home/bt308081/build/eigen-3.4.0/unsupported/test /home/bt308081/build/eigen-3.4.0/build /home/bt308081/build/eigen-3.4.0/build/unsupported/test /home/bt308081/build/eigen-3.4.0/build/unsupported/test/CMakeFiles/autodiff_scalar_2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unsupported/test/CMakeFiles/autodiff_scalar_2.dir/depend

