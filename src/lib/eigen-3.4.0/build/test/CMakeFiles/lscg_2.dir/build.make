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
include test/CMakeFiles/lscg_2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/lscg_2.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/lscg_2.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/lscg_2.dir/flags.make

test/CMakeFiles/lscg_2.dir/lscg.cpp.o: test/CMakeFiles/lscg_2.dir/flags.make
test/CMakeFiles/lscg_2.dir/lscg.cpp.o: ../test/lscg.cpp
test/CMakeFiles/lscg_2.dir/lscg.cpp.o: test/CMakeFiles/lscg_2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/lscg_2.dir/lscg.cpp.o"
	cd /home/bt308081/build/eigen-3.4.0/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/lscg_2.dir/lscg.cpp.o -MF CMakeFiles/lscg_2.dir/lscg.cpp.o.d -o CMakeFiles/lscg_2.dir/lscg.cpp.o -c /home/bt308081/build/eigen-3.4.0/test/lscg.cpp

test/CMakeFiles/lscg_2.dir/lscg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscg_2.dir/lscg.cpp.i"
	cd /home/bt308081/build/eigen-3.4.0/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bt308081/build/eigen-3.4.0/test/lscg.cpp > CMakeFiles/lscg_2.dir/lscg.cpp.i

test/CMakeFiles/lscg_2.dir/lscg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscg_2.dir/lscg.cpp.s"
	cd /home/bt308081/build/eigen-3.4.0/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bt308081/build/eigen-3.4.0/test/lscg.cpp -o CMakeFiles/lscg_2.dir/lscg.cpp.s

# Object files for target lscg_2
lscg_2_OBJECTS = \
"CMakeFiles/lscg_2.dir/lscg.cpp.o"

# External object files for target lscg_2
lscg_2_EXTERNAL_OBJECTS =

test/lscg_2: test/CMakeFiles/lscg_2.dir/lscg.cpp.o
test/lscg_2: test/CMakeFiles/lscg_2.dir/build.make
test/lscg_2: test/CMakeFiles/lscg_2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lscg_2"
	cd /home/bt308081/build/eigen-3.4.0/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lscg_2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/lscg_2.dir/build: test/lscg_2
.PHONY : test/CMakeFiles/lscg_2.dir/build

test/CMakeFiles/lscg_2.dir/clean:
	cd /home/bt308081/build/eigen-3.4.0/build/test && $(CMAKE_COMMAND) -P CMakeFiles/lscg_2.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/lscg_2.dir/clean

test/CMakeFiles/lscg_2.dir/depend:
	cd /home/bt308081/build/eigen-3.4.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bt308081/build/eigen-3.4.0 /home/bt308081/build/eigen-3.4.0/test /home/bt308081/build/eigen-3.4.0/build /home/bt308081/build/eigen-3.4.0/build/test /home/bt308081/build/eigen-3.4.0/build/test/CMakeFiles/lscg_2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/lscg_2.dir/depend
