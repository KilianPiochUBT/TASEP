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
include failtest/CMakeFiles/ref_4_ok.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include failtest/CMakeFiles/ref_4_ok.dir/compiler_depend.make

# Include the progress variables for this target.
include failtest/CMakeFiles/ref_4_ok.dir/progress.make

# Include the compile flags for this target's objects.
include failtest/CMakeFiles/ref_4_ok.dir/flags.make

failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.o: failtest/CMakeFiles/ref_4_ok.dir/flags.make
failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.o: ../failtest/ref_4.cpp
failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.o: failtest/CMakeFiles/ref_4_ok.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.o"
	cd /home/bt308081/build/eigen-3.4.0/build/failtest && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.o -MF CMakeFiles/ref_4_ok.dir/ref_4.cpp.o.d -o CMakeFiles/ref_4_ok.dir/ref_4.cpp.o -c /home/bt308081/build/eigen-3.4.0/failtest/ref_4.cpp

failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ref_4_ok.dir/ref_4.cpp.i"
	cd /home/bt308081/build/eigen-3.4.0/build/failtest && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bt308081/build/eigen-3.4.0/failtest/ref_4.cpp > CMakeFiles/ref_4_ok.dir/ref_4.cpp.i

failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ref_4_ok.dir/ref_4.cpp.s"
	cd /home/bt308081/build/eigen-3.4.0/build/failtest && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bt308081/build/eigen-3.4.0/failtest/ref_4.cpp -o CMakeFiles/ref_4_ok.dir/ref_4.cpp.s

# Object files for target ref_4_ok
ref_4_ok_OBJECTS = \
"CMakeFiles/ref_4_ok.dir/ref_4.cpp.o"

# External object files for target ref_4_ok
ref_4_ok_EXTERNAL_OBJECTS =

failtest/ref_4_ok: failtest/CMakeFiles/ref_4_ok.dir/ref_4.cpp.o
failtest/ref_4_ok: failtest/CMakeFiles/ref_4_ok.dir/build.make
failtest/ref_4_ok: failtest/CMakeFiles/ref_4_ok.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ref_4_ok"
	cd /home/bt308081/build/eigen-3.4.0/build/failtest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ref_4_ok.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
failtest/CMakeFiles/ref_4_ok.dir/build: failtest/ref_4_ok
.PHONY : failtest/CMakeFiles/ref_4_ok.dir/build

failtest/CMakeFiles/ref_4_ok.dir/clean:
	cd /home/bt308081/build/eigen-3.4.0/build/failtest && $(CMAKE_COMMAND) -P CMakeFiles/ref_4_ok.dir/cmake_clean.cmake
.PHONY : failtest/CMakeFiles/ref_4_ok.dir/clean

failtest/CMakeFiles/ref_4_ok.dir/depend:
	cd /home/bt308081/build/eigen-3.4.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bt308081/build/eigen-3.4.0 /home/bt308081/build/eigen-3.4.0/failtest /home/bt308081/build/eigen-3.4.0/build /home/bt308081/build/eigen-3.4.0/build/failtest /home/bt308081/build/eigen-3.4.0/build/failtest/CMakeFiles/ref_4_ok.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : failtest/CMakeFiles/ref_4_ok.dir/depend

