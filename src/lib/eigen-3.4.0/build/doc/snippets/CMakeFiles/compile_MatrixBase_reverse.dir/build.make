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
include doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compiler_depend.make

# Include the progress variables for this target.
include doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/progress.make

# Include the compile flags for this target's objects.
include doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/flags.make

doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o: doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/flags.make
doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o: doc/snippets/compile_MatrixBase_reverse.cpp
doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o: ../doc/snippets/MatrixBase_reverse.cpp
doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o: doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o"
	cd /home/bt308081/build/eigen-3.4.0/build/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o -MF CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o.d -o CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o -c /home/bt308081/build/eigen-3.4.0/build/doc/snippets/compile_MatrixBase_reverse.cpp

doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.i"
	cd /home/bt308081/build/eigen-3.4.0/build/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bt308081/build/eigen-3.4.0/build/doc/snippets/compile_MatrixBase_reverse.cpp > CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.i

doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.s"
	cd /home/bt308081/build/eigen-3.4.0/build/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bt308081/build/eigen-3.4.0/build/doc/snippets/compile_MatrixBase_reverse.cpp -o CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.s

# Object files for target compile_MatrixBase_reverse
compile_MatrixBase_reverse_OBJECTS = \
"CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o"

# External object files for target compile_MatrixBase_reverse
compile_MatrixBase_reverse_EXTERNAL_OBJECTS =

doc/snippets/compile_MatrixBase_reverse: doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/compile_MatrixBase_reverse.cpp.o
doc/snippets/compile_MatrixBase_reverse: doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/build.make
doc/snippets/compile_MatrixBase_reverse: doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bt308081/build/eigen-3.4.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_MatrixBase_reverse"
	cd /home/bt308081/build/eigen-3.4.0/build/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_MatrixBase_reverse.dir/link.txt --verbose=$(VERBOSE)
	cd /home/bt308081/build/eigen-3.4.0/build/doc/snippets && ./compile_MatrixBase_reverse >/home/bt308081/build/eigen-3.4.0/build/doc/snippets/MatrixBase_reverse.out

# Rule to build all files generated by this target.
doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/build: doc/snippets/compile_MatrixBase_reverse
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/build

doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/clean:
	cd /home/bt308081/build/eigen-3.4.0/build/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_MatrixBase_reverse.dir/cmake_clean.cmake
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/clean

doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/depend:
	cd /home/bt308081/build/eigen-3.4.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bt308081/build/eigen-3.4.0 /home/bt308081/build/eigen-3.4.0/doc/snippets /home/bt308081/build/eigen-3.4.0/build /home/bt308081/build/eigen-3.4.0/build/doc/snippets /home/bt308081/build/eigen-3.4.0/build/doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_reverse.dir/depend
