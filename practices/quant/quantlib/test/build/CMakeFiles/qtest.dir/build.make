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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.22.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.22.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jungjaeyong/projects/practices/quant/quantlib/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jungjaeyong/projects/practices/quant/quantlib/test/build

# Include any dependencies generated for this target.
include CMakeFiles/qtest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/qtest.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/qtest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/qtest.dir/flags.make

CMakeFiles/qtest.dir/main.cpp.o: CMakeFiles/qtest.dir/flags.make
CMakeFiles/qtest.dir/main.cpp.o: ../main.cpp
CMakeFiles/qtest.dir/main.cpp.o: CMakeFiles/qtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jungjaeyong/projects/practices/quant/quantlib/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/qtest.dir/main.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/qtest.dir/main.cpp.o -MF CMakeFiles/qtest.dir/main.cpp.o.d -o CMakeFiles/qtest.dir/main.cpp.o -c /Users/jungjaeyong/projects/practices/quant/quantlib/test/main.cpp

CMakeFiles/qtest.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/qtest.dir/main.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jungjaeyong/projects/practices/quant/quantlib/test/main.cpp > CMakeFiles/qtest.dir/main.cpp.i

CMakeFiles/qtest.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/qtest.dir/main.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jungjaeyong/projects/practices/quant/quantlib/test/main.cpp -o CMakeFiles/qtest.dir/main.cpp.s

# Object files for target qtest
qtest_OBJECTS = \
"CMakeFiles/qtest.dir/main.cpp.o"

# External object files for target qtest
qtest_EXTERNAL_OBJECTS =

qtest: CMakeFiles/qtest.dir/main.cpp.o
qtest: CMakeFiles/qtest.dir/build.make
qtest: /usr/local/lib/libQuantLib.dylib
qtest: CMakeFiles/qtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/jungjaeyong/projects/practices/quant/quantlib/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable qtest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/qtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/qtest.dir/build: qtest
.PHONY : CMakeFiles/qtest.dir/build

CMakeFiles/qtest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/qtest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/qtest.dir/clean

CMakeFiles/qtest.dir/depend:
	cd /Users/jungjaeyong/projects/practices/quant/quantlib/test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jungjaeyong/projects/practices/quant/quantlib/test /Users/jungjaeyong/projects/practices/quant/quantlib/test /Users/jungjaeyong/projects/practices/quant/quantlib/test/build /Users/jungjaeyong/projects/practices/quant/quantlib/test/build /Users/jungjaeyong/projects/practices/quant/quantlib/test/build/CMakeFiles/qtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/qtest.dir/depend

