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
CMAKE_SOURCE_DIR = /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build

# Include any dependencies generated for this target.
include CMakeFiles/example.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/example.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/example.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example.dir/flags.make

CMakeFiles/example.dir/main.cpp.o: CMakeFiles/example.dir/flags.make
CMakeFiles/example.dir/main.cpp.o: ../main.cpp
CMakeFiles/example.dir/main.cpp.o: CMakeFiles/example.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/example.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/example.dir/main.cpp.o -MF CMakeFiles/example.dir/main.cpp.o.d -o CMakeFiles/example.dir/main.cpp.o -c /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/main.cpp

CMakeFiles/example.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/example.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/main.cpp > CMakeFiles/example.dir/main.cpp.i

CMakeFiles/example.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/example.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/main.cpp -o CMakeFiles/example.dir/main.cpp.s

# Object files for target example
example_OBJECTS = \
"CMakeFiles/example.dir/main.cpp.o"

# External object files for target example
example_EXTERNAL_OBJECTS =

example: CMakeFiles/example.dir/main.cpp.o
example: CMakeFiles/example.dir/build.make
example: /Users/jungjaeyong/opt/anaconda3/lib/python3.7/site-packages/torch/lib/libc10.dylib
example: /usr/local/lib/libkineto.a
example: /Users/jungjaeyong/opt/anaconda3/lib/python3.7/site-packages/torch/lib/libtorch.dylib
example: /Users/jungjaeyong/opt/anaconda3/lib/python3.7/site-packages/torch/lib/libtorch_cpu.dylib
example: /Users/jungjaeyong/opt/anaconda3/lib/python3.7/site-packages/torch/lib/libc10.dylib
example: CMakeFiles/example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable example"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example.dir/build: example
.PHONY : CMakeFiles/example.dir/build

CMakeFiles/example.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example.dir/clean

CMakeFiles/example.dir/depend:
	cd /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build /Users/jungjaeyong/projects/practices/C++/libtorch_study/libtorch_example/build/CMakeFiles/example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/example.dir/depend

