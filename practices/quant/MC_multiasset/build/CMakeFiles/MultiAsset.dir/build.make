# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.12.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.12.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jungjaeyong/projects/practices/quant/MC_multiasset

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jungjaeyong/projects/practices/quant/MC_multiasset/build

# Include any dependencies generated for this target.
include CMakeFiles/MultiAsset.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MultiAsset.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MultiAsset.dir/flags.make

CMakeFiles/MultiAsset.dir/main.cpp.o: CMakeFiles/MultiAsset.dir/flags.make
CMakeFiles/MultiAsset.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jungjaeyong/projects/practices/quant/MC_multiasset/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MultiAsset.dir/main.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MultiAsset.dir/main.cpp.o -c /Users/jungjaeyong/projects/practices/quant/MC_multiasset/main.cpp

CMakeFiles/MultiAsset.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MultiAsset.dir/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jungjaeyong/projects/practices/quant/MC_multiasset/main.cpp > CMakeFiles/MultiAsset.dir/main.cpp.i

CMakeFiles/MultiAsset.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MultiAsset.dir/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jungjaeyong/projects/practices/quant/MC_multiasset/main.cpp -o CMakeFiles/MultiAsset.dir/main.cpp.s

CMakeFiles/MultiAsset.dir/utils.cpp.o: CMakeFiles/MultiAsset.dir/flags.make
CMakeFiles/MultiAsset.dir/utils.cpp.o: ../utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jungjaeyong/projects/practices/quant/MC_multiasset/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MultiAsset.dir/utils.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MultiAsset.dir/utils.cpp.o -c /Users/jungjaeyong/projects/practices/quant/MC_multiasset/utils.cpp

CMakeFiles/MultiAsset.dir/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MultiAsset.dir/utils.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jungjaeyong/projects/practices/quant/MC_multiasset/utils.cpp > CMakeFiles/MultiAsset.dir/utils.cpp.i

CMakeFiles/MultiAsset.dir/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MultiAsset.dir/utils.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jungjaeyong/projects/practices/quant/MC_multiasset/utils.cpp -o CMakeFiles/MultiAsset.dir/utils.cpp.s

CMakeFiles/MultiAsset.dir/models.cpp.o: CMakeFiles/MultiAsset.dir/flags.make
CMakeFiles/MultiAsset.dir/models.cpp.o: ../models.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jungjaeyong/projects/practices/quant/MC_multiasset/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MultiAsset.dir/models.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MultiAsset.dir/models.cpp.o -c /Users/jungjaeyong/projects/practices/quant/MC_multiasset/models.cpp

CMakeFiles/MultiAsset.dir/models.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MultiAsset.dir/models.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jungjaeyong/projects/practices/quant/MC_multiasset/models.cpp > CMakeFiles/MultiAsset.dir/models.cpp.i

CMakeFiles/MultiAsset.dir/models.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MultiAsset.dir/models.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jungjaeyong/projects/practices/quant/MC_multiasset/models.cpp -o CMakeFiles/MultiAsset.dir/models.cpp.s

# Object files for target MultiAsset
MultiAsset_OBJECTS = \
"CMakeFiles/MultiAsset.dir/main.cpp.o" \
"CMakeFiles/MultiAsset.dir/utils.cpp.o" \
"CMakeFiles/MultiAsset.dir/models.cpp.o"

# External object files for target MultiAsset
MultiAsset_EXTERNAL_OBJECTS =

MultiAsset: CMakeFiles/MultiAsset.dir/main.cpp.o
MultiAsset: CMakeFiles/MultiAsset.dir/utils.cpp.o
MultiAsset: CMakeFiles/MultiAsset.dir/models.cpp.o
MultiAsset: CMakeFiles/MultiAsset.dir/build.make
MultiAsset: CMakeFiles/MultiAsset.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/jungjaeyong/projects/practices/quant/MC_multiasset/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable MultiAsset"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MultiAsset.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MultiAsset.dir/build: MultiAsset

.PHONY : CMakeFiles/MultiAsset.dir/build

CMakeFiles/MultiAsset.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MultiAsset.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MultiAsset.dir/clean

CMakeFiles/MultiAsset.dir/depend:
	cd /Users/jungjaeyong/projects/practices/quant/MC_multiasset/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jungjaeyong/projects/practices/quant/MC_multiasset /Users/jungjaeyong/projects/practices/quant/MC_multiasset /Users/jungjaeyong/projects/practices/quant/MC_multiasset/build /Users/jungjaeyong/projects/practices/quant/MC_multiasset/build /Users/jungjaeyong/projects/practices/quant/MC_multiasset/build/CMakeFiles/MultiAsset.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MultiAsset.dir/depend

