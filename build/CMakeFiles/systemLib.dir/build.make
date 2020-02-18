# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /afs/crc.nd.edu/x86_64_linux/c/cmake/3.13.2/bin/cmake

# The command to remove a file.
RM = /afs/crc.nd.edu/x86_64_linux/c/cmake/3.13.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build

# Include any dependencies generated for this target.
include CMakeFiles/systemLib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/systemLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/systemLib.dir/flags.make

CMakeFiles/systemLib.dir/System.cu.o: CMakeFiles/systemLib.dir/flags.make
CMakeFiles/systemLib.dir/System.cu.o: /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src/System.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/systemLib.dir/System.cu.o"
	/afs/crc.nd.edu/x86_64_linux/c/cuda/10.0/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src/System.cu -o CMakeFiles/systemLib.dir/System.cu.o

CMakeFiles/systemLib.dir/System.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/systemLib.dir/System.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/systemLib.dir/System.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/systemLib.dir/System.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/systemLib.dir/SystemBuilder.cpp.o: CMakeFiles/systemLib.dir/flags.make
CMakeFiles/systemLib.dir/SystemBuilder.cpp.o: /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src/SystemBuilder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/systemLib.dir/SystemBuilder.cpp.o"
	/opt/crc/g/gcc/7.1.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/systemLib.dir/SystemBuilder.cpp.o -c /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src/SystemBuilder.cpp

CMakeFiles/systemLib.dir/SystemBuilder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/systemLib.dir/SystemBuilder.cpp.i"
	/opt/crc/g/gcc/7.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src/SystemBuilder.cpp > CMakeFiles/systemLib.dir/SystemBuilder.cpp.i

CMakeFiles/systemLib.dir/SystemBuilder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/systemLib.dir/SystemBuilder.cpp.s"
	/opt/crc/g/gcc/7.1.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src/SystemBuilder.cpp -o CMakeFiles/systemLib.dir/SystemBuilder.cpp.s

# Object files for target systemLib
systemLib_OBJECTS = \
"CMakeFiles/systemLib.dir/System.cu.o" \
"CMakeFiles/systemLib.dir/SystemBuilder.cpp.o"

# External object files for target systemLib
systemLib_EXTERNAL_OBJECTS =

libsystemLib.a: CMakeFiles/systemLib.dir/System.cu.o
libsystemLib.a: CMakeFiles/systemLib.dir/SystemBuilder.cpp.o
libsystemLib.a: CMakeFiles/systemLib.dir/build.make
libsystemLib.a: CMakeFiles/systemLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libsystemLib.a"
	$(CMAKE_COMMAND) -P CMakeFiles/systemLib.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/systemLib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/systemLib.dir/build: libsystemLib.a

.PHONY : CMakeFiles/systemLib.dir/build

CMakeFiles/systemLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/systemLib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/systemLib.dir/clean

CMakeFiles/systemLib.dir/depend:
	cd /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/src /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build /afs/crc.nd.edu/user/k/ktsai/BuddingCode_nucleus5/budding_cmake/build/CMakeFiles/systemLib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/systemLib.dir/depend

