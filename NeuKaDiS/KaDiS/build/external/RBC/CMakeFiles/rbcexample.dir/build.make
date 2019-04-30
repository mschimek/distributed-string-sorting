# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_SOURCE_DIR = /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build

# Include any dependencies generated for this target.
include external/RBC/CMakeFiles/rbcexample.dir/depend.make

# Include the progress variables for this target.
include external/RBC/CMakeFiles/rbcexample.dir/progress.make

# Include the compile flags for this target's objects.
include external/RBC/CMakeFiles/rbcexample.dir/flags.make

external/RBC/CMakeFiles/rbcexample.dir/example/rbc_example.cpp.o: external/RBC/CMakeFiles/rbcexample.dir/flags.make
external/RBC/CMakeFiles/rbcexample.dir/example/rbc_example.cpp.o: ../external/RBC/example/rbc_example.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/RBC/CMakeFiles/rbcexample.dir/example/rbc_example.cpp.o"
	cd /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC && /usr/lib64/openmpi/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rbcexample.dir/example/rbc_example.cpp.o -c /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC/example/rbc_example.cpp

external/RBC/CMakeFiles/rbcexample.dir/example/rbc_example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rbcexample.dir/example/rbc_example.cpp.i"
	cd /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC && /usr/lib64/openmpi/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC/example/rbc_example.cpp > CMakeFiles/rbcexample.dir/example/rbc_example.cpp.i

external/RBC/CMakeFiles/rbcexample.dir/example/rbc_example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rbcexample.dir/example/rbc_example.cpp.s"
	cd /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC && /usr/lib64/openmpi/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC/example/rbc_example.cpp -o CMakeFiles/rbcexample.dir/example/rbc_example.cpp.s

# Object files for target rbcexample
rbcexample_OBJECTS = \
"CMakeFiles/rbcexample.dir/example/rbc_example.cpp.o"

# External object files for target rbcexample
rbcexample_EXTERNAL_OBJECTS =

external/RBC/rbcexample: external/RBC/CMakeFiles/rbcexample.dir/example/rbc_example.cpp.o
external/RBC/rbcexample: external/RBC/CMakeFiles/rbcexample.dir/build.make
external/RBC/rbcexample: external/RBC/librbc.a
external/RBC/rbcexample: external/RBC/external/tlx/tlx/libtlx.a
external/RBC/rbcexample: external/RBC/CMakeFiles/rbcexample.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable rbcexample"
	cd /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rbcexample.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/RBC/CMakeFiles/rbcexample.dir/build: external/RBC/rbcexample

.PHONY : external/RBC/CMakeFiles/rbcexample.dir/build

external/RBC/CMakeFiles/rbcexample.dir/clean:
	cd /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC && $(CMAKE_COMMAND) -P CMakeFiles/rbcexample.dir/cmake_clean.cmake
.PHONY : external/RBC/CMakeFiles/rbcexample.dir/clean

external/RBC/CMakeFiles/rbcexample.dir/depend:
	cd /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/external/RBC /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC /home/matthias/Dokumente/Studium_Karlsruhe/Semester11/Masterarbeit/DistributedStringSorting/distributed-string-sorting/KaDiS/build/external/RBC/CMakeFiles/rbcexample.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/RBC/CMakeFiles/rbcexample.dir/depend
