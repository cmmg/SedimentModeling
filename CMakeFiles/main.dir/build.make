# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns"

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cc.o: main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cc.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cc.o -c "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns/main.cc"

CMakeFiles/main.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cc.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns/main.cc" > CMakeFiles/main.dir/main.cc.i

CMakeFiles/main.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cc.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns/main.cc" -o CMakeFiles/main.dir/main.cc.s

CMakeFiles/main.dir/main.cc.o.requires:

.PHONY : CMakeFiles/main.dir/main.cc.o.requires

CMakeFiles/main.dir/main.cc.o.provides: CMakeFiles/main.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cc.o.provides

CMakeFiles/main.dir/main.cc.o.provides.build: CMakeFiles/main.dir/main.cc.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cc.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/main.cc.o
main: CMakeFiles/main.dir/build.make
main: /home/kbhagat2/software/dealiiInstall/deal.II-v9.0.0/lib/libdeal_II.g.so.9.0.0
main: /home/kbhagat2/software/dealiiInstall/p4est-2.0/DEBUG/lib/libp4est.so
main: /home/kbhagat2/software/dealiiInstall/p4est-2.0/DEBUG/lib/libsc.so
main: /usr/lib/x86_64-linux-gnu/libz.so
main: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
main: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
main: /usr/lib/x86_64-linux-gnu/libboost_system.so
main: /usr/lib/x86_64-linux-gnu/libboost_thread.so
main: /usr/lib/x86_64-linux-gnu/libboost_regex.so
main: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
main: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
main: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libmuelu-adapters.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libmuelu-interface.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libmuelu.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteko.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libstratimikos.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libstratimikosbelos.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libstratimikosaztecoo.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libstratimikosamesos.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libstratimikosml.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libstratimikosifpack.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libifpack2-adapters.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libifpack2.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libanasazitpetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libModeLaplace.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libanasaziepetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libanasazi.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libamesos2.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libbelostpetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libbelosepetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libbelos.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libml.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libifpack.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libzoltan2.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libpamgen_extras.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libpamgen.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libamesos.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libgaleri-xpetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libgaleri-epetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libaztecoo.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libisorropia.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libxpetra-sup.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libxpetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libthyratpetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libthyraepetraext.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libthyraepetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libthyracore.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libepetraext.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetraext.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetrainout.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libkokkostsqr.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetrakernels.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetraclassiclinalg.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetraclassicnodeapi.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtpetraclassic.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libtriutils.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libzoltan.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libepetra.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libsacado.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/librtop.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchoskokkoscomm.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchoskokkoscompat.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchosremainder.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchosnumerics.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchoscomm.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchosparameterlist.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libteuchoscore.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libkokkosalgorithms.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libkokkoscontainers.so
main: /home/kbhagat2/software/dealiiInstall/trilinos-release-12-10-1/lib/libkokkoscore.so
main: /usr/lib/x86_64-linux-gnu/libumfpack.so
main: /usr/lib/x86_64-linux-gnu/libcholmod.so
main: /usr/lib/x86_64-linux-gnu/libccolamd.so
main: /usr/lib/x86_64-linux-gnu/libcolamd.so
main: /usr/lib/x86_64-linux-gnu/libcamd.so
main: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
main: /usr/lib/x86_64-linux-gnu/libamd.so
main: /home/kbhagat2/software/dealiiInstall/hdf5-1.10.1/lib/libhdf5_hl.so
main: /home/kbhagat2/software/dealiiInstall/hdf5-1.10.1/lib/libhdf5.so
main: /home/kbhagat2/software/dealiiInstall/slepc-3.7.3/lib/libslepc.so
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libpetsc.so
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libcmumps.a
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libdmumps.a
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libsmumps.a
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libzmumps.a
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libmumps_common.a
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libpord.a
main: /home/kbhagat2/software/dealiiInstall/parmetis-4.0.3/lib/libparmetis.so
main: /home/kbhagat2/software/dealiiInstall/parmetis-4.0.3/lib/libmetis.so
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libHYPRE.a
main: /home/kbhagat2/software/dealiiInstall/petsc-3.7.6/lib/libscalapack.a
main: /usr/lib/x86_64-linux-gnu/liblapack.so
main: /usr/lib/x86_64-linux-gnu/libblas.so
main: /usr/lib/x86_64-linux-gnu/libhwloc.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempi_ignore_tkr.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cc.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns" "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns" "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns" "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns" "/home/kbhagat2/workspace/Dendritic_solidifcation_2 Eqns/CMakeFiles/main.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

