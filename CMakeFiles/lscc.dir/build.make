# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/user/geant4/lscc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/user/geant4/lscc

# Include any dependencies generated for this target.
include CMakeFiles/lscc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lscc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lscc.dir/flags.make

CMakeFiles/lscc.dir/lscc.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/lscc.cc.o: lscc.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lscc.dir/lscc.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/lscc.cc.o -c /home/user/geant4/lscc/lscc.cc

CMakeFiles/lscc.dir/lscc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/lscc.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/lscc.cc > CMakeFiles/lscc.dir/lscc.cc.i

CMakeFiles/lscc.dir/lscc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/lscc.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/lscc.cc -o CMakeFiles/lscc.dir/lscc.cc.s

CMakeFiles/lscc.dir/lscc.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/lscc.cc.o.requires

CMakeFiles/lscc.dir/lscc.cc.o.provides: CMakeFiles/lscc.dir/lscc.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/lscc.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/lscc.cc.o.provides

CMakeFiles/lscc.dir/lscc.cc.o.provides.build: CMakeFiles/lscc.dir/lscc.cc.o


CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o: src/PrimaryGeneratorAction_GPS.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o -c /home/user/geant4/lscc/src/PrimaryGeneratorAction_GPS.cc

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/PrimaryGeneratorAction_GPS.cc > CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.i

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/PrimaryGeneratorAction_GPS.cc -o CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.s

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.requires

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.provides: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.provides

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.provides.build: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o


CMakeFiles/lscc.dir/src/EventAction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/EventAction.cc.o: src/EventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/lscc.dir/src/EventAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/EventAction.cc.o -c /home/user/geant4/lscc/src/EventAction.cc

CMakeFiles/lscc.dir/src/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/EventAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/EventAction.cc > CMakeFiles/lscc.dir/src/EventAction.cc.i

CMakeFiles/lscc.dir/src/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/EventAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/EventAction.cc -o CMakeFiles/lscc.dir/src/EventAction.cc.s

CMakeFiles/lscc.dir/src/EventAction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/EventAction.cc.o.requires

CMakeFiles/lscc.dir/src/EventAction.cc.o.provides: CMakeFiles/lscc.dir/src/EventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/EventAction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/EventAction.cc.o.provides

CMakeFiles/lscc.dir/src/EventAction.cc.o.provides.build: CMakeFiles/lscc.dir/src/EventAction.cc.o


CMakeFiles/lscc.dir/src/TrackingAction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/TrackingAction.cc.o: src/TrackingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/lscc.dir/src/TrackingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/TrackingAction.cc.o -c /home/user/geant4/lscc/src/TrackingAction.cc

CMakeFiles/lscc.dir/src/TrackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/TrackingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/TrackingAction.cc > CMakeFiles/lscc.dir/src/TrackingAction.cc.i

CMakeFiles/lscc.dir/src/TrackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/TrackingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/TrackingAction.cc -o CMakeFiles/lscc.dir/src/TrackingAction.cc.s

CMakeFiles/lscc.dir/src/TrackingAction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/TrackingAction.cc.o.requires

CMakeFiles/lscc.dir/src/TrackingAction.cc.o.provides: CMakeFiles/lscc.dir/src/TrackingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/TrackingAction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/TrackingAction.cc.o.provides

CMakeFiles/lscc.dir/src/TrackingAction.cc.o.provides.build: CMakeFiles/lscc.dir/src/TrackingAction.cc.o


CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o: src/DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o -c /home/user/geant4/lscc/src/DetectorConstruction.cc

CMakeFiles/lscc.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/DetectorConstruction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/DetectorConstruction.cc > CMakeFiles/lscc.dir/src/DetectorConstruction.cc.i

CMakeFiles/lscc.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/DetectorConstruction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/DetectorConstruction.cc -o CMakeFiles/lscc.dir/src/DetectorConstruction.cc.s

CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.requires

CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.provides: CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.provides

CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.provides.build: CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o


CMakeFiles/lscc.dir/src/SteppingAction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/SteppingAction.cc.o: src/SteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/lscc.dir/src/SteppingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/SteppingAction.cc.o -c /home/user/geant4/lscc/src/SteppingAction.cc

CMakeFiles/lscc.dir/src/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/SteppingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/SteppingAction.cc > CMakeFiles/lscc.dir/src/SteppingAction.cc.i

CMakeFiles/lscc.dir/src/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/SteppingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/SteppingAction.cc -o CMakeFiles/lscc.dir/src/SteppingAction.cc.s

CMakeFiles/lscc.dir/src/SteppingAction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/SteppingAction.cc.o.requires

CMakeFiles/lscc.dir/src/SteppingAction.cc.o.provides: CMakeFiles/lscc.dir/src/SteppingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/SteppingAction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/SteppingAction.cc.o.provides

CMakeFiles/lscc.dir/src/SteppingAction.cc.o.provides.build: CMakeFiles/lscc.dir/src/SteppingAction.cc.o


CMakeFiles/lscc.dir/src/PhysicsList.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/PhysicsList.cc.o: src/PhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/lscc.dir/src/PhysicsList.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/PhysicsList.cc.o -c /home/user/geant4/lscc/src/PhysicsList.cc

CMakeFiles/lscc.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/PhysicsList.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/PhysicsList.cc > CMakeFiles/lscc.dir/src/PhysicsList.cc.i

CMakeFiles/lscc.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/PhysicsList.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/PhysicsList.cc -o CMakeFiles/lscc.dir/src/PhysicsList.cc.s

CMakeFiles/lscc.dir/src/PhysicsList.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/PhysicsList.cc.o.requires

CMakeFiles/lscc.dir/src/PhysicsList.cc.o.provides: CMakeFiles/lscc.dir/src/PhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/PhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/PhysicsList.cc.o.provides

CMakeFiles/lscc.dir/src/PhysicsList.cc.o.provides.build: CMakeFiles/lscc.dir/src/PhysicsList.cc.o


CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o: src/PrimaryGeneratorAction_PG.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o -c /home/user/geant4/lscc/src/PrimaryGeneratorAction_PG.cc

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/PrimaryGeneratorAction_PG.cc > CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.i

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/PrimaryGeneratorAction_PG.cc -o CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.s

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.requires

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.provides: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.provides

CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.provides.build: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o


CMakeFiles/lscc.dir/src/CCRunAction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/CCRunAction.cc.o: src/CCRunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/lscc.dir/src/CCRunAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/CCRunAction.cc.o -c /home/user/geant4/lscc/src/CCRunAction.cc

CMakeFiles/lscc.dir/src/CCRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/CCRunAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/CCRunAction.cc > CMakeFiles/lscc.dir/src/CCRunAction.cc.i

CMakeFiles/lscc.dir/src/CCRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/CCRunAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/CCRunAction.cc -o CMakeFiles/lscc.dir/src/CCRunAction.cc.s

CMakeFiles/lscc.dir/src/CCRunAction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/CCRunAction.cc.o.requires

CMakeFiles/lscc.dir/src/CCRunAction.cc.o.provides: CMakeFiles/lscc.dir/src/CCRunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/CCRunAction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/CCRunAction.cc.o.provides

CMakeFiles/lscc.dir/src/CCRunAction.cc.o.provides.build: CMakeFiles/lscc.dir/src/CCRunAction.cc.o


CMakeFiles/lscc.dir/src/PhotCntSD.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/PhotCntSD.cc.o: src/PhotCntSD.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/lscc.dir/src/PhotCntSD.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/PhotCntSD.cc.o -c /home/user/geant4/lscc/src/PhotCntSD.cc

CMakeFiles/lscc.dir/src/PhotCntSD.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/PhotCntSD.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/PhotCntSD.cc > CMakeFiles/lscc.dir/src/PhotCntSD.cc.i

CMakeFiles/lscc.dir/src/PhotCntSD.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/PhotCntSD.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/PhotCntSD.cc -o CMakeFiles/lscc.dir/src/PhotCntSD.cc.s

CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.requires

CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.provides: CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.provides

CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.provides.build: CMakeFiles/lscc.dir/src/PhotCntSD.cc.o


CMakeFiles/lscc.dir/src/DEPosSD.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/DEPosSD.cc.o: src/DEPosSD.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/lscc.dir/src/DEPosSD.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/DEPosSD.cc.o -c /home/user/geant4/lscc/src/DEPosSD.cc

CMakeFiles/lscc.dir/src/DEPosSD.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/DEPosSD.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/DEPosSD.cc > CMakeFiles/lscc.dir/src/DEPosSD.cc.i

CMakeFiles/lscc.dir/src/DEPosSD.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/DEPosSD.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/DEPosSD.cc -o CMakeFiles/lscc.dir/src/DEPosSD.cc.s

CMakeFiles/lscc.dir/src/DEPosSD.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/DEPosSD.cc.o.requires

CMakeFiles/lscc.dir/src/DEPosSD.cc.o.provides: CMakeFiles/lscc.dir/src/DEPosSD.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/DEPosSD.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/DEPosSD.cc.o.provides

CMakeFiles/lscc.dir/src/DEPosSD.cc.o.provides.build: CMakeFiles/lscc.dir/src/DEPosSD.cc.o


CMakeFiles/lscc.dir/src/DEPosHit.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/DEPosHit.cc.o: src/DEPosHit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/lscc.dir/src/DEPosHit.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/DEPosHit.cc.o -c /home/user/geant4/lscc/src/DEPosHit.cc

CMakeFiles/lscc.dir/src/DEPosHit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/DEPosHit.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/DEPosHit.cc > CMakeFiles/lscc.dir/src/DEPosHit.cc.i

CMakeFiles/lscc.dir/src/DEPosHit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/DEPosHit.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/DEPosHit.cc -o CMakeFiles/lscc.dir/src/DEPosHit.cc.s

CMakeFiles/lscc.dir/src/DEPosHit.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/DEPosHit.cc.o.requires

CMakeFiles/lscc.dir/src/DEPosHit.cc.o.provides: CMakeFiles/lscc.dir/src/DEPosHit.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/DEPosHit.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/DEPosHit.cc.o.provides

CMakeFiles/lscc.dir/src/DEPosHit.cc.o.provides.build: CMakeFiles/lscc.dir/src/DEPosHit.cc.o


CMakeFiles/lscc.dir/src/RunAction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/RunAction.cc.o: src/RunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/lscc.dir/src/RunAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/RunAction.cc.o -c /home/user/geant4/lscc/src/RunAction.cc

CMakeFiles/lscc.dir/src/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/RunAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/RunAction.cc > CMakeFiles/lscc.dir/src/RunAction.cc.i

CMakeFiles/lscc.dir/src/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/RunAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/RunAction.cc -o CMakeFiles/lscc.dir/src/RunAction.cc.s

CMakeFiles/lscc.dir/src/RunAction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/RunAction.cc.o.requires

CMakeFiles/lscc.dir/src/RunAction.cc.o.provides: CMakeFiles/lscc.dir/src/RunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/RunAction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/RunAction.cc.o.provides

CMakeFiles/lscc.dir/src/RunAction.cc.o.provides.build: CMakeFiles/lscc.dir/src/RunAction.cc.o


CMakeFiles/lscc.dir/src/ActionInitialization.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/ActionInitialization.cc.o: src/ActionInitialization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/lscc.dir/src/ActionInitialization.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/ActionInitialization.cc.o -c /home/user/geant4/lscc/src/ActionInitialization.cc

CMakeFiles/lscc.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/ActionInitialization.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/ActionInitialization.cc > CMakeFiles/lscc.dir/src/ActionInitialization.cc.i

CMakeFiles/lscc.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/ActionInitialization.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/ActionInitialization.cc -o CMakeFiles/lscc.dir/src/ActionInitialization.cc.s

CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.requires

CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.provides: CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.provides

CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.provides.build: CMakeFiles/lscc.dir/src/ActionInitialization.cc.o


CMakeFiles/lscc.dir/src/PhotCntHit.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/PhotCntHit.cc.o: src/PhotCntHit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/lscc.dir/src/PhotCntHit.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/PhotCntHit.cc.o -c /home/user/geant4/lscc/src/PhotCntHit.cc

CMakeFiles/lscc.dir/src/PhotCntHit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/PhotCntHit.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/PhotCntHit.cc > CMakeFiles/lscc.dir/src/PhotCntHit.cc.i

CMakeFiles/lscc.dir/src/PhotCntHit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/PhotCntHit.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/PhotCntHit.cc -o CMakeFiles/lscc.dir/src/PhotCntHit.cc.s

CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.requires

CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.provides: CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.provides

CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.provides.build: CMakeFiles/lscc.dir/src/PhotCntHit.cc.o


CMakeFiles/lscc.dir/src/CCEventAction.cc.o: CMakeFiles/lscc.dir/flags.make
CMakeFiles/lscc.dir/src/CCEventAction.cc.o: src/CCEventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object CMakeFiles/lscc.dir/src/CCEventAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lscc.dir/src/CCEventAction.cc.o -c /home/user/geant4/lscc/src/CCEventAction.cc

CMakeFiles/lscc.dir/src/CCEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lscc.dir/src/CCEventAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/geant4/lscc/src/CCEventAction.cc > CMakeFiles/lscc.dir/src/CCEventAction.cc.i

CMakeFiles/lscc.dir/src/CCEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lscc.dir/src/CCEventAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/geant4/lscc/src/CCEventAction.cc -o CMakeFiles/lscc.dir/src/CCEventAction.cc.s

CMakeFiles/lscc.dir/src/CCEventAction.cc.o.requires:

.PHONY : CMakeFiles/lscc.dir/src/CCEventAction.cc.o.requires

CMakeFiles/lscc.dir/src/CCEventAction.cc.o.provides: CMakeFiles/lscc.dir/src/CCEventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/lscc.dir/build.make CMakeFiles/lscc.dir/src/CCEventAction.cc.o.provides.build
.PHONY : CMakeFiles/lscc.dir/src/CCEventAction.cc.o.provides

CMakeFiles/lscc.dir/src/CCEventAction.cc.o.provides.build: CMakeFiles/lscc.dir/src/CCEventAction.cc.o


# Object files for target lscc
lscc_OBJECTS = \
"CMakeFiles/lscc.dir/lscc.cc.o" \
"CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o" \
"CMakeFiles/lscc.dir/src/EventAction.cc.o" \
"CMakeFiles/lscc.dir/src/TrackingAction.cc.o" \
"CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/lscc.dir/src/SteppingAction.cc.o" \
"CMakeFiles/lscc.dir/src/PhysicsList.cc.o" \
"CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o" \
"CMakeFiles/lscc.dir/src/CCRunAction.cc.o" \
"CMakeFiles/lscc.dir/src/PhotCntSD.cc.o" \
"CMakeFiles/lscc.dir/src/DEPosSD.cc.o" \
"CMakeFiles/lscc.dir/src/DEPosHit.cc.o" \
"CMakeFiles/lscc.dir/src/RunAction.cc.o" \
"CMakeFiles/lscc.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/lscc.dir/src/PhotCntHit.cc.o" \
"CMakeFiles/lscc.dir/src/CCEventAction.cc.o"

# External object files for target lscc
lscc_EXTERNAL_OBJECTS =

lscc: CMakeFiles/lscc.dir/lscc.cc.o
lscc: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o
lscc: CMakeFiles/lscc.dir/src/EventAction.cc.o
lscc: CMakeFiles/lscc.dir/src/TrackingAction.cc.o
lscc: CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o
lscc: CMakeFiles/lscc.dir/src/SteppingAction.cc.o
lscc: CMakeFiles/lscc.dir/src/PhysicsList.cc.o
lscc: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o
lscc: CMakeFiles/lscc.dir/src/CCRunAction.cc.o
lscc: CMakeFiles/lscc.dir/src/PhotCntSD.cc.o
lscc: CMakeFiles/lscc.dir/src/DEPosSD.cc.o
lscc: CMakeFiles/lscc.dir/src/DEPosHit.cc.o
lscc: CMakeFiles/lscc.dir/src/RunAction.cc.o
lscc: CMakeFiles/lscc.dir/src/ActionInitialization.cc.o
lscc: CMakeFiles/lscc.dir/src/PhotCntHit.cc.o
lscc: CMakeFiles/lscc.dir/src/CCEventAction.cc.o
lscc: CMakeFiles/lscc.dir/build.make
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4Tree.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4GMocren.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4visHepRep.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4RayTracer.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4VRML.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4OpenGL.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4gl2ps.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4interfaces.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4persistency.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4analysis.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4error_propagation.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4readout.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4physicslists.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4parmodels.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4FR.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4vis_management.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4modeling.so
lscc: /usr/lib/x86_64-linux-gnu/libXm.so
lscc: /usr/lib/x86_64-linux-gnu/libSM.so
lscc: /usr/lib/x86_64-linux-gnu/libICE.so
lscc: /usr/lib/x86_64-linux-gnu/libX11.so
lscc: /usr/lib/x86_64-linux-gnu/libXext.so
lscc: /usr/lib/x86_64-linux-gnu/libXt.so
lscc: /usr/lib/x86_64-linux-gnu/libXmu.so
lscc: /usr/lib/x86_64-linux-gnu/libGLU.so
lscc: /usr/lib/x86_64-linux-gnu/libGL.so
lscc: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
lscc: /usr/lib/x86_64-linux-gnu/libQtGui.so
lscc: /usr/lib/x86_64-linux-gnu/libQtCore.so
lscc: /usr/lib/x86_64-linux-gnu/libxerces-c.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4run.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4event.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4tracking.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4processes.so
lscc: /usr/lib/x86_64-linux-gnu/libz.so
lscc: /usr/lib/x86_64-linux-gnu/libexpat.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4digits_hits.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4track.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4particles.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4geometry.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4materials.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4graphics_reps.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4intercoms.so
lscc: /usr/local/src/geant4/geant4_1003/lib/libG4global.so
lscc: /usr/local/lib/libCLHEP.so
lscc: CMakeFiles/lscc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/user/geant4/lscc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Linking CXX executable lscc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lscc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lscc.dir/build: lscc

.PHONY : CMakeFiles/lscc.dir/build

CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/lscc.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_GPS.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/EventAction.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/TrackingAction.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/DetectorConstruction.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/SteppingAction.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/PhysicsList.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/PrimaryGeneratorAction_PG.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/CCRunAction.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/PhotCntSD.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/DEPosSD.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/DEPosHit.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/RunAction.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/ActionInitialization.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/PhotCntHit.cc.o.requires
CMakeFiles/lscc.dir/requires: CMakeFiles/lscc.dir/src/CCEventAction.cc.o.requires

.PHONY : CMakeFiles/lscc.dir/requires

CMakeFiles/lscc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lscc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lscc.dir/clean

CMakeFiles/lscc.dir/depend:
	cd /home/user/geant4/lscc && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/user/geant4/lscc /home/user/geant4/lscc /home/user/geant4/lscc /home/user/geant4/lscc /home/user/geant4/lscc/CMakeFiles/lscc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lscc.dir/depend
