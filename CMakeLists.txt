#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(qshields)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# This allows CMake to find your module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Find ROOT
find_package(ROOT)
# Add ROOTs header paths
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include($ENV{G4INSTALL}/cmake/Modules/FindROOT.cmake)

#Check here for qshieldsENV_GPS_USE
#Use GPS if the environment variable is set, and the user has *not* supplied the CMake 
#command line argument. Otherwise, the argument of the command line has the priority.
#
if("$ENV{qshieldsENV_GPS_USE}" AND NOT DEFINED qshieldsENV_GPS_USE)
  set(qshieldsENV_GPS_USE_DEFAULT ON)
else()
  set(qshieldsENV_GPS_USE_DEFAULT ${qshieldsENV_GPS_USE})
endif()

option(qshieldsENV_GPS_USE "Build qshields with support for GPS" ${qshieldsENV_GPS_USE_DEFAULT})
mark_as_advanced(qshieldsENV_GPS_USE)

if(qshieldsENV_GPS_USE)
  add_definitions(-DqshieldsENV_GPS_USE)
endif()

if($ENV{G4IS10})
  add_definitions(-DG4GT10)
endif()

if($ENV{G4NEWVIS})
  add_definitions(-DG4NEWVIS)
endif()

add_definitions(-O0 -g)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include	$ENV{QSRC}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
if($ENV{G4IS10})
   file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc        $ENV{QSRC}/src/G*.cc)
else()
   file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc        $ENV{QSRC}/src/*.cc)
endif()
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh	$ENV{QSRC}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(qshields qshields.cc ${sources} ${headers})
target_link_libraries(qshields ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build qshields. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(qshields_SCRIPTS)

foreach(_script ${qshields_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
#add_custom_target(qshields DEPENDS qshields)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS qshields DESTINATION bin)

