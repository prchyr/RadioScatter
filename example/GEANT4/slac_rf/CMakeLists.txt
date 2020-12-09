# $Id: CMakeLists.txt 100923 2016-11-03 10:50:34Z gcosmo $

#



#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(slac)



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



unset(CMAKE_MODULE_PATH)



#use the ROOT things
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
include(${ROOT_USE_FILE})
#add_definitions(${ROOT_CXX_FLAGS})

if(DEFINED ENV{RS_SOURCE_DIR})
  set(RS_INCLUDE_DIR $ENV{RS_SOURCE_DIR})
else()
  set(RS_INCLUDE_DIR .)
endif()

if(DEFINED ENV{RS_INSTALL_DIR})
  set(RS_INSTALL_DIR $ENV{RS_INSTALL_DIR}/include)
  set(LIB_INSTALL_DIR $ENV{RS_INSTALL_DIR}/lib)
else()
  set(RS_INSTALL_DIR /usr/local/include/RadioScatter)
  set(LIB_INSTALL_DIR /usr/local/lib)
endif()
# # Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#



include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include  ${ROOT_INCLUDE_DIRS} ${RS_INCLUDE_DIR} ${RS_INCLUDE_DIR}/../)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${RS_INCLUDE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${RS_INCLUDE_DIR}/*.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

add_executable(slac slac.cc ${sources} ${headers})
#add_executable(threshscan threshscan.cc ${sources} ${headers})

target_link_libraries(slac  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} RadioScatter)
#target_link_libraries(threshscan ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} RadioScatterEvent)
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SLAC_SCRIPTS
  slac.out
  slac.in
  gui.mac
  init_vis.mac
  run1.mac
  vis.mac
  )

foreach(_script ${SLAC_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS slac DESTINATION ${PROJECT_BINARY_DIR}/bin)

