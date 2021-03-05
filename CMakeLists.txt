# $Id: CMakeLists.txt 100923 2016-11-03 10:50:34Z gcosmo $

#



#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(rs)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
unset(CMAKE_MODULE_PATH)
#use the CLHEP things
#will generate a bunch of library modules CLHEP::stuff that can be linked to
#find_package(CLHEP REQUIRED)
find_package(GSL REQUIRED)

#use the ROOT things
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#list(APPEND CMAKE_PREFIX_PATH $ENV{CLHEP_ROOT_DIR})
find_package(ROOT REQUIRED COMPONENTS RIO Tree MathMore)
include(${ROOT_USE_FILE})
#add_definitions(${ROOT_CXX_FLAGS})

if(DEFINED ENV{RS_SOURCE_DIR})
  set(RS_INCLUDE_DIR $ENV{RS_SOURCE_DIR}/include)
  set(RS_SOURCE_DIR $ENV{RS_SOURCE_DIR}/src)
else()
  set(RS_INCLUDE_DIR ./include)
  set(RS_INCLUDE_DIR ./src)
endif()

if(DEFINED ENV{RS_INSTALL_DIR})
  set(RS_INSTALL_DIR $ENV{RS_INSTALL_DIR}/include/RadioScatter)
  set(LIB_INSTALL_DIR $ENV{RS_INSTALL_DIR}/lib)
else()
  set(RS_INSTALL_DIR /usr/local/include/RadioScatter)
  set(LIB_INSTALL_DIR /usr/local/lib)
endif()


include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${RS_INCLUDE_DIR} ${RS_INCLUDE_DIR}/../ ${GSL_INCLUDE_DIRS} )

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${RS_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${RS_INCLUDE_DIR}/*.hh)

file(GLOB SOURCE_FILES ${RS_SOURCE_DIR}/*.cc )
file(GLOB HEADER_FILES  ${RS_INCLUDE_DIR}/*.hh ${GSL_INCLUDE_DIRS}/gsl/*.h)


#make the ROOT dictionaries and shared object library for reading the event tree

ROOT_GENERATE_DICTIONARY(G__RadioScatter ${HEADER_FILES} LINKDEF ${PROJECT_SOURCE_DIR}/LinkDef.h)


add_library(RadioScatter SHARED ${SOURCE_FILES} G__RadioScatter.cxx)



target_link_libraries(RadioScatter ${ROOT_LIBRARIES} GSL::gsl GSL::gslcblas)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS slac DESTINATION ${PROJECT_BINARY_DIR}/bin)
install (FILES ${RS_INCLUDE_DIR}/RadioScatter.hh DESTINATION ${RS_INSTALL_DIR})
#install (FILES ${RS_SOURCE_DIR}/RadioScatter.cc DESTINATION ${RS_INSTALL_DIR})
install (FILES ${RS_INCLUDE_DIR}/RadioScatterEvent.hh DESTINATION ${RS_INSTALL_DIR})
#install (FILES ${RS_SOURCE_DIR}/RadioScatterEvent.cc DESTINATION ${RS_INSTALL_DIR})
install (FILES ${PROJECT_BINARY_DIR}/libRadioScatter_rdict.pcm DESTINATION ${LIB_INSTALL_DIR})

install (FILES ${RS_INCLUDE_DIR}/RSEventSummary.hh DESTINATION ${RS_INSTALL_DIR})
#install (FILES ${RS_SOURCE_DIR}/RSEventSummary.cc DESTINATION ${RS_INSTALL_DIR})
#install (FILES ${PROJECT_BINARY_DIR}/libRSEventSummary_rdict.pcm DESTINATION ${LIB_INSTALL_DIR})
install (FILES ${RS_INCLUDE_DIR}/TUtilRadioScatter.hh DESTINATION ${RS_INSTALL_DIR})
install (FILES ${RS_INCLUDE_DIR}/IceRayTracing.hh DESTINATION ${RS_INSTALL_DIR})
install (FILES ${RS_SOURCE_DIR}/IceRayTracing.cc DESTINATION ${RS_INSTALL_DIR})

install (TARGETS RadioScatter LIBRARY DESTINATION ${LIB_INSTALL_DIR})
#install (TARGETS RSEventSummary LIBRARY DESTINATION ${LIB_INSTALL_DIR})
