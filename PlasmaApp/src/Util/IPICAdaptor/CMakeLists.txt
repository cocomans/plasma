cmake_minimum_required(VERSION 2.8.8)
project(IPICAdaptor)



find_package(ParaView REQUIRED COMPONENTS vtkPVCatalyst vtkPVPythonCatalyst)
include(${PARAVIEW_USE_FILE})

option(BUILD_SHARED_LIBS "Build ParaView using shared libraries" ON)


include_directories(
  ../../../src
  ../../Util
)

set(IPICAdaptor_SRCS
  IPICAdaptor.cxx
)
add_library(IPICAdaptor ${IPICAdaptor_SRCS})
target_link_libraries(
  IPICAdaptor vtkPVCatalyst vtkPVPythonCatalyst
)
