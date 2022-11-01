
message("CMAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%")
# Require out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(STATUS "The project requires an out-of-source build.")
  message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
  message(STATUS "  CMakeCache.txt")
  message(STATUS "  CMakeFiles")
  message(STATUS "Once these files are removed, create a separate directory")
  message(STATUS "and run CMake from there, pointing it to:")
  message(STATUS "  ${CMAKE_SOURCE_DIR}")
  message(FATAL_ERROR "in-source build detected")
endif()

# Project setup
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(snfee_examples_read_red VERSION 0.1.0)

# - Load Builtin/Custom Modules
include(GNUInstallDirs)
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)
message(STATUS "[info] CMAKE_MODULE_PATH = '${CMAKE_MODULE_PATH}'")

# - Dependencies
find_package(SNFrontEndElectronics REQUIRED)
include_directories(${SNFrontEndElectronics_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${Event3D_PATH}/include/)

# - Executable:
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
file(GLOB_RECURSE SOURCES2 . $ENV{ROOTSYS}/lib/*.so ${Event3D_PATH}/lib/*.so) # 
foreach(source ${SOURCES})
    file(RELATIVE_PATH cxx_filename ${CMAKE_CURRENT_SOURCE_DIR} ${source})
    string(REPLACE ".cxx" "" exe_filename ${cxx_filename})
    message(STATUS "adding executable ${exe_filename}")
    add_executable(${exe_filename} ${cxx_filename} )
    target_link_libraries(${exe_filename} SNFrontEndElectronics::snfee ${SOURCES2})
  #  _library(${exe_filename} SHARED  )  
endforeach(source ${SOURCES})
