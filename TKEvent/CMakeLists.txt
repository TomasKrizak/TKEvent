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
cmake_minimum_required(VERSION 3.3)
project(TKEvent VERSION "2.0.0")
set(CMAKE_CXX_STANDARD 14)

#Finds ROOT package and defines cmake variables
find_package(ROOT REQUIRED)
include_directories(include ${ROOT_INCLUDE_DIRS})

file(GLOB SOURCES "./src/*.cpp" "./src/dicts/*.cpp")
add_library(TKEvent SHARED ${SOURCES})
target_link_libraries(TKEvent PUBLIC ${ROOT_LIBRARIES} ROOT::Geom)
