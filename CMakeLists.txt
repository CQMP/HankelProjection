cmake_minimum_required(VERSION 3.19)

project(LegendreHankel C CXX)
set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories("include")
include_directories(${EIGEN3_INCLUDE_DIR})
add_subdirectory(src)

SET(CMAKE_CXX_FLAGS "-O3 -g -fno-omit-frame-pointer")

add_executable(dyproject src/dyproject.cpp src/aux.cpp src/project.cpp)

# Installation commands
install(TARGETS dyproject DESTINATION bin) # Install executable
install(DIRECTORY include/ DESTINATION include) # Install headers