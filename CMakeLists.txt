cmake_minimum_required(VERSION 3.5)

project(inferPhysicsProject)

# Some flags for Freeglut and GLUI.
#SET(CMAKE_CXX_COMPILER clang++)
#ADD_DEFINITIONS("-std=c++14 -g -O0 -Wall -Wuninitialized")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")


#find_package(OpenGL REQUIRED)
#find_package(SDL2 REQUIRED)
find_package(EIGEN3 REQUIRED)



include_directories (
  ${EIGEN3_INCLUDE_DIR}
  )

set (SOURCES
  main.cpp
  )

add_executable(inferPhysics
  ${SOURCES}
  )

target_link_libraries ( 
  inferPhysics
  ) 


