cmake_minimum_required(VERSION 3.5)

project(inferPhysicsProject)

# Some flags for Freeglut and GLUI.
#SET(CMAKE_CXX_COMPILER clang++)
#ADD_DEFINITIONS("-std=c++14 -g -O0 -Wall -Wuninitialized")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")


#find_package(OpenGL REQUIRED)
#find_package(SDL2 REQUIRED)
#find_package(EIGEN3 REQUIRED)

#find_file(Data SpringForce.json)
#find_file(Data ProjectileMotion.json)
#find_file(Data actualPath)
#find_file(Data predictedPath)



SET(EIGEN_INCLUDE_DIR EIGEN_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to eigen")


message("${CMAKE_CXX_COMPILER_ID}")
IF(APPLE)
  ADD_DEFINITIONS("-std=c++1z -g -O3 -Wall")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  ADD_DEFINITIONS("-std=c++1z -g -O3 -Wall")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  ADD_DEFINITIONS("-std=c++14 -g -O3 -Wall")
elseif(MSVC)
  ADD_DEFINITIONS("/std:c++14 /Wall")
endif()

include_directories (
  ${EIGEN_INCLUDE_DIR}
  )


set (SOURCES
  main.cpp
  Artist.cpp
  jsoncpp/dist/jsoncpp
  CreatePath.cpp
  )

add_executable(inferPhysics
  ${SOURCES}
  )

target_link_libraries ( 
  inferPhysics
  ) 


