# Install script for directory: /home/bt308081/build/eigen-3.4.0/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/AdolcForward"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/AlignedVector3"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/ArpackSupport"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/AutoDiff"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/BVH"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/EulerAngles"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/FFT"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/IterativeSolvers"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/KroneckerProduct"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/LevenbergMarquardt"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/MatrixFunctions"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/MoreVectorization"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/MPRealSupport"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/NonLinearOptimization"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/NumericalDiff"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/OpenGLSupport"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/Polynomials"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/Skyline"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/SparseExtra"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/SpecialFunctions"
    "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/home/bt308081/build/eigen-3.4.0/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/bt308081/build/eigen-3.4.0/build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

