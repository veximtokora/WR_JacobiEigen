cmake_minimum_required(VERSION 3.17)
project(WR_Jacobi_Eigen)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3 -pthread")

add_executable(WR_Jacobi_Eigen src/main.cpp 
                               src/utils.h
			       src/testMatrices.h)
