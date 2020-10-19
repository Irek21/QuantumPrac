cmake_minimum_required(VERSION 2.8)
project(unitEvolution)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
set (CMAKE_CXX_STANDARD 11)
set(SRC unitEvolution.cpp main.cpp)
add_executable(unitEvolution ${SRC})
target_link_libraries(unitEvolution scalapack lapack blas gfortran  ${MPI_C_LIBRARIES})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/H
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ro
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
