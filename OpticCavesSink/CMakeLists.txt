cmake_minimum_required(VERSION 2.8)
project(unitEvolution)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
set (CMAKE_CXX_STANDARD 11)
set(SRC Lib/Htc.cpp Lib/unitEvolution.cpp main.cpp)
add_executable(opticCaves ${SRC})
target_link_libraries(opticCaves scalapack lapack blas gfortran  ${MPI_C_LIBRARIES})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Data/a.txt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Data/w.txt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Data/phi.txt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
