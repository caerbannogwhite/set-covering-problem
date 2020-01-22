cmake_minimum_required(VERSION 3.14)

project(set_covering C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER /usr/bin/gcc)

add_executable(set_covering
        src/main.c src/main.h
        src/sc.c src/sc.h
        src/aux.c src/aux.h
        src/balas_dense.c src/balas_dense.h
        src/callbacks.c src/callbacks.h
        src/preprocessing.c src/preprocessing.h src/balas_sparse.c src/balas_sparse.h)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /home/macs/coding/optimization/set_covering)
find_package(CPLEX)

#message(">>> CPLEX FOUND: ${CPLEX_FOUND}")
#message(">>> CPLEX_INCLUDE_DIR: ${CPLEX_INCLUDE_DIR}")

target_include_directories(set_covering PRIVATE ${CPLEX_INCLUDE_DIR}/ilcplex)
target_include_directories(set_covering PRIVATE /usr/local/Cellar/openblas/0.3.7/include)

target_link_libraries(set_covering cplex-library)
target_link_libraries(set_covering m)
target_link_libraries(set_covering dl)
target_link_libraries(set_covering pthread)
target_link_libraries(set_covering gfortran)
target_link_libraries(set_covering /opt/OpenBLAS/lib/libopenblas.a)

add_compile_options(O2)
