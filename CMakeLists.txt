cmake_minimum_required(VERSION 3.22)
project(N_body_simulation)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

aux_source_directory(include INC_FILES)

add_executable(NaiveMethod src/NaiveMethod.cpp ${INC_FILES})
add_executable(BarnesHut src/BarnesHut.cpp ${INC_FILES})
add_executable(FMM src/FMM.cpp ${INC_FILES})
add_executable(DataGenerator src/DataGenerator.cpp ${INC_FILES})