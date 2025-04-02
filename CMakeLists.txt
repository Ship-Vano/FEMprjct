cmake_minimum_required(VERSION 3.27)
project(FEMprjct)

set(CMAKE_CXX_STANDARD 23)

include_directories(${CMAKE_SOURCE_DIR}/InputData)
include_directories(${CMAKE_SOURCE_DIR}/OutputData)
include_directories(${CMAKE_SOURCE_DIR}/src/include)
include_directories(${CMAKE_SOURCE_DIR}/src/geometry)

# Find OpenMP
find_package(OpenMP REQUIRED)



add_executable(FEMprjct
        src/fem/FEMSolver2D.cpp
        src/geometry/NetGeometry.cpp
        src/fem/solver.cpp
        src/include/json/jsoncpp.cpp)

#  Geometry executable
add_executable(Geometry
        src/geometry/worldGen.cpp
        src/geometry/NetGeometry.cpp
        src/include/json/jsoncpp.cpp
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Geometry PRIVATE OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()
target_compile_options(Geometry PRIVATE -fopenmp)
target_link_libraries(Geometry PRIVATE -fopenmp)
