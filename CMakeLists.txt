cmake_minimum_required(VERSION 3.20)
project(gravity0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-math-errno -ffinite-math-only")
find_package(HDF5 COMPONENTS CXX HL)
find_package(OpenMP)
find_package(VTK COMPONENTS IOXdmf3)
add_executable(direct_verlet src/solvers/direct_verlet.cc)
add_executable(direct_leapfrog src/solvers/direct_leapfrog.cc)
include_directories(
    ${HDF5_CXX_INCLUDE_DIRS}
)
target_link_libraries(
    direct_verlet
    OpenMP::OpenMP_CXX
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${VTK_LIBRARIES}
)
target_link_libraries(
    direct_leapfrog
    OpenMP::OpenMP_CXX
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${VTK_LIBRARIES}
)
