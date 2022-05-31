cmake_minimum_required(VERSION 3.20)

project(gravity0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_BUILD_TYPE Release)

find_package(HDF5   REQUIRED COMPONENTS CXX HL)
find_package(OpenMP REQUIRED)
find_package(VTK    REQUIRED COMPONENTS IOXML IOImage)
find_package(PkgConfig REQUIRED)
# https://github.com/FFTW/fftw3/issues/130
pkg_check_modules(FFTW IMPORTED_TARGET REQUIRED fftw3)

add_executable(direct_verlet   src/solvers/direct_verlet.cc)
add_executable(direct_leapfrog src/solvers/direct_leapfrog.cc)
add_executable(field_periodic  src/solvers/field_periodic.cc)

target_include_directories(direct_verlet   PRIVATE src/solvers src/storage ${HDF5_CXX_INCLUDE_DIRS} ${fmt_INCLUDE_DIRS})
target_include_directories(direct_leapfrog PRIVATE src/solvers src/storage ${HDF5_CXX_INCLUDE_DIRS} ${fmt_INCLUDE_DIRS})
target_include_directories(field_periodic  PRIVATE src/solvers src/storage ${HDF5_CXX_INCLUDE_DIRS} ${fmt_INCLUDE_DIRS} PkgConfig::FFTW)

target_link_libraries(direct_verlet   PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${VTK_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(direct_leapfrog PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${VTK_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(field_periodic  PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${VTK_LIBRARIES} PkgConfig::FFTW)
