cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(gravity0)

find_package(HDF5 REQUIRED COMPONENTS CXX HL)

add_executable(gravity0 main.cc)
set_property(TARGET gravity0 PROPERTY CXX_STANDARD 17)

target_include_directories(gravity0 PRIVATE ${HDF5_CXX_INCLUDE_DIRS})
target_link_libraries(gravity0 PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
