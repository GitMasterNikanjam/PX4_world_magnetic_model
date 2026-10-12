cmake_minimum_required(VERSION 3.16)

project(PX4_world_magnetic_model LANGUAGES CXX)

add_library(world_magnetic_model geo_mag_declination.cpp)
add_library(PX4::world_magnetic_model ALIAS world_magnetic_model)
target_compile_features(world_magnetic_model PUBLIC cxx_std_11)
target_include_directories(world_magnetic_model PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
