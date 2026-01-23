cmake_minimum_required(VERSION 3.10)
project(safest_planner)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(libsafest_planner)
target_include_directories(libsafest_planner PUBLIC inc)