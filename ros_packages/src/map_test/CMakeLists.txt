cmake_minimum_required(VERSION 3.8)
project(map_test)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_subdirectory(../../../safe_planner
                 ${CMAKE_BINARY_DIR}/safe_planner_build)

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(ament_cmake_auto REQUIRED)
find_package(PCL REQUIRED COMPONENTS common)

ament_auto_find_build_dependencies()
ament_auto_add_executable(map_test src/map_test.cpp)

target_include_directories(map_test PRIVATE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> 
  $<INSTALL_INTERFACE:include>)
target_link_libraries(map_test safe_planner)

ament_auto_package(USE_SCOPED_HEADER_INSTALL_DIR)
