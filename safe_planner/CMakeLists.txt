cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


project(safest_planner)


file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/**.cpp)
add_library(safe_planner ${srcs})

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED COMPONENTS common)


add_definitions(${PCL_DEFINITIONS})

target_include_directories(safe_planner PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc> 
    $<INSTALL_INTERFACE:include>
    ${PCL_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS}
)

target_link_libraries(safe_planner PRIVATE 
    ${PCL_LIBRARIES}
)
install(TARGETS safe_planner EXPORT safe_plannerTargets DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY inc/ DESTINATION include/safe_planner)
