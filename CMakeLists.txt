cmake_minimum_required(VERSION 3.16)
project(VecSIMD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(VecSIMD INTERFACE)
target_include_directories(VecSIMD INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Tests for the watertight ray/triangle intersection (RayTriangleTest.cpp).
enable_testing()

add_executable(ray_triangle_test tests/ray_triangle_test.cpp)
target_link_libraries(ray_triangle_test PRIVATE VecSIMD)

if(MSVC)
  target_compile_options(ray_triangle_test PRIVATE /W4)
else()
  target_compile_options(ray_triangle_test PRIVATE -Wall -Wextra)
endif()

add_test(NAME ray_triangle_test COMMAND ray_triangle_test)
