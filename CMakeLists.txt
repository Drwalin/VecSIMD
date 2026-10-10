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
target_compile_options(ray_triangle_test PRIVATE -Wall -Wextra)

add_test(NAME ray_triangle_test COMMAND ray_triangle_test)

# Tests for simd::aabb4 (Aabb.hpp). accept [c1]: builds clean with -Wall -Wextra;
# accept [c2]: runs and exits 0 via ctest.
add_executable(aabb4_test tests/aabb4_test.cpp)
target_link_libraries(aabb4_test PRIVATE VecSIMD)
target_compile_options(aabb4_test PRIVATE -Wall -Wextra)

add_test(NAME aabb4_test COMMAND aabb4_test)
