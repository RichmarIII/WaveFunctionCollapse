# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

# Project Name and Version
project(WaveFunctionCollapse VERSION 1.0 LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable optimizations for release builds
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add the executable
add_executable(wfc main.cpp WFC.hpp stb_image.h stb_image_write.h)

# Include directories (current directory)
target_include_directories(wfc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Compiler Flags for Optimization
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(wfc PRIVATE -O3)
endif()

# Optionally, enable OpenMP for parallel processing
# Uncomment the following lines to enable OpenMP
# find_package(OpenMP REQUIRED)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(wfc PRIVATE OpenMP::OpenMP_CXX)
#     target_compile_options(wfc PRIVATE ${OpenMP_CXX_FLAGS})
# endif()
