cmake_minimum_required(VERSION 3.20)

project(healthcare_analytics_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/or-tools" "/opt/homebrew")

find_package(ortools CONFIG REQUIRED)

add_executable(solver
  src/main.cpp
  src/ipinstance.cpp
  src/timer.cpp
)

target_include_directories(solver PRIVATE src)
target_link_libraries(solver PRIVATE ortools::ortools)
