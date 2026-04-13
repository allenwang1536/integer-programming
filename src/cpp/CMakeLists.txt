cmake_minimum_required(VERSION 3.20)
project(ortools_shim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(ortools CONFIG REQUIRED)

add_library(ortools_shim SHARED ortools_shim.cc)
target_include_directories(ortools_shim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ortools_shim PRIVATE ortools::ortools)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS ortools_shim LIBRARY DESTINATION lib)
