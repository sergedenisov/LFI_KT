cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME demo)
project(${PROJECT_NAME})

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)
find_package(DOLFIN REQUIRED COMPONENTS ufc)
include(${DOLFIN_USE_FILE})

add_executable(poisson simulate.cpp)
target_link_libraries(poisson dolfin)


