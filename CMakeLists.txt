cmake_minimum_required(VERSION 3.20)
project(ThermalPlate LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(include)

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

add_executable(plate src/main.cpp src/func.cpp)

target_include_directories(plate PRIVATE
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS})

target_link_libraries(plate PRIVATE ${Python3_LIBRARIES})

target_compile_definitions(plate PRIVATE MATPLOTLIBCPP_PYTHON_HEADER_ONLY)
