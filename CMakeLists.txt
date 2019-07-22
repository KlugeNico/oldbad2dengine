cmake_minimum_required(VERSION 3.7)

project(engine)

#set(CMAKE_CXX_STANDARD 14)

find_package(SDL2 REQUIRED)
find_package(OpenGL REQUIRED)

include_directories(${SDL2_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} )

include_directories(.)

add_executable(engine
        kengine.cpp
        kengine.hpp
        lodepng.cpp
        lodepng.h
        main.cpp)

target_link_libraries(engine ${SDL2_LIBRARIES} ${OPENGL_LIBRARIES} )