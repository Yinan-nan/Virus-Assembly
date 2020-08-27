cmake_minimum_required(VERSION 3.7)
project(shellgrow)
#add_subdirectory(source)

set(CMAKE_CXX_STANDARD 11)

#file(GLOB SOURCES
#    *.h
#    *.cpp
#)

set(SOURCE_FILES source/main.cpp source/shell.cpp source/shell.h source/parameters.h source/energy.cpp source/allvectors.cpp source/Arrow.cpp source/Arrow.h source/grow.cpp source/Tensor.cpp source/Tensor.h)


add_executable(shellgrow ${SOURCE_FILES} ${SOURCES})
target_link_libraries(shellgrow gsl gslcblas m)

#add_custom_target(run
#   COMMAND binary
#   DEPENDS binary
#   WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
#)