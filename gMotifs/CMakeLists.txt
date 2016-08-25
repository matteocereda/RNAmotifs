project(RBM)
cmake_minimum_required(VERSION 2.6)

include_directories(${RBM_SOURCE_DIR})

ADD_EXECUTABLE(tetramer tetramer.cpp anyoption.cpp inc.cpp geco_array.cpp geco_base.cpp)
ADD_EXECUTABLE(counting counting.cpp anyoption.cpp inc.cpp geco_array.cpp geco_base.cpp)

INSTALL(TARGETS tetramer DESTINATION bin)
INSTALL(TARGETS counting DESTINATION bin)
