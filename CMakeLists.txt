cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Wall")

project(try)

include_directories( ${usr/local/include}/src )

set(SOURCE_FILES
	basic_funcs.hpp
	empirical_distr.hpp
	resampling.hpp
	equal_distr_statistics.hpp
	equal_distr_test.hpp
	critical_vector_comparison.cpp)

add_executable(critical_vector_comparison ${SOURCE_FILES})