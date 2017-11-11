cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)

project(try)

set(SOURCE_FILES
	basic_funcs.hpp
	empirical_distr.hpp
	resampling.hpp
	equal_distr_statistics.hpp
	equal_distr_test.hpp
	main.cpp)

add_executable(equal_test ${SOURCE_FILES})