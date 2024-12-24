cmake_minimum_required(VERSION 3.10)
project(fluid)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(DEFINED TYPES)
	add_compile_definitions("-DTYPES=${TYPES}")
	if(DEFINED SIZES)
	  add_compile_definitions("-DSIZES=${SIZES}")
	endif(DEFINED SIZES)
	add_executable(main main.cpp)
else(DEFINED TYPES)
	message(SEND_ERROR "TYPES UNDEFINED")
endif(DEFINED TYPES)
