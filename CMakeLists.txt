#  Cypress -- C++ Spiking Neural Network Simulation Framework
#  Copyright (C) 2019  Christoph Ostrau
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# CMake version and project metadata
cmake_minimum_required(VERSION 3.0)
project(euter_example LANGUAGES CXX)

# Enable C++17 support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# Necessary for pymarocco
find_package(PythonLibs 2.7 REQUIRED )
find_package(Boost REQUIRED serialization system filesystem mpi python27)

# Enable all warnings on MSVC and GCC/Clang/Intel
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
	endif()
elseif(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
		("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
		("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# Setup all includes, we use current github version for includes, which might break at some point -> BS do not deploy headers yet
execute_process(COMMAND "../setup_includes.sh"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set(BS_INCLUDE_PATH ${CMAKE_BINARY_DIR}/bs_includes) 
include_directories(
    PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${Boost_INCLUDE_DIRS} 
        ${PYTHON_INCLUDE_DIR}
        ${BS_INCLUDE_PATH}/euter/include
        ${BS_INCLUDE_PATH}/lib-boost-patches
        ${BS_INCLUDE_PATH}/marocco
        ${BS_INCLUDE_PATH}/pywrap/src
        ${BS_INCLUDE_PATH}/halbe
        ${BS_INCLUDE_PATH}/halco
        ${BS_INCLUDE_PATH}/ztl
        ${BS_INCLUDE_PATH}/rant
        ${BS_INCLUDE_PATH}/redman/include
        ${BS_INCLUDE_PATH}/sthal
        ${BS_INCLUDE_PATH}/pyublas/pyublas/include
        ${BS_INCLUDE_PATH}/calibtic/include
        ${BS_INCLUDE_PATH}/hicann-system/units/communication/source
        ${BS_INCLUDE_PATH}/hwdb
        ${BS_INCLUDE_PATH}/ester
)

# Set in BS Software, if not set, link will not work
add_compile_options(-D BOOST_VARIANT_DO_NOT_USE_VARIADIC_TEMPLATES)

# Manually set library path, because it is not included in LIBRARY_PATH environment variable
link_directories("/wang/environment/software/container/meta-nmpm-software/current/lib")

# Actually compile and link 
add_executable(euter_test euter_test.cpp)
target_link_libraries(euter_test
    ester-submit
    euter
    marocco
    marocco_runtime
    marocco_parameters
    marocco_coordinates
    marocco_results
    pymarocco
    halco_common
    halco_hicann_v2
    sthal
    halbe_container
    log4cxx.so.10
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
)
