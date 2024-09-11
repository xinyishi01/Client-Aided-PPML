cmake_minimum_required (VERSION 2.8.11)
project (emp-sh2pc)
set(NAME "emp-sh2pc")
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/semihonest)
include_directories("/usr/include/eigen3")

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} cmake/)

FIND_PACKAGE(emp-ot)
INCLUDE_DIRECTORIES(${EMP-OT_INCLUDE_DIRS})

FIND_PACKAGE(emp-sh2pc)
INCLUDE_DIRECTORIES(${EMP-SH2PC_INCLUDE_DIRS})

FIND_PACKAGE(OpenSSL)
INCLUDE_DIRECTORIES(${OPENSSL_INCLUDE_DIR})

set (CMAKE_C_FLAGS "-pthread -Wall -march=native -O3 -maes")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -ggdb")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++0x")

macro (add_test _name)
	add_executable("test_${_name}" "test/${_name}.cpp" ${basics})
	target_link_libraries("test_${_name}"  relic  ${OPENSSL_LIBRARIES} gmp ${EMP-TOOL_LIBRARIES}) 
endmacro()

macro (add_bench _name)
	add_executable("bench_${_name}" "bench/${_name}.cpp" ${basics})
	target_link_libraries("bench_${_name}"  relic  ${OPENSSL_LIBRARIES} gmp ${EMP-TOOL_LIBRARIES}) 
endmacro()

macro (add_example _name)
	add_executable("example_${_name}" "example/${_name}.cpp" ${basics})
	target_link_libraries("example_${_name}"  relic  ${OPENSSL_LIBRARIES} gmp ${EMP-TOOL_LIBRARIES}) 
endmacro()


set(basics
)

add_executable(server main-server.cpp utils.hpp Layer.hpp Layer1_Server.hpp Layer_Server.hpp NN_Server.hpp)
target_link_libraries(server relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})


add_executable(client main-client.cpp utils.hpp Layerc.hpp Layer1_Client.hpp Layer_Client.hpp NN_Client.hpp)
target_link_libraries(client relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})
