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

#circuit tests
#add_test (bit)
#add_test (int)
#add_test (float)

#add_test (circuit_file)

#add_example (semi_honest)

#install(FILES emp-sh2pc-config.cmake DESTINATION cmake/)

#file(GLOB  headers 
#"semihonest/*.h"
#"emp-sh2pc.h"
#)
#install(FILES ${headers}  DESTINATION include/${NAME})
#add_executable(generate_MT generate_MT.cpp)
#target_link_libraries(generate_MT  relic  ${OPENSSL_LIBRARIES} m ntl gmp ${EMP-TOOL_LIBRARIES})


#add_executable(logistic_nontl logistic_nontl.cpp)
#target_link_libraries(logistic_nontl  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

#add_executable(otmt_nontl otmt_nontl.cpp)
#target_link_libraries(otmt_nontl  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

add_executable(otmt_2op otmt_2op.cpp)
target_link_libraries(otmt_2op  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

#add_executable(Read_MT_nontl Read_MT_nontl.cpp)
#target_link_libraries(Read_MT_nontl  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

#add_executable(clientmt clientmt.cpp)
#target_link_libraries(clientmt  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

#add_executable(linear_clientMT linear_clientMT.cpp)
#target_link_libraries(linear_clientMT  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

#add_executable(logistic_clientMT logistic_clientMT.cpp)
#target_link_libraries(logistic_clientMT  relic  ${OPENSSL_LIBRARIES} m gmp ${EMP-TOOL_LIBRARIES})

