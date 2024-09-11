cmake_minimum_required (VERSION 2.8.11)
project (pailler)

set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} /usr/local")
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories("/usr/include/eigen3")

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} cmake/)
FIND_PACKAGE(emp-tool)

FIND_PACKAGE(OpenSSL)
INCLUDE_DIRECTORIES(${OPENSSL_INCLUDE_DIR})

set (CMAKE_C_FLAGS "-pthread -Wall -march=native -O3 -maes")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -ggdb")
set (CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++0x")

#add_executable(tmp test.cpp)
#target_link_libraries(tmp   ${OPENSSL_LIBRARIES} gmp )

#add_executable(hemt_test hemt_test.cpp)
#target_link_libraries(hemt_test   ${OPENSSL_LIBRARIES} gmp )

add_executable(dgkmt dgkmt.cpp dgk.cpp powmod.cpp)
target_link_libraries(dgkmt   ${OPENSSL_LIBRARIES} gmp )
