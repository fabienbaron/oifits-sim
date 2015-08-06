cmake_minimum_required(VERSION 2.8)

project (oifits-sim)

set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -ftree-vectorize -march=native -ggdb -fomit-frame-pointer")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "The compiler ${CMAKE_C_COMPILER} has no OpenMP support. Please use a different C++ compiler if you want to use OpenMP capabilities.")
endif()

find_package(GLIB2 REQUIRED)
include_directories(${GLIB2_INCLUDE_DIRS})
find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})
find_package(FFTW REQUIRED)

include_directories(../lib/oifitslib)
include_directories(${FFTW_INCLUDE_DIRS})
include_directories(../lib/textio/src)

file(GLOB_RECURSE CPP_SOURCE . *.cpp)
file(GLOB_RECURSE C_SOURCE . *.c)

add_executable(oifits-sim ${CPP_SOURCE} ${C_SOURCE})

#SET_TARGET_PROPERTIES(oifits-sim PROPERTIES LINKER_LANGUAGE CPP)

target_link_libraries(oifits-sim textio oifits ${CFITSIO_LIBRARIES} ${FFTW_LIBRARIES} ${GLIB2_LIBRARIES})

file(GLOB CONFIG_FILES ${CMAKE_SOURCE_DIR}/etc/*)
file(COPY ${CONFIG_FILES} DESTINATION ${EXECUTABLE_OUTPUT_PATH}/etc/)
